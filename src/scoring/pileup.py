"""Weighted CpG pileup from a modkit extract read table.

Replaces ``modkit pileup`` + ``bedtools intersect`` when reads carry weights.
Every call is counted with the weight of its read, so a read with weight 0.2
contributes 0.2 to the coverage and, if methylated, 0.2 to the methylated
count. With all weights equal to 1 the percentages equal modkit's.

Mirrors the ``traditional`` pileup preset: CpG context, strands combined at
the plus-strand C, calls below the pass threshold counted as failed. The
threshold is either given or, like modkit's default, the 10th percentile of
the call confidences. Only sites of the annotation are kept, and the output
is written in the layout ``bedtools intersect -wa -wb`` would produce from a
bedMethyl and the annotation, which is what ``nanoflux infer`` reads:
columns 0-2 chromosome/start/end, column 9 the space-separated
``Nvalid percent Nmod Ncanon Nother Ndelete Nfail Ndiff Nnocall`` block and
column 13 the probe id.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.scoring.calls import _NON_PRIMARY
from src.utils.log import logger

FILTER_PERCENTILE = 0.10  # modkit default --filter-percentile
_PROB_BINS = 1024


def load_annotation(path: str | Path) -> pd.DataFrame:
    anno = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3],
                       names=["chrom", "start", "end", "probe"],
                       dtype={"chrom": str, "start": np.int64, "end": np.int64, "probe": str})
    return anno.reset_index(drop=True)


def load_weights(path: str | Path, *, column: str = "weight") -> pd.Series:
    df = pd.read_parquet(path, columns=["read_id", column])
    if df["read_id"].duplicated().any():
        raise ValueError(f"duplicate read ids in {path}")
    return pd.Series(df[column].to_numpy(np.float64), index=pd.Index(df["read_id"].astype(str)))


class _SiteIndex:
    """Annotation lookup: plus-strand CpG [c, c+2) overlapping a probe interval.

    Probes are short intervals (2 bp for EPIC), so a CpG matches a probe when
    ``c < probe.end`` and ``c + 2 > probe.start``. Keyed per chromosome by the
    C positions the probe can start from, which for 2 bp probes is c-1, c, c+1.
    """

    def __init__(self, anno: pd.DataFrame) -> None:
        self.chrom_codes = {c: i for i, c in enumerate(pd.unique(anno["chrom"]))}
        codes = anno["chrom"].map(self.chrom_codes).to_numpy(np.int64)
        keys, rows = [], []
        for row, (code, s, e) in enumerate(zip(codes, anno["start"].to_numpy(), anno["end"].to_numpy())):
            for c in range(max(s - 1, 0), e):  # every C position whose [c, c+2) overlaps [s, e)
                keys.append((code << 32) | c)
                rows.append(row)
        keys, rows = np.asarray(keys, np.int64), np.asarray(rows, np.int64)
        order = np.argsort(keys, kind="stable")
        keys, rows = keys[order], rows[order]
        first = np.concatenate([[True], keys[1:] != keys[:-1]])  # one probe per C (keep first)
        self.index = pd.Index(keys[first])
        self.rows = rows[first]

    def lookup(self, chrom: np.ndarray, c_pos: np.ndarray) -> np.ndarray:
        cat = pd.Categorical(chrom)
        per_cat = np.asarray([self.chrom_codes.get(str(c), -1) for c in cat.categories], np.int64)
        codes = per_cat[cat.codes] if len(per_cat) else np.full(len(cat), -1, np.int64)
        known = codes >= 0
        key = (np.where(known, codes, 0) << 32) | np.asarray(c_pos, np.int64)
        loc = self.index.get_indexer(key)
        out = np.full(len(key), -1, np.int64)
        ok = known & (loc >= 0)
        out[ok] = self.rows[loc[ok]]
        return out


def weighted_pileup(
    reads_tsv: str | Path,
    annotation: str | Path,
    weights: pd.Series,
    output: str | Path,
    *,
    filter_threshold: float | None = None,
    unscored_weight: float = 1.0,
    chunksize: int = 2_000_000,
) -> dict:
    """Write the annotated weighted pileup to ``output``; return run statistics."""
    anno = load_annotation(annotation)
    sites = _SiteIndex(anno)
    usecols = ["read_id", "ref_position", "chrom", "ref_strand", "call_prob", "call_code", "flag"]

    # pass over the read table: keep calls on annotation sites, histogram all confidences
    kept, hist = [], np.zeros(_PROB_BINS, np.int64)
    stats = {"n_rows": 0, "n_non_primary": 0, "n_other_mod": 0, "n_calls": 0,
             "n_calls_on_probes": 0, "n_reads_seen": 0, "n_reads_unscored": 0}
    seen, unscored = set(), set()
    reader = pd.read_csv(reads_tsv, sep="\t", usecols=usecols, chunksize=chunksize,
                         dtype={"read_id": "string", "ref_position": np.int64, "chrom": "category",
                                "ref_strand": "category", "call_prob": np.float64,
                                "call_code": "category", "flag": np.int64})
    for chunk in reader:
        stats["n_rows"] += len(chunk)
        primary = (chunk["flag"].to_numpy() & _NON_PRIMARY) == 0
        code = chunk["call_code"].astype(str).to_numpy()
        is_m, is_can = code == "m", code == "-"
        keep = primary & (is_m | is_can)
        stats["n_non_primary"] += int((~primary).sum())
        stats["n_other_mod"] += int((primary & ~keep).sum())
        stats["n_calls"] += int(keep.sum())
        sub = chunk.loc[keep]
        prob = sub["call_prob"].to_numpy(np.float64)
        hist += np.bincount(np.clip((prob * _PROB_BINS).astype(np.int64), 0, _PROB_BINS - 1),
                            minlength=_PROB_BINS)

        # plus-strand read: ref_position is the C; minus-strand read: it is the G, the C is one left
        plus = sub["ref_strand"].astype(str).to_numpy() == "+"
        c_pos = sub["ref_position"].to_numpy(np.int64) - (~plus).astype(np.int64)
        site = sites.lookup(sub["chrom"].astype(str).to_numpy(), c_pos)
        on = site >= 0
        if not on.any():
            continue
        rid = sub["read_id"].astype(str).to_numpy()[on]
        seen.update(np.unique(rid).tolist())
        w = weights.reindex(rid).to_numpy(np.float64)
        missing = np.isnan(w)
        unscored.update(np.unique(rid[missing]).tolist())
        w[missing] = unscored_weight
        kept.append(pd.DataFrame({"site": site[on], "w": w, "prob": prob[on],
                                  "is_m": is_m[keep][on]}))
    stats["n_reads_seen"], stats["n_reads_unscored"] = len(seen), len(unscored)
    if not kept:
        raise ValueError(f"no calls of {reads_tsv} fall on annotation sites")
    calls = pd.concat(kept, ignore_index=True)
    stats["n_calls_on_probes"] = int(len(calls))

    # pass threshold: given, or the modkit default percentile of all call confidences
    if filter_threshold is None:
        cdf = np.cumsum(hist) / max(hist.sum(), 1)
        filter_threshold = float(np.searchsorted(cdf, FILTER_PERCENTILE) / _PROB_BINS)
        stats["filter_threshold_source"] = f"{FILTER_PERCENTILE:.0%} percentile of call confidence"
    else:
        stats["filter_threshold_source"] = "given"
    stats["filter_threshold"] = filter_threshold

    passed = calls["prob"].to_numpy() >= filter_threshold
    w = calls["w"].to_numpy()
    site = calls["site"].to_numpy()
    n_sites = len(anno)
    valid = np.bincount(site, weights=w * passed, minlength=n_sites)
    mod = np.bincount(site, weights=w * (passed & calls["is_m"].to_numpy()), minlength=n_sites)
    fail = np.bincount(site, weights=w * ~passed, minlength=n_sites)
    covered = valid > 0
    stats["n_probes_covered"] = int(covered.sum())
    stats["n_probes_total"] = int(n_sites)

    a = anno.loc[covered]
    v, mo, f = valid[covered], mod[covered], fail[covered]
    block = [f"{vi:.4g} {100 * mi / vi:.2f} {mi:.4g} {vi - mi:.4g} 0 0 {fi:.4g} 0 0"
             for vi, mi, fi in zip(v, mo, f)]
    out = pd.DataFrame({
        0: a["chrom"].to_numpy(), 1: a["start"].to_numpy(), 2: a["start"].to_numpy() + 1,
        3: "m", 4: np.round(v).astype(np.int64), 5: ".", 6: a["start"].to_numpy(),
        7: a["start"].to_numpy() + 1, 8: "255,0,0", 9: block,
        10: a["chrom"].to_numpy(), 11: a["start"].to_numpy(), 12: a["end"].to_numpy(),
        13: a["probe"].to_numpy(),
    })
    out.to_csv(output, sep="\t", header=False, index=False)
    logger.info(
        f"Weighted pileup: {stats['n_calls_on_probes']:,} calls on {stats['n_probes_covered']:,} "
        f"probes from {stats['n_reads_seen']:,} reads ({stats['n_reads_unscored']:,} without a "
        f"score, weight {unscored_weight}); pass threshold {filter_threshold:.3f} "
        f"({stats['filter_threshold_source']})"
    )
    return stats
