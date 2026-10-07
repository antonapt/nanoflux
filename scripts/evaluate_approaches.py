#!/usr/bin/env python
"""Compare ways of using per-read scores before classification, on a cohort with known labels.

    python scripts/evaluate_approaches.py --sheet eval_sheet.tsv --out results/ [--ref chm13v2] [--threads 8]

Sheet (tab-separated, one row per sample):
    sample_id     name
    category      control / tumor_low / tumor_high (only used for grouping the output)
    true_label    class name exactly as the classifier reports it, e.g. "CONTR, INFLAM", "GBM", "MB, G3 G4"
    reads_tsv     modkit read-calls table on the classifier's genome build (reads.tsv of `nanoflux prepare --extract-reads`)
    scores        read_scores.parquet of `nanoflux score`; with --tumor-atlas it holds all five scores
                  (z_control, llr_control_random_per_call, z_tumor, llr_tumor_random_per_call,
                  llr_tumor_control_per_call); old column names z / llr_per_call / llr_tumor_per_call are accepted
    scores_tumor  optional (old layout): a second parquet whose tumour columns are merged in

For every score and sample, reads are ranked by tumour-likeness within the sample (u in (0, 1],
1 = most tumour-like; control-reference scores are flipped so that "unlike control" ranks high).
Using ranks gives every score exactly the same weight distribution, so scores can be compared with
each other: the only difference is WHICH reads they put on top. It also settles the threshold
direction once: with a control reference the methods remove what looks like control, with a tumour
reference they keep what looks like tumour, and both are "keep high u".

Approaches
    baseline          all reads, weight 1
    weighted pileup   per probe: sum(w * methylated) / sum(w),  w = u^k            k in RANK_POWERS
    site weights      probe feature (+1 / -1) scaled by the mean w of its reads, rescaled to mean 1
    read filter       keep only the top q most tumour-like reads, then an unweighted pileup
    probe filter      keep a probe only if max w >= tau (its most tumour-like read) or sum w >= tau, w = u

Calls below the sample's pass threshold (10th percentile of call confidence, like modkit) are dropped.
Reads without a score keep weight 1 in the pileup, are ignored for site weights and dropped by the filters.

Outputs in --out:  approaches_long.csv (one row per score, approach, variant, sample),
                   approaches_summary.csv (correct calls and mean probability of the true label per group),
                   tumor_low_deltas.csv / tumor_low_summary.csv (paired change against baseline for the
                   primary test category), approaches_table.txt (the printed tables)
"""
from __future__ import annotations

import argparse
import io
from contextlib import redirect_stdout
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd

from src.scoring.calls import _NON_PRIMARY
from src.scoring.pileup import FILTER_PERCENTILE, _SiteIndex, load_annotation
from src.scoring.scores import ALL_SCORES, TUMOR_LIKE_IS_HIGH, rename_score_columns

RANK_POWERS = [1, 2, 4, 8]
KEEP_TOP = [0.5, 0.3, 0.1, 0.05, 0.02]
PROBE_MAX_W = [0.5, 0.7, 0.9]
PROBE_SUM_W = [0.6, 1.0, 1.5]
SITE_POWERS = [1, 2, 4]
CONTROL = "CONTR, INFLAM"
PRIMARY_CATEGORY = "tumor_low"


def load_probe_calls(reads_tsv: Path, sites: _SiteIndex, probes: np.ndarray, chunksize: int) -> pd.DataFrame:
    """Passing calls on annotation probes: probe, is_m, read_id."""
    parts, hist = [], np.zeros(1024, np.int64)
    usecols = ["read_id", "ref_position", "chrom", "ref_strand", "call_prob", "call_code", "flag"]
    for chunk in pd.read_csv(reads_tsv, sep="\t", chunksize=chunksize, usecols=usecols,
                             dtype={"read_id": "string", "chrom": "category", "ref_strand": "category",
                                    "call_code": "category", "flag": np.int64}):
        code = chunk["call_code"].astype(str).to_numpy()
        keep = ((chunk["flag"].to_numpy() & _NON_PRIMARY) == 0) & ((code == "m") | (code == "-"))
        sub = chunk.loc[keep]
        prob = sub["call_prob"].to_numpy(np.float64)
        hist += np.bincount(np.clip((prob * 1024).astype(np.int64), 0, 1023), minlength=1024)
        plus = sub["ref_strand"].astype(str).to_numpy() == "+"
        c_pos = sub["ref_position"].to_numpy(np.int64) - (~plus).astype(np.int64)
        site = sites.lookup(sub["chrom"].astype(str).to_numpy(), c_pos)
        on = site >= 0
        parts.append(pd.DataFrame({"probe": probes[site[on]], "is_m": (code[keep][on] == "m").astype(float),
                                   "prob": prob[on], "read_id": sub["read_id"].astype(str).to_numpy()[on]}))
    calls = pd.concat(parts, ignore_index=True)
    threshold = np.searchsorted(np.cumsum(hist) / max(hist.sum(), 1), FILTER_PERCENTILE) / 1024
    return calls[calls["prob"] >= threshold].drop(columns="prob").reset_index(drop=True)


def tumor_rank(scores: pd.DataFrame, column: str) -> pd.Series:
    """Within-sample rank in (0, 1] over all scored reads; 1 = most tumour-like.

    Control-reference scores (z_control, llr_control_random_per_call) are low when a
    read does not look like control, so they are flipped; tumour-reference scores are
    high when a read looks like tumour. After ranking, the direction question is settled
    once for every enrichment method.
    """
    s = scores.set_index("read_id")[column].dropna()
    a = s if column in TUMOR_LIKE_IS_HIGH else -s
    return a.rank(method="average") / len(a)


def load_scores(row) -> pd.DataFrame:
    """Read score table(s) of one sample; merges an optional scores_tumor table (old layout)."""
    table = rename_score_columns(pd.read_parquet(row.scores))
    extra = getattr(row, "scores_tumor", "")
    if extra:
        t2 = rename_score_columns(pd.read_parquet(extra))
        new = [c for c in t2.columns if c in ALL_SCORES and c not in table.columns]
        if new:
            table = table.merge(t2[["read_id"] + new], on="read_id", how="outer")
    return table


def binarise(frac: pd.Series) -> pd.Series:
    return pd.Series(np.where(frac >= 0.7, 1.0, np.where(frac < 0.3, -1.0, np.nan)), index=frac.index)


class Classifier:
    def __init__(self, threads: int) -> None:
        from src.model.ensemble_model import EnsembleWrapper
        from src.tasks.infer import load_features_list
        from src.utils.data_utils import normalize_X

        self.features = load_features_list()
        self.normalize = normalize_X
        self.ensemble = EnsembleWrapper(models_dir=files("data") / "models", n_threads=threads)

    def predict(self, x: pd.Series, name: str) -> tuple[pd.Series, int]:
        xp = pd.Series(0.0, index=self.features, name=name)
        common = x.index.intersection(xp.index)
        xp.loc[common] = x.loc[common].astype(np.float32)
        xp[xp.isna()] = 0
        probas, _ = self.ensemble.predict(self.normalize(xp.to_frame().T))
        return probas.iloc[0].sort_values(ascending=False), int((xp != 0).sum())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--ref", default="chm13v2", help="annotation build of the classifier input (default chm13v2)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    sheet = pd.read_csv(args.sheet, sep="\t", dtype=str).fillna("")
    for col in ("sample_id", "category", "true_label", "reads_tsv", "scores"):
        if col not in sheet.columns:
            raise SystemExit(f"sheet needs a '{col}' column; found {list(sheet.columns)}")
    if "scores_tumor" not in sheet.columns:
        sheet["scores_tumor"] = ""

    anno = load_annotation(Path(str(files("data") / "features" / f"EPIC_{args.ref.strip('-as')}.bed")))
    sites, probes = _SiteIndex(anno), anno["probe"].to_numpy()
    clf = Classifier(args.threads)
    rows: list[dict] = []

    def record(x, s, score, approach, variant, extra=None):
        p, n_measured = clf.predict(x, s.sample_id)
        rows.append({"score": score, "approach": approach, "variant": variant, "sample": s.sample_id,
                     "category": s.category, "true_label": s.true_label, "top_class": p.index[0],
                     "top_prob": round(float(p.iloc[0]), 3), "correct": p.index[0] == s.true_label,
                     "p_true_label": round(float(p.get(s.true_label, np.nan)), 3),
                     "p_control": round(float(p.get(CONTROL, np.nan)), 3), "n_measured": n_measured,
                     **(extra or {})})

    for s in sheet.itertuples(index=False):
        calls = load_probe_calls(Path(s.reads_tsv), sites, probes, args.chunksize)
        state = binarise(calls.groupby("probe")["is_m"].mean())
        n_probes = int(state.notna().sum())
        record(state, s, "-", "baseline", "all reads, weight 1")
        table = load_scores(s)
        for score in ALL_SCORES:
            if score not in table.columns or table[score].notna().sum() == 0:
                continue
            u = tumor_rank(table, score).reindex(calls["read_id"]).to_numpy()
            scored = ~np.isnan(u)
            u0 = np.nan_to_num(u)
            # weighted pileup and site weights: w = u^k
            for k in sorted(set(RANK_POWERS) | set(SITE_POWERS)):
                w = np.where(scored, u0 ** k, 1.0)
                c = calls.assign(w=w, wm=w * calls["is_m"].to_numpy())
                g = c.groupby("probe")
                if k in RANK_POWERS:
                    record(binarise(g["wm"].sum() / g["w"].sum()), s, score, "weighted pileup", f"w = u^{k}")
                if k in SITE_POWERS:
                    omega = c[scored].groupby("probe")["w"].mean().reindex(state.index)
                    omega = (omega / omega[state.notna()].mean()).fillna(1.0)
                    record(state * omega, s, score, "site weights", f"w = u^{k}")
            # read filter: keep the top q most tumour-like reads, unweighted pileup
            for q in KEEP_TOP:
                keep = scored & (u0 > 1 - q)
                record(binarise(calls[keep].groupby("probe")["is_m"].mean()), s, score, "read filter",
                       f"keep top {q:.0%}", {"calls_kept": round(float(keep.sum() / max(scored.sum(), 1)), 3)})
            # probe filter: keep a probe if its most tumour-like read (max w) or its summed w is high enough
            per_probe = calls[scored].assign(w=u0[scored]).groupby("probe")["w"].agg(["max", "sum"])
            for tau in PROBE_MAX_W:
                keep_probes = per_probe.index[per_probe["max"] >= tau]
                record(state.reindex(keep_probes).dropna(), s, score, "probe filter", f"max w >= {tau}",
                       {"probes_kept": round(float(len(keep_probes) / max(n_probes, 1)), 3)})
            for tau in PROBE_SUM_W:
                keep_probes = per_probe.index[per_probe["sum"] >= tau]
                record(state.reindex(keep_probes).dropna(), s, score, "probe filter", f"sum w >= {tau}",
                       {"probes_kept": round(float(len(keep_probes) / max(n_probes, 1)), 3)})
        print(f"[done] {s.sample_id}", flush=True)

    tab = pd.DataFrame(rows)
    tab.to_csv(args.out / "approaches_long.csv", index=False)

    tab["row"] = tab["score"] + " | " + tab["approach"] + " | " + tab["variant"]
    row_order = list(dict.fromkeys(tab["row"]))
    summary = tab.groupby("row").agg(n_correct=("correct", "sum"), n_samples=("correct", "size")).reindex(row_order)
    by_cat = tab.pivot_table(index="row", columns="category", values="p_true_label", aggfunc="mean").reindex(row_order)
    by_cat.columns = [f"mean p(true label), {c}" for c in by_cat.columns]
    acc_cat = tab.pivot_table(index="row", columns="category", values="correct", aggfunc="sum").reindex(row_order)
    acc_cat.columns = [f"correct, {c}" for c in acc_cat.columns]
    summary = pd.concat([summary, acc_cat, by_cat.round(3)], axis=1)
    summary.to_csv(args.out / "approaches_summary.csv")

    # primary read-out: per-sample change against the baseline, paired, for the category of interest
    base = tab[tab["approach"] == "baseline"].set_index("sample")
    prim = tab[tab["category"] == PRIMARY_CATEGORY].copy()
    prim["d_p_true_label"] = prim["p_true_label"].to_numpy() - base["p_true_label"].reindex(prim["sample"]).to_numpy()
    prim["d_p_control"] = prim["p_control"].to_numpy() - base["p_control"].reindex(prim["sample"]).to_numpy()
    prim["d_n_measured"] = prim["n_measured"].to_numpy() - base["n_measured"].reindex(prim["sample"]).to_numpy()
    prim[["score", "approach", "variant", "sample", "true_label", "top_class", "correct", "p_true_label",
          "d_p_true_label", "p_control", "d_p_control", "n_measured", "d_n_measured"]].to_csv(
        args.out / f"{PRIMARY_CATEGORY}_deltas.csv", index=False)
    prim["row"] = prim["score"] + " | " + prim["approach"] + " | " + prim["variant"]
    deltas = prim[prim["approach"] != "baseline"].groupby("row").agg(
        n=("d_p_true_label", "size"),
        n_improved=("d_p_true_label", lambda d: int((d > 0).sum())),
        n_worse=("d_p_true_label", lambda d: int((d < 0).sum())),
        mean_d_p_true_label=("d_p_true_label", "mean"),
        median_d_p_true_label=("d_p_true_label", "median"),
        mean_d_p_control=("d_p_control", "mean"),
        n_correct=("correct", "sum"),
        min_n_measured=("n_measured", "min"),
    ).reindex([r for r in row_order if r in set(prim.loc[prim["approach"] != "baseline", "row"])]).round(3)
    deltas.to_csv(args.out / f"{PRIMARY_CATEGORY}_summary.csv")

    buf = io.StringIO()
    with redirect_stdout(buf):
        order = sheet["sample_id"].tolist()
        tab["cell"] = tab.apply(lambda x: f"{x.top_class} {x.top_prob:.2f}" + (" *" if x.correct else ""), axis=1)
        print("true labels:")
        print(sheet.set_index("sample_id")[["category", "true_label"]].T.to_string())
        print("\nsummary (correct calls, mean probability of the true label per group):")
        print(summary.to_string())
        if len(deltas):
            print(f"\n{PRIMARY_CATEGORY}: change against baseline per variant (paired over samples; "
                  "n_improved / n_worse = samples whose p(true label) rose / fell):")
            print(deltas.to_string())
        print("\ntop class and probability (* = correct):")
        print(tab.pivot(index="row", columns="sample", values="cell").reindex(row_order)[order].to_string())
        print("\nprobability of the true label:")
        print(tab.pivot(index="row", columns="sample", values="p_true_label").reindex(row_order)[order].to_string())
    (args.out / "approaches_table.txt").write_text(buf.getvalue())
    print(buf.getvalue())


if __name__ == "__main__":
    main()
