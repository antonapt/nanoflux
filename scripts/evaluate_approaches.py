#!/usr/bin/env python
"""Compare ways of using per-read scores before classification, on a cohort with known labels.

    python scripts/evaluate_approaches.py --sheet eval_sheet.tsv --out results/ [--ref chm13v2] [--threads 8]

Sheet (tab-separated, one row per sample):
    sample_id     name
    category      control / tumor_low / tumor_high (only used for grouping the output)
    true_label    class name exactly as the classifier reports it, e.g. "CONTR, INFLAM", "GBM", "MB, G3 G4"
    reads_tsv     modkit read-calls table on the classifier's genome build (reads.tsv of `nanoflux prepare --extract-reads`)
    scores        read_scores.parquet of `nanoflux score` (z, llr_per_call)
    scores_tumor  optional: read_scores.parquet of a run with --tumor-calibration (llr_tumor_per_call)

For every score and sample, reads are ranked by abnormality within the sample (u in (0, 1], 1 = least
like the control atlas). Using ranks gives every score exactly the same weight distribution, so scores
can be compared with each other: the only difference is WHICH reads they put on top.

Approaches
    baseline          all reads, weight 1
    weighted pileup   per probe: sum(w * methylated) / sum(w),  w = u^k
    site weights      probe feature (+1 / -1) scaled by the mean w of its reads, rescaled to mean 1
    read filter       keep only the top q most abnormal reads, then an unweighted pileup

Calls below the sample's pass threshold (10th percentile of call confidence, like modkit) are dropped.
Reads without a score keep weight 1 in the pileup, are ignored for site weights and dropped by the filter.

Outputs in --out:  approaches_long.csv (one row per score, approach, variant, sample),
                   approaches_summary.csv (correct calls and mean probability of the true label per group),
                   approaches_table.txt (the printed tables)
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

HIGHER_IS_ABNORMAL = {"llr_tumor_per_call"}
SCORE_SOURCE = {"z": "scores", "llr_per_call": "scores", "llr_tumor_per_call": "scores_tumor"}
RANK_POWERS = [1, 2, 4]
KEEP_TOP = [0.5, 0.3, 0.1]
CONTROL = "CONTR, INFLAM"


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


def abnormality_rank(scores: pd.DataFrame, column: str) -> pd.Series:
    """Within-sample rank in (0, 1] over all scored reads; 1 = least like the control atlas."""
    s = scores.set_index("read_id")[column].dropna()
    a = s if column in HIGHER_IS_ABNORMAL else -s
    return a.rank(method="average") / len(a)


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
        record(state, s, "-", "baseline", "all reads, weight 1")
        for score, source in SCORE_SOURCE.items():
            path = getattr(s, source)
            if not path:
                continue
            table = pd.read_parquet(path)
            if score not in table.columns:
                continue
            u = abnormality_rank(table, score).reindex(calls["read_id"]).to_numpy()
            scored = ~np.isnan(u)
            for k in RANK_POWERS:
                w = np.where(scored, np.nan_to_num(u) ** k, 1.0)
                c = calls.assign(w=w, wm=w * calls["is_m"].to_numpy())
                g = c.groupby("probe")
                record(binarise(g["wm"].sum() / g["w"].sum()), s, score, "weighted pileup", f"w = u^{k}")
                omega = c[scored].groupby("probe")["w"].mean().reindex(state.index)
                omega = (omega / omega[state.notna()].mean()).fillna(1.0)
                record(state * omega, s, score, "site weights", f"w = u^{k}")
            for q in KEEP_TOP:
                keep = scored & (u > 1 - q)
                record(binarise(calls[keep].groupby("probe")["is_m"].mean()), s, score, "read filter",
                       f"keep top {q:.0%}", {"calls_kept": round(float(keep.sum() / max(scored.sum(), 1)), 3)})
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

    buf = io.StringIO()
    with redirect_stdout(buf):
        order = sheet["sample_id"].tolist()
        tab["cell"] = tab.apply(lambda x: f"{x.top_class} {x.top_prob:.2f}" + (" *" if x.correct else ""), axis=1)
        print("true labels:")
        print(sheet.set_index("sample_id")[["category", "true_label"]].T.to_string())
        print("\nsummary (correct calls, mean probability of the true label per group):")
        print(summary.to_string())
        print("\ntop class and probability (* = correct):")
        print(tab.pivot(index="row", columns="sample", values="cell").reindex(row_order)[order].to_string())
        print("\nprobability of the true label:")
        print(tab.pivot(index="row", columns="sample", values="p_true_label").reindex(row_order)[order].to_string())
    (args.out / "approaches_table.txt").write_text(buf.getvalue())
    print(buf.getvalue())


if __name__ == "__main__":
    main()
