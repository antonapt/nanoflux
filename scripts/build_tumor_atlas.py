#!/usr/bin/env python
"""Build a tumour model from read tables: pooled per-site counts + a tumour calibration table.

    python scripts/build_tumor_atlas.py --control-atlas reference.feather \
        --samples s_16/extracted.tsv s_20/extracted.tsv s_21/extracted.tsv --out tumor_atlas/

Writes
    <out>/tumor_atlas.feather        chrom_x, start_genomic (1-based plus-strand C), methylated_pooled, total_pooled
    <out>/tumor_calibration.json     Calibration table fitted on the tumour reads against the control atlas:
                                     per (control p bin, control coverage stratum) the share of methylated tumour calls
    <out>/build_info.json

Use with:  nanoflux score ... --tumor-atlas <out>/tumor_atlas.feather --tumor-calibration <out>/tumor_calibration.json
For a sample that is part of the pool, build a second model without it (leave-one-out).
A call is counted methylated when the read's probability of methylation exceeds 0.5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.scoring.atlas import Atlas
from src.scoring.calibration import Calibration
from src.scoring.calls import FORMATS


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control-atlas", required=True)
    ap.add_argument("--samples", nargs="+", required=True, help="read tables aligned to the atlas build")
    ap.add_argument("--format", choices=list(FORMATS), default="modkit")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    atlas = Atlas.from_feather(args.control_atlas)
    calib = Calibration()
    counts: list[pd.DataFrame] = []
    info = {"samples": [], "control_atlas": str(args.control_atlas)}
    for path in args.samples:
        n_calls = n_matched = 0
        for calls, stats in FORMATS[args.format](path, chunksize=args.chunksize):
            n_calls += len(calls)
            m, n, cnt = atlas.lookup(calls["chrom"], calls["start"], calls["end"])
            found = cnt > 0
            n_matched += int(found.sum())
            calib.add(m[found], n[found], calls["r"].to_numpy()[found])
            meth = (calls["r"].to_numpy() > 0.5).astype(np.int64)
            counts.append(pd.DataFrame({"chrom_x": calls["chrom"].to_numpy(), "start_genomic": calls["start"].to_numpy(),
                                        "methylated_pooled": meth, "total_pooled": 1})
                          .groupby(["chrom_x", "start_genomic"], sort=False).sum().reset_index())
        info["samples"].append({"path": str(path), "n_calls": n_calls, "n_matched_control_atlas": n_matched})
        print(f"{path}: {n_calls:,} calls, {n_matched:,} on control-atlas sites", flush=True)

    pooled = pd.concat(counts).groupby(["chrom_x", "start_genomic"], sort=True).sum().reset_index()
    pooled.to_feather(args.out / "tumor_atlas.feather")
    calib.save(args.out / "tumor_calibration.json")
    cov = pooled["total_pooled"]
    info.update({"n_sites": int(len(pooled)), "coverage_median": float(cov.median()),
                 "share_sites_cov_ge_3": float((cov >= 3).mean()), "calibration_calls": calib.n_calls,
                 "mean_methylated_fraction": float(pooled["methylated_pooled"].sum() / cov.sum())})
    (args.out / "build_info.json").write_text(json.dumps(info, indent=2))
    print(json.dumps({k: v for k, v in info.items() if k != "samples"}, indent=2))


if __name__ == "__main__":
    main()
