#!/usr/bin/env python
"""Build a tumour model from read tables: pooled per-site counts + a tumour calibration table.

    python scripts/build_tumor_atlas.py --control-atlas reference.feather \
        --samples s_16/extracted.tsv s_20/extracted.tsv s_21/extracted.tsv --out tumor_atlas/

Writes
    <out>/tumor_atlas.feather        chrom_x, start_genomic (1-based plus-strand C), methylated_pooled, total_pooled
    <out>/tumor_prior.json           Beta prior (alpha, beta) fitted on the well-covered tumour atlas sites, like
                                     the control-side prior; with the sites used and whether --prior-min-cov held
    <out>/tumor_calibration.json     Calibration table fitted on the tumour reads against the control atlas:
                                     per (control p bin, control coverage stratum) the share of methylated tumour calls
    <out>/build_info.json            sample call counts, coverage report (median, share of sites >= 3/10/20) and the prior

Use with (deep atlas, Beta prior):   nanoflux score ... --tumor-atlas <out>/tumor_atlas.feather --tumor-prior auto
Use with (sparse atlas, calibration): nanoflux score ... --tumor-atlas <out>/tumor_atlas.feather --tumor-calibration <out>/tumor_calibration.json
For a sample that is part of the pool, build a second model without it (leave-one-out), unless the biased
in-pool evaluation is intended.
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
from src.scoring.prior import fit_beta_prior


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control-atlas", required=True)
    ap.add_argument("--samples", nargs="+", required=True, help="read tables aligned to the atlas build")
    ap.add_argument("--format", choices=list(FORMATS), default="modkit")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    ap.add_argument("--prior-min-cov", type=int, default=20,
                    help="pooled tumour coverage a site needs to enter the Beta prior fit (default 20)")
    ap.add_argument("--prior-min-sites", type=int, default=1000,
                    help="sites required at --prior-min-cov; fewer is an error unless --allow-shallow-prior (default 1000)")
    ap.add_argument("--allow-shallow-prior", action="store_true",
                    help="only warn when fewer than --prior-min-sites sites reach --prior-min-cov")
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
    cov = pooled["total_pooled"].to_numpy()
    meth = pooled["methylated_pooled"].to_numpy()
    info.update({"n_sites": int(len(pooled)), "coverage_median": float(np.median(cov)),
                 "share_sites_cov_ge_3": float((cov >= 3).mean()),
                 "share_sites_cov_ge_10": float((cov >= 10).mean()),
                 "share_sites_cov_ge_20": float((cov >= 20).mean()),
                 "calibration_calls": calib.n_calls,
                 "mean_methylated_fraction": float(meth.sum() / cov.sum())})

    # Beta prior on the tumour atlas, mirroring the control-side prior (nanoflux score --prior auto)
    n_hi = int((cov >= args.prior_min_cov).sum())
    prior = fit_beta_prior(meth, cov, min_cov=args.prior_min_cov, min_sites=args.prior_min_sites)
    prior_info = {"alpha": prior.alpha, "beta": prior.beta, "mean": prior.mean,
                  "requested_min_cov": args.prior_min_cov, "n_sites_at_requested_min_cov": n_hi,
                  "fit_ok": n_hi >= args.prior_min_sites,
                  "note": "use with: nanoflux score --tumor-atlas tumor_atlas.feather --tumor-prior ALPHA BETA "
                          "(or --tumor-prior auto --tumor-prior-min-cov N to refit)"}
    (args.out / "tumor_prior.json").write_text(json.dumps(prior_info, indent=2))
    info["prior"] = prior_info
    (args.out / "build_info.json").write_text(json.dumps(info, indent=2))
    print(json.dumps({k: v for k, v in info.items() if k != "samples"}, indent=2))
    if not prior_info["fit_ok"]:
        msg = (f"only {n_hi:,} sites have coverage >= {args.prior_min_cov} (need {args.prior_min_sites:,}); "
               f"the prior was fitted at a lower coverage. Lower --prior-min-cov deliberately or add samples.")
        if args.allow_shallow_prior:
            print(f"[WARN] {msg}")
        else:
            raise SystemExit(f"[ERROR] {msg} Atlas and calibration were written; pass --allow-shallow-prior to accept.")


if __name__ == "__main__":
    main()
