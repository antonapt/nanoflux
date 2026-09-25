"""``nanoflux score``: score every read against a pooled CpG atlas and derive a weight."""

from __future__ import annotations

import json
import logging
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from src.scoring.atlas import Atlas
from src.scoring.calibration import Calibration
from src.scoring.calls import FORMATS
from src.scoring.prior import BetaPrior, MomentStats, PriorModel, fit_beta_prior
from src.scoring.scores import finalize_reads, score_calls
from src.scoring.weights import N_BIN_LABELS, bin_quantiles, compute_weights, parse_thresholds
from src.utils.filehandling import prepare_location
from src.utils.log import logger

OUTPUT_COLUMNS = ["read_id", "n_calls", "n_cpg", "n_bin", "z", "llr_per_call", "weight"]


def parse_prior(values: list[str]) -> BetaPrior | None:
    """``["auto"]`` -> None (fit from the atlas); ``["a", "b"]`` -> BetaPrior."""
    if len(values) == 1 and values[0].lower() == "auto":
        return None
    if len(values) == 2:
        try:
            alpha, beta = float(values[0]), float(values[1])
        except ValueError:
            raise ValueError(f"--prior expects 'auto' or two numbers, got {values}")
        if alpha <= 0 or beta <= 0:
            raise ValueError("--prior alpha and beta must be positive")
        return BetaPrior(alpha, beta)
    raise ValueError(f"--prior expects 'auto' or two numbers, got {values}")


def iter_matched_calls(input_file: Path, atlas: Atlas, *, fmt: str, chunksize: int, tau: float,
                       totals: Counter, reads_seen: set):
    """Yield ``(calls, m, n, count)`` for atlas-matched calls, chunk by chunk."""
    for calls, stats in FORMATS[fmt](input_file, chunksize=chunksize, tau=tau):
        reads_seen.update(calls["read_id"].unique().tolist())
        m, n, count = atlas.lookup(calls["chrom"], calls["start"], calls["end"])
        found = count > 0
        stats["n_matched"] = int(found.sum())
        totals.update(stats)
        yield calls.loc[found], m[found], n[found], count[found]


def build_site_model(args, atlas: Atlas, output_dir: Path, iterate) -> tuple[object, dict]:
    """Choose the site model: calibration table (given or fitted) or Beta prior + moments."""
    info: dict = {}
    if args.calibration:
        model = Calibration.load(args.calibration)
        info = {"site_model": "calibration table", "calibration": str(args.calibration),
                "n_calls": model.n_calls}
        logger.info(f"Loaded calibration table from {args.calibration} ({model.n_calls:,} calls)")
        return model, info

    if args.fit_calibration:
        logger.warning(
            "Fitting the calibration table on this sample. Only meaningful if these reads are "
            "NOT part of the atlas and mostly atlas-like."
        )
        model = Calibration()
        for calls, m, n, _ in iterate():
            model.add(m, n, calls["r"].to_numpy())
        model.save(output_dir / "calibration.json")
        info = {"site_model": "calibration table", "calibration": "fitted on sample",
                "n_calls": model.n_calls}
        logger.info(f"Calibration table fitted on {model.n_calls:,} calls")
        return model, info

    prior = parse_prior(args.prior)
    if prior is None:
        prior = fit_beta_prior(atlas.methylated, atlas.total, min_cov=args.prior_min_cov)
        prior_source = f"fitted on atlas sites with coverage >= {args.prior_min_cov}"
    else:
        prior_source = "given"
        logger.info(f"Using Beta prior alpha={prior.alpha}, beta={prior.beta}")

    if args.moments:
        moments = PriorModel.load(args.moments).moments
        moments_source = str(args.moments)
        logger.info(f"Loaded caller moments from {args.moments}")
    else:
        stats = MomentStats()
        for calls, m, n, _ in iterate():
            stats.add(prior.p(m, n), calls["r"].to_numpy())
        moments = stats.fit()
        moments_source = "fitted on sample"
        logger.info(
            f"Caller moments fitted on {moments.n_calls:,} calls: mu1={moments.mu1:.3f} "
            f"v1={moments.v1:.4f} mu0={moments.mu0:.3f} v0={moments.v0:.4f}"
        )
    model = PriorModel(prior, moments)
    model.save(output_dir / "moments.json")
    info = {"site_model": "beta prior + moments", "prior": {"alpha": prior.alpha, "beta": prior.beta,
            "source": prior_source}, "moments": {**model.to_dict()["moments"], "source": moments_source}}
    return model, info


def summarize_bins(reads: pd.DataFrame, cutoffs: dict, medians: dict) -> list[dict]:
    rows = []
    for b in N_BIN_LABELS:
        sel = reads["n_bin"].astype(str) == b
        z = reads.loc[sel, "z"].dropna()
        w = reads.loc[sel, "weight"].dropna()
        rows.append({
            "n_bin": b,
            "n_reads": int(sel.sum()),
            "cutoff": cutoffs[b],
            "median_z": medians[b],
            "z_q01": float(z.quantile(0.01)) if len(z) else np.nan,
            "z_q05": float(z.quantile(0.05)) if len(z) else np.nan,
            "z_sd": float(z.std()) if len(z) else np.nan,
            "mean_weight": float(w.mean()) if len(w) else np.nan,
            "share_weight_above_0.5": float((w > 0.5).mean()) if len(w) else np.nan,
        })
    return rows


def main(args):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    input_file = Path(args.input).resolve()
    output_dir = Path(args.output).resolve()
    out_scores = output_dir / "read_scores.parquet"
    out_summary = output_dir / "score_summary.json"
    for p in (out_scores, out_summary):
        prepare_location(p, args.create_dir)

    logger.info("Running [blue bold]nanoflux score[/]")
    atlas = Atlas.from_feather(args.atlas, columns=tuple(args.atlas_columns), min_cov=args.min_atlas_cov)
    fmt, chunksize, tau = args.format, args.chunksize, args.tau

    def iterate(totals: Counter | None = None, seen: set | None = None):
        return iter_matched_calls(input_file, atlas, fmt=fmt, chunksize=chunksize, tau=tau,
                                  totals=Counter() if totals is None else totals,
                                  reads_seen=set() if seen is None else seen)

    model, model_info = build_site_model(args, atlas, output_dir, iterate)

    # ---- scoring ----
    partials, totals, seen = [], Counter(), set()
    for calls, m, n, count in iterate(totals, seen):
        partials.append(score_calls(calls, m, n, count, model))
    reads = finalize_reads(partials)
    logger.info(
        f"Read {totals['n_rows']:,} rows: {totals['n_non_primary']:,} non-primary and "
        f"{totals['n_other_mod']:,} other-modification rows dropped, "
        f"{totals['n_matched']:,} of {totals['n_calls']:,} calls matched the atlas"
    )
    logger.info(f"Scored {len(reads):,} of {len(seen):,} reads (the rest have no atlas CpG)")

    # ---- weights ----
    z, n_bin = reads["z"].to_numpy(), reads["n_bin"]
    medians = bin_quantiles(z, n_bin, 0.5)
    if args.weight_thresholds:
        cutoffs = parse_thresholds(args.weight_thresholds)
        cutoff_source = "explicit"
    else:
        cutoffs = bin_quantiles(z, n_bin, args.weight_quantile)
        cutoff_source = f"sample quantile {args.weight_quantile}"
    reads["weight"] = compute_weights(z, n_bin, cutoffs, medians, mode=args.weight_mode,
                                      temperature=args.weight_temperature, w_min=args.weight_min)

    out = reads.reset_index()
    out["n_bin"] = out["n_bin"].astype(str)
    pq.write_table(pa.Table.from_pandas(out[OUTPUT_COLUMNS], preserve_index=False),
                   out_scores, compression="zstd")

    bins = summarize_bins(reads, cutoffs, medians)
    summary = {
        "input": str(input_file), "format": fmt, "atlas": str(Path(args.atlas).resolve()),
        "atlas_info": atlas.describe(), "tau": tau, **model_info,
        "weight": {"mode": args.weight_mode, "cutoff_source": cutoff_source,
                   "temperature": args.weight_temperature, "w_min": args.weight_min},
        "counts": dict(totals), "n_reads_in_file": len(seen), "n_reads_scored": int(len(reads)),
        "bins": bins,
    }
    out_summary.write_text(json.dumps(summary, indent=2, default=str))

    logger.info("Per CpG-count bin: reads, cutoff, median z, 1% z, z SD, mean weight")
    for row in bins:
        logger.info(f"  {row['n_bin']:>4}: {row['n_reads']:>9,}  cutoff {row['cutoff']:6.2f}  "
                    f"median {row['median_z']:5.2f}  q01 {row['z_q01']:6.2f}  sd {row['z_sd']:4.2f}  "
                    f"weight {row['mean_weight']:.3f}")
    logger.info("[blue bold]nanoflux score[/] has successfully run!")
    logger.info(f"Per-read scores saved to: {out_scores}")
