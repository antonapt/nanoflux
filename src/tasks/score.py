"""``nanoflux score``: score every read against a pooled CpG atlas and derive a weight."""

from __future__ import annotations

import json
import logging
from collections import Counter
from importlib.resources import files
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
from src.scoring.weights import N_BIN_LABELS, compute_weights, parse_thresholds
from src.utils.filehandling import prepare_location
from src.utils.log import logger

OUTPUT_COLUMNS = ["read_id", "n_calls", "n_cpg", "n_bin", "z", "llr_per_call", "weight"]
TUMOR_COLUMNS = ["llr_tumor", "llr_tumor_per_call"]
DEFAULT_MOMENTS = "csf_controls.json"
# per score: (center, temperature) of the sigmoid weight when not given
WEIGHT_DEFAULTS = {"z": (0.0, 1.0), "llr_per_call": (0.0, 0.25), "llr_tumor_per_call": (0.0, 0.1)}
# scores where a HIGHER value means less like the control atlas (the weight rises with the score)
HIGHER_IS_ABNORMAL = {"llr_tumor_per_call"}


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


def default_moments_path() -> Path:
    return Path(str(files("data") / "moments" / DEFAULT_MOMENTS))


def iter_matched_calls(path: Path, atlas: Atlas, *, fmt: str, chunksize: int, tau: float,
                       totals: Counter | None = None, reads_seen: set | None = None):
    """Yield ``(calls, m, n, count)`` for atlas-matched calls, chunk by chunk."""
    for calls, stats in FORMATS[fmt](path, chunksize=chunksize, tau=tau):
        if reads_seen is not None:
            reads_seen.update(calls["read_id"].unique().tolist())
        m, n, count = atlas.lookup(calls["chrom"], calls["start"], calls["end"])
        found = count > 0
        stats["n_matched"] = int(found.sum())
        if totals is not None:
            totals.update(stats)
        yield calls.loc[found], m[found], n[found], count[found]


def fit_moments(paths: list[Path], atlas: Atlas, prior: BetaPrior, *, min_cov: int, fmt: str,
                chunksize: int, tau: float):
    """Caller moments pooled over the given read tables, at atlas sites with coverage >= min_cov."""
    stats = MomentStats()
    for path in paths:
        for calls, m, n, _ in iter_matched_calls(path, atlas, fmt=fmt, chunksize=chunksize, tau=tau):
            sel = n >= min_cov
            stats.add(prior.p(m[sel], n[sel]), calls["r"].to_numpy()[sel])
        logger.info(f"Moment statistics accumulated from {path.name}: {stats.n:,} calls so far")
    return stats.fit()


def build_site_model(args, atlas: Atlas, input_file: Path, output_dir: Path, *, fmt: str,
                     chunksize: int, tau: float) -> tuple[object, dict]:
    """Choose the site model: calibration table (given or fitted) or Beta prior + moments."""
    if args.calibration:
        model = Calibration.load(args.calibration)
        logger.info(f"Loaded calibration table from {args.calibration} ({model.n_calls:,} calls)")
        return model, {"site_model": "calibration table", "calibration": str(args.calibration),
                       "n_calls": model.n_calls}

    if args.fit_calibration:
        logger.warning("Fitting the calibration table on this sample. Only meaningful if these reads "
                       "are NOT part of the atlas and mostly atlas-like.")
        model = Calibration()
        for calls, m, n, _ in iter_matched_calls(input_file, atlas, fmt=fmt, chunksize=chunksize, tau=tau):
            model.add(m, n, calls["r"].to_numpy())
        model.save(output_dir / "calibration.json")
        logger.info(f"Calibration table fitted on {model.n_calls:,} calls")
        return model, {"site_model": "calibration table", "calibration": "fitted on sample",
                       "n_calls": model.n_calls}

    prior = parse_prior(args.prior)
    if prior is None:
        prior = fit_beta_prior(atlas.methylated, atlas.total, min_cov=args.prior_min_cov)
        prior_source = f"fitted on atlas sites with coverage >= {args.prior_min_cov}"
    else:
        prior_source = "given"
        logger.info(f"Using Beta prior alpha={prior.alpha}, beta={prior.beta}")

    fit_kwargs = dict(min_cov=args.moments_min_cov, fmt=fmt, chunksize=chunksize, tau=tau)
    if args.moments_from:
        paths = [Path(p).resolve() for p in args.moments_from]
        moments = fit_moments(paths, atlas, prior, **fit_kwargs)
        moments_source = f"fitted on {', '.join(p.name for p in paths)} at atlas coverage >= {args.moments_min_cov}"
    elif args.moments == "fit":
        moments = fit_moments([input_file], atlas, prior, **fit_kwargs)
        moments_source = f"fitted on this sample at atlas coverage >= {args.moments_min_cov}"
    else:
        path = default_moments_path() if args.moments is None else Path(args.moments)
        if not path.exists():
            raise SystemExit(f"moments file not found: {path}. Pass --moments fit, --moments <file> or --moments-from ...")
        loaded = PriorModel.load(path)
        moments = loaded.moments
        moments_source = f"{path}" + (f" ({loaded.info.get('description')})" if loaded.info.get("description") else "")
    logger.info(f"Caller moments: mu1={moments.mu1:.3f} v1={moments.v1:.4f} mu0={moments.mu0:.3f} "
                f"v0={moments.v0:.4f} ({moments_source})", extra={"markup": False})
    model = PriorModel(prior, moments, {"prior_source": prior_source, "moments_source": moments_source})
    model.save(output_dir / "moments.json")
    return model, {"site_model": "beta prior + moments",
                   "prior": {"alpha": prior.alpha, "beta": prior.beta, "source": prior_source},
                   "moments": {**model.to_dict()["moments"], "source": moments_source}}


def summarize_bins(reads: pd.DataFrame, score: str) -> list[dict]:
    rows = []
    for b in N_BIN_LABELS:
        sel = reads["n_bin"].astype(str) == b
        s = reads.loc[sel, score].dropna()
        w = reads.loc[sel, "weight"].dropna()
        rows.append({
            "n_bin": b,
            "n_reads": int(sel.sum()),
            "median_score": float(s.median()) if len(s) else np.nan,
            "score_q01": float(s.quantile(0.01)) if len(s) else np.nan,
            "score_q05": float(s.quantile(0.05)) if len(s) else np.nan,
            "score_sd": float(s.std()) if len(s) else np.nan,
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

    model, model_info = build_site_model(args, atlas, input_file, output_dir, fmt=fmt,
                                         chunksize=chunksize, tau=tau)

    tumor = None
    if getattr(args, "tumor_calibration", None):
        from src.scoring.tumor import TumorModel

        t_atlas = Atlas.from_feather(args.tumor_atlas) if getattr(args, "tumor_atlas", None) else None
        tumor = TumorModel(Calibration.load(args.tumor_calibration), t_atlas, strength=args.tumor_prior_strength)
        model_info["tumor_model"] = {**tumor.describe(), "calibration": str(args.tumor_calibration),
                                     "atlas": str(args.tumor_atlas) if args.tumor_atlas else None}
        logger.info(f"Tumour model: {tumor.describe()}")
    elif args.weight_score == "llr_tumor_per_call":
        raise SystemExit("--weight-score llr_tumor_per_call needs --tumor-calibration (and optionally --tumor-atlas)")

    # ---- scoring ----
    partials, totals, seen = [], Counter(), set()
    for calls, m, n, count in iter_matched_calls(input_file, atlas, fmt=fmt, chunksize=chunksize,
                                                 tau=tau, totals=totals, reads_seen=seen):
        p_tumor = tumor.p(m, n, calls["chrom"], calls["start"], calls["end"]) if tumor else None
        partials.append(score_calls(calls, m, n, count, model, p_tumor=p_tumor))
    reads = finalize_reads(partials)
    logger.info(
        f"Read {totals['n_rows']:,} rows: {totals['n_non_primary']:,} non-primary and "
        f"{totals['n_other_mod']:,} other-modification rows dropped, "
        f"{totals['n_matched']:,} of {totals['n_calls']:,} calls matched the atlas"
    )
    logger.info(f"Scored {len(reads):,} of {len(seen):,} reads (the rest have no atlas CpG)")

    # ---- weights: a plain function of the chosen score, no data-derived cutoffs ----
    score = args.weight_score
    default_center, default_temperature = WEIGHT_DEFAULTS[score]
    center = default_center if args.weight_center is None else args.weight_center
    temperature = default_temperature if args.weight_temperature is None else args.weight_temperature
    thresholds = parse_thresholds(args.weight_thresholds) if args.weight_thresholds else None
    if args.weight_mode == "hard" and thresholds is None:
        raise SystemExit("--weight-mode hard needs --weight-thresholds")
    # compute_weights treats lower values as more abnormal; flip scores that run the other way
    signed = reads[score].to_numpy()
    if score in HIGHER_IS_ABNORMAL:
        signed, center = -signed, -center
    reads["weight"] = compute_weights(signed, reads["n_bin"], mode=args.weight_mode,
                                      center=center, temperature=temperature, w_min=args.weight_min,
                                      thresholds=thresholds)
    if score in HIGHER_IS_ABNORMAL:
        center = -center

    out = reads.reset_index()
    out["n_bin"] = out["n_bin"].astype(str)
    columns = OUTPUT_COLUMNS + [c for c in TUMOR_COLUMNS if c in out.columns]
    pq.write_table(pa.Table.from_pandas(out[columns], preserve_index=False),
                   out_scores, compression="zstd")

    bins = summarize_bins(reads, score)
    summary = {
        "input": str(input_file), "format": fmt, "atlas": str(Path(args.atlas).resolve()),
        "atlas_info": atlas.describe(), "tau": tau, **model_info,
        "weight": {"score": score, "mode": args.weight_mode, "center": center,
                   "temperature": temperature, "w_min": args.weight_min, "thresholds": thresholds},
        "counts": dict(totals), "n_reads_in_file": len(seen), "n_reads_scored": int(len(reads)),
        "bins": bins,
    }
    out_summary.write_text(json.dumps(summary, indent=2, default=str))

    formula = (f"sigmoid(({score} - {center}) / {temperature})" if score in HIGHER_IS_ABNORMAL
               else f"sigmoid(({center} - {score}) / {temperature})")
    logger.info(f"Weight = {formula}; per CpG-count bin: reads, median, 5%, 1%, SD of {score}, mean weight")
    for row in bins:
        logger.info(f"  {row['n_bin']:>4}: {row['n_reads']:>9,}  median {row['median_score']:5.2f}  "
                    f"q05 {row['score_q05']:6.2f}  q01 {row['score_q01']:6.2f}  sd {row['score_sd']:4.2f}  "
                    f"weight {row['mean_weight']:.3f}")
    logger.info("[blue bold]nanoflux score[/] has successfully run!")
    logger.info(f"Per-read scores saved to: {out_scores}")
