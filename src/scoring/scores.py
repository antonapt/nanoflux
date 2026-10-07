"""Per-read scores against the control atlas and, optionally, a tumour model.

Per call, with the site methylation probability ``p`` and the expected soft
call ``e_r`` with variance ``v_r`` (both from the site model: a Beta prior
plus caller moments, or a calibration table) and the read's soft call ``r``:

    obs = log(1 - p) + r   * logit(p)      observed log-likelihood under the reference
    exp = log(1 - p) + e_r * logit(p)      its expectation for a read that follows the reference
    var = v_r * logit(p)**2                its variance
    llr = log(p e^L + 1 - p) - log(p_alt e^L + 1 - p_alt),  L = logit(r)

Summed over the calls of one read, with p = p_control,

    z_control                    = (sum obs - sum exp) / sqrt(sum var)
    llr_control_random_per_call  = sum llr / n_calls              (p_alt = 0.5)

and, when a tumour model gives p_tumour per call (same caller moments),

    z_tumor                      = the same z with p_tumour
    llr_tumor_random_per_call    = sum [log(p_t e^L + 1 - p_t) - log(0.5 e^L + 0.5)] / n_calls
    llr_tumor_control_per_call   = sum [log(p_t e^L + 1 - p_t) - log(p_c e^L + 1 - p_c)] / n_calls

Direction: the two control scores are LOW for a read that does not look like the
control reference; the three tumour scores are HIGH for a read that looks like
tumour. ``TUMOR_LIKE_IS_HIGH`` records this. Old column names (``z``,
``llr_per_call``, ``llr_tumor_per_call``) are accepted as aliases.

Sums are accumulated per chunk so memory scales with the number of reads, not calls.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np
import pandas as pd

from src.scoring.weights import assign_bins


class SiteModel(Protocol):
    def expected(self, m: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(p, e_r, v_r)`` per call."""
        ...


EPS = 1e-6

CONTROL_SCORES = ["z_control", "llr_control_random_per_call"]
TUMOR_SCORES = ["z_tumor", "llr_tumor_random_per_call", "llr_tumor_control_per_call"]
ALL_SCORES = CONTROL_SCORES + TUMOR_SCORES
# scores where a HIGHER value means more tumour-like; the others are LOWER = less like control
TUMOR_LIKE_IS_HIGH = set(TUMOR_SCORES)
# old column names -> new
SCORE_ALIASES = {"z": "z_control", "llr_per_call": "llr_control_random_per_call",
                 "llr_tumor_per_call": "llr_tumor_control_per_call"}


def canonical_score(name: str) -> str:
    return SCORE_ALIASES.get(name, name)


def rename_score_columns(table: pd.DataFrame) -> pd.DataFrame:
    """Give an old read_scores table the current column names (no-op on new ones)."""
    ren = {k: v for k, v in SCORE_ALIASES.items() if k in table.columns and v not in table.columns}
    return table.rename(columns=ren) if ren else table


_AGG = {
    "n_calls": ("obs_c", "size"),
    "n_cpg": ("n_cpg", "sum"),
    "obs_c": ("obs_c", "sum"),
    "exp_c": ("exp_c", "sum"),
    "var_c": ("var_c", "sum"),
    "llr_c": ("llr_c", "sum"),
}
_AGG_TUMOR = {
    "obs_t": ("obs_t", "sum"),
    "exp_t": ("exp_t", "sum"),
    "var_t": ("var_t", "sum"),
    "llr_t": ("llr_t", "sum"),
}


def _logit(p: np.ndarray) -> np.ndarray:
    return np.log(p) - np.log1p(-p)


def moments_expected(p: np.ndarray, moments) -> tuple[np.ndarray, np.ndarray]:
    """``(e_r, v_r)`` of the soft call at sites with probability p, from the caller moments."""
    e_r = p * moments.mu1 + (1 - p) * moments.mu0
    v_r = p * (moments.v1 + moments.mu1**2) + (1 - p) * (moments.v0 + moments.mu0**2) - e_r**2
    return e_r, np.maximum(v_r, 0.0)


def _terms(p: np.ndarray, e_r: np.ndarray, v_r: np.ndarray, r: np.ndarray, L: np.ndarray,
           p_alt: float) -> dict[str, np.ndarray]:
    p = np.clip(p, EPS, 1 - EPS)
    lp = _logit(p)
    ll = np.logaddexp(np.log(p) + L, np.log1p(-p))
    return {
        "obs": np.log1p(-p) + r * lp,
        "exp": np.log1p(-p) + e_r * lp,
        "var": v_r * lp**2,
        "ll": ll,
        "llr": ll - np.logaddexp(np.log(p_alt) + L, np.log1p(-p_alt)),
    }


def score_calls(
    calls: pd.DataFrame,
    m: np.ndarray,
    n: np.ndarray,
    n_cpg: np.ndarray,
    model: SiteModel,
    *,
    p_alt: float = 0.5,
    p_tumor: np.ndarray | None = None,
    tumor_moments=None,
) -> pd.DataFrame:
    """Per-read partial sums for one chunk of atlas-matched calls.

    ``p_tumor`` (per call) adds the tumour terms. ``tumor_moments`` (caller
    moments) gives the expectation and variance of the soft call under the
    tumour reference, needed for z_tumor; without them z_tumor is NaN.
    """
    p_c, e_c, v_c = model.expected(m, n)
    r = np.clip(calls["r"].to_numpy(np.float64), EPS, 1 - EPS)
    L = _logit(r)
    c = _terms(p_c, e_c, v_c, r, L, p_alt)
    per_call = pd.DataFrame(
        {
            "read_id": calls["read_id"].to_numpy(),
            "n_cpg": np.asarray(n_cpg, dtype=np.int64),
            "obs_c": c["obs"], "exp_c": c["exp"], "var_c": c["var"], "llr_c": c["llr"],
        }
    )
    agg = dict(_AGG)
    if p_tumor is not None:
        p_t = np.clip(np.asarray(p_tumor, dtype=np.float64), EPS, 1 - EPS)
        if tumor_moments is not None:
            e_t, v_t = moments_expected(p_t, tumor_moments)
        else:
            e_t, v_t = np.full_like(p_t, np.nan), np.full_like(p_t, np.nan)
        t = _terms(p_t, e_t, v_t, r, L, p_alt)
        per_call["obs_t"], per_call["exp_t"], per_call["var_t"] = t["obs"], t["exp"], t["var"]
        per_call["llr_t"] = t["llr"]                      # tumour vs random
        per_call["llr_tc"] = t["ll"] - c["ll"]            # tumour vs control
        agg.update(_AGG_TUMOR)
        agg["llr_tc"] = ("llr_tc", "sum")
    return per_call.groupby("read_id", sort=False).agg(**agg)


def finalize_reads(partials: list[pd.DataFrame]) -> pd.DataFrame:
    """Combine chunk partials into one row per read with the scores and the CpG bin."""
    if not partials:
        raise ValueError("no scored calls")
    reads = pd.concat(partials).groupby(level=0, sort=False).sum(min_count=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        reads["z_control"] = (reads["obs_c"] - reads["exp_c"]) / np.sqrt(reads["var_c"])
        reads["llr_control_random_per_call"] = reads["llr_c"] / reads["n_calls"]
        if "llr_t" in reads:
            reads["z_tumor"] = (reads["obs_t"] - reads["exp_t"]) / np.sqrt(reads["var_t"])
            reads["llr_tumor_random_per_call"] = reads["llr_t"] / reads["n_calls"]
            reads["llr_tumor_control_per_call"] = reads["llr_tc"] / reads["n_calls"]
    reads["n_bin"] = assign_bins(reads["n_cpg"].to_numpy())
    reads.index.name = "read_id"
    cols = ["n_calls", "n_cpg", "n_bin"] + CONTROL_SCORES
    if "llr_t" in reads:
        cols += TUMOR_SCORES
    return reads[cols]
