"""Per-read scores against the atlas.

Per call, with the site methylation probability ``p`` and the expected soft
call ``e_r`` with variance ``v_r`` (both from the site model: a Beta prior
plus caller moments, or a calibration table) and the read's soft call ``r``:

    obs = log(1 - p) + r   * logit(p)      observed log-likelihood under the atlas
    exp = log(1 - p) + e_r * logit(p)      its expectation for a read that follows the atlas
    var = v_r * logit(p)**2                its variance
    llr = log(p e^L + 1 - p) - log(p_alt e^L + 1 - p_alt),  L = logit(r)

Summed over the calls of one read,

    z            = (sum obs - sum exp) / sqrt(sum var)
    llr_per_call = sum llr / n_calls

``z`` is about 0 for a read that looks like the atlas and strongly negative for
one that does not. Sums are accumulated per chunk so memory scales with the
number of reads, not calls.
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

_AGG = {
    "n_calls": ("obs", "size"),
    "n_cpg": ("n_cpg", "sum"),
    "obs": ("obs", "sum"),
    "exp": ("exp", "sum"),
    "var": ("var", "sum"),
    "llr": ("llr", "sum"),
}


def _logit(p: np.ndarray) -> np.ndarray:
    return np.log(p) - np.log1p(-p)


def score_calls(
    calls: pd.DataFrame,
    m: np.ndarray,
    n: np.ndarray,
    n_cpg: np.ndarray,
    model: SiteModel,
    *,
    p_alt: float = 0.5,
) -> pd.DataFrame:
    """Per-read partial sums for one chunk of atlas-matched calls."""
    p, e_r, v_r = model.expected(m, n)
    p = np.clip(p, EPS, 1 - EPS)
    lp = _logit(p)
    r = np.clip(calls["r"].to_numpy(np.float64), EPS, 1 - EPS)
    L = _logit(r)
    per_call = pd.DataFrame(
        {
            "read_id": calls["read_id"].to_numpy(),
            "n_cpg": np.asarray(n_cpg, dtype=np.int64),
            "obs": np.log1p(-p) + r * lp,
            "exp": np.log1p(-p) + e_r * lp,
            "var": v_r * lp**2,
            "llr": np.logaddexp(np.log(p) + L, np.log1p(-p))
            - np.logaddexp(np.log(p_alt) + L, np.log1p(-p_alt)),
        }
    )
    return per_call.groupby("read_id", sort=False).agg(**_AGG)


def finalize_reads(partials: list[pd.DataFrame]) -> pd.DataFrame:
    """Combine chunk partials into one row per read with z, llr_per_call and the CpG bin."""
    if not partials:
        raise ValueError("no scored calls")
    reads = pd.concat(partials).groupby(level=0, sort=False).sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        reads["z"] = (reads["obs"] - reads["exp"]) / np.sqrt(reads["var"])
    reads["llr_per_call"] = reads["llr"] / reads["n_calls"]
    reads["n_bin"] = assign_bins(reads["n_cpg"].to_numpy())
    reads.index.name = "read_id"
    return reads[["n_calls", "n_cpg", "n_bin", "z", "llr_per_call"]]
