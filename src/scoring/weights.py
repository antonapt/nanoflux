"""Read weights from read scores.

Reads are binned by their number of atlas CpGs, because a short read cannot
reach an extreme score. Per bin a cutoff ``c`` marks abnormal reads (lower z is
more abnormal). Weights follow one family for both modes:

    soft:  w = w_min + (1 - w_min) * sigmoid((c - z) / s)
    hard:  w = 1 if z < c else w_min

with the scale ``s = temperature * (median_z - c)`` of the bin, so that at
temperature 1 a typical read (z at the bin median) gets about 0.27, a read at
the cutoff 0.5 and a clearly abnormal read 1. Temperature 0 is the hard mode.
``w_min`` is the floor for reads that look like the atlas; with ``w_min = 0``
the hard mode is a plain filter.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

N_BIN_EDGES = [0, 1, 2, 4, 9, np.inf]
N_BIN_LABELS = ["1", "2", "3-4", "5-9", "10+"]
_MIN_SCALE = 1e-3


def assign_bins(n_cpg: np.ndarray) -> pd.Categorical:
    return pd.cut(np.asarray(n_cpg), N_BIN_EDGES, labels=N_BIN_LABELS)


def parse_thresholds(text: str) -> dict[str, float]:
    """``"1:-2.5,2:-3,3-4:-3.5,5-9:-4,10+:-6"`` -> per-bin cutoff. All bins are required."""
    out: dict[str, float] = {}
    for item in text.split(","):
        if ":" not in item:
            raise ValueError(f"bad threshold entry {item!r}; expected BIN:VALUE")
        label, value = item.rsplit(":", 1)
        label = label.strip()
        if label not in N_BIN_LABELS:
            raise ValueError(f"unknown bin {label!r}; bins are {N_BIN_LABELS}")
        out[label] = float(value)
    missing = [b for b in N_BIN_LABELS if b not in out]
    if missing:
        raise ValueError(f"thresholds missing for bins {missing}")
    return out


def bin_quantiles(z: np.ndarray, n_bin, q: float) -> dict[str, float]:
    s = pd.Series(np.asarray(z, dtype=np.float64))
    labels = pd.Series(np.asarray(n_bin).astype(str))
    return {b: float(s[labels == b].dropna().quantile(q)) if (labels == b).any() else np.nan
            for b in N_BIN_LABELS}


def compute_weights(
    z: np.ndarray,
    n_bin,
    cutoffs: dict[str, float],
    medians: dict[str, float],
    *,
    mode: str = "soft",
    temperature: float = 1.0,
    w_min: float = 0.0,
) -> np.ndarray:
    if mode not in ("soft", "hard"):
        raise ValueError(f"unknown weight mode {mode!r}")
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min must be in [0, 1]")
    z = np.asarray(z, dtype=np.float64)
    labels = np.asarray(n_bin).astype(str)
    w = np.full(len(z), np.nan)
    for b in N_BIN_LABELS:
        sel = labels == b
        if not sel.any():
            continue
        c, med = cutoffs[b], medians[b]
        if np.isnan(c) or np.isnan(med):
            continue
        zb = z[sel]
        if mode == "hard" or temperature <= 0:
            wb = np.where(zb < c, 1.0, w_min)
        else:
            scale = max(temperature * (med - c), _MIN_SCALE)
            wb = w_min + (1.0 - w_min) / (1.0 + np.exp(-(c - zb) / scale))
        w[sel] = wb
    w[np.isnan(z)] = np.nan
    return w
