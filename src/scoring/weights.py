"""Read weights from read scores.

The weight is a plain decreasing function of z (lower z = less like the
reference = higher weight). No cutoffs are derived from the data:

    sigmoid:  w = w_min + (1 - w_min) * sigmoid((center - z) / temperature)
    hard:     w = 1 if z < threshold[bin] else w_min   (explicit per-bin thresholds)

With the defaults (center 0, temperature 1) a read that fits the reference
(z near 0) gets about 0.5, a read at z = -2 about 0.88 and at z = -4 about
0.98. A smaller temperature sharpens the contrast. In a weighted pileup only
the ratio of weights between reads at the same site matters, so the absolute
level is irrelevant; ``w_min`` keeps the most reference-like reads from
vanishing entirely.

Reads are binned by their number of atlas CpGs for reporting and for the
hard mode, because a short read cannot reach an extreme score.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

N_BIN_EDGES = [0, 1, 2, 4, 9, np.inf]
N_BIN_LABELS = ["1", "2", "3-4", "5-9", "10+"]
WEIGHT_MODES = ("sigmoid", "hard")


def assign_bins(n_cpg: np.ndarray) -> pd.Categorical:
    return pd.cut(np.asarray(n_cpg), N_BIN_EDGES, labels=N_BIN_LABELS)


def parse_thresholds(text: str) -> dict[str, float]:
    """``"1:-2.5,2:-3,3-4:-3.5,5-9:-4,10+:-6"`` -> per-bin z threshold. All bins are required."""
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
    n_bin=None,
    *,
    mode: str = "sigmoid",
    center: float = 0.0,
    temperature: float = 1.0,
    w_min: float = 0.0,
    thresholds: dict[str, float] | None = None,
) -> np.ndarray:
    if mode not in WEIGHT_MODES:
        raise ValueError(f"unknown weight mode {mode!r}; choose from {WEIGHT_MODES}")
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min must be in [0, 1]")
    z = np.asarray(z, dtype=np.float64)
    if mode == "sigmoid":
        if temperature <= 0:
            raise ValueError("temperature must be positive")
        w = w_min + (1.0 - w_min) / (1.0 + np.exp((z - center) / temperature))
    else:
        if thresholds is None or n_bin is None:
            raise ValueError("hard mode needs per-bin thresholds and the CpG bins")
        labels = np.asarray(n_bin).astype(str)
        thr = np.array([thresholds[b] for b in labels], dtype=np.float64)
        w = np.where(z < thr, 1.0, w_min)
    w[np.isnan(z)] = np.nan
    return w
