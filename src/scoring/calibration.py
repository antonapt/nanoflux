"""Calibration table: what a read calls at a site, given the atlas counts there.

Optional refinement of the Beta prior in ``prior.py``. The atlas gives pooled
counts (m methylated of n reads) per CpG. With shallow pooled coverage the
naive estimate m / n is a biased predictor of what a new read will call:
low-coverage sites are noisy and the genome is mostly methylated, so a low
observed fraction at a 3-read site usually overstates how unmethylated the
site is. A Beta prior fitted to the atlas removes most of that bias; this
table removes the rest, because it assumes no shape for the site distribution.

It must be learned from reads that are NOT part of the atlas (for a control
sample: score it against a reference built without it), otherwise the counts
already contain the reads and the correction cancels out exactly. Calls
are binned by the Laplace-smoothed fraction ``(m + 1) / (n + 2)`` and by a
coverage stratum of n. Per cell it keeps additive sums of the soft call r, so
tables from several samples can be added. From the sums it derives, per cell,

    e_r    E[r]              expected soft call
    v_r    Var[r]            its variance
    p_hat  P(call is methylated)   the calibrated methylation probability

Cells with fewer than ``min_calls`` calls fall back to the p-bin pooled over
all coverage strata, and then to the global values.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

N_P_BINS = 20
COV_EDGES = np.array([5, 10, 20, 50, 100])  # strata: <=5, 6-10, 11-20, 21-50, 51-100, >100
N_COV = len(COV_EDGES) + 1
N_CELLS = N_P_BINS * N_COV
P_HAT_CLIP = 1e-3


def cell_index(m: np.ndarray, n: np.ndarray) -> np.ndarray:
    p = (np.asarray(m, dtype=np.float64) + 1.0) / (np.asarray(n, dtype=np.float64) + 2.0)
    p_bin = np.clip((p * N_P_BINS).astype(np.int64), 0, N_P_BINS - 1)
    cov_bin = np.searchsorted(COV_EDGES, np.asarray(n, dtype=np.float64), side="left")
    return p_bin * N_COV + cov_bin


class Calibration:
    def __init__(
        self,
        counts: np.ndarray | None = None,
        sum_r: np.ndarray | None = None,
        sum_r2: np.ndarray | None = None,
        sum_hard: np.ndarray | None = None,
        *,
        min_calls: int = 200,
    ) -> None:
        z = lambda: np.zeros(N_CELLS, dtype=np.float64)  # noqa: E731
        self.counts = z() if counts is None else np.asarray(counts, dtype=np.float64)
        self.sum_r = z() if sum_r is None else np.asarray(sum_r, dtype=np.float64)
        self.sum_r2 = z() if sum_r2 is None else np.asarray(sum_r2, dtype=np.float64)
        self.sum_hard = z() if sum_hard is None else np.asarray(sum_hard, dtype=np.float64)
        self.min_calls = min_calls
        self._tables: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None

    # ------------------------------------------------------------------ accumulation
    def add(self, m: np.ndarray, n: np.ndarray, r: np.ndarray) -> None:
        c = cell_index(m, n)
        r = np.asarray(r, dtype=np.float64)
        self.counts += np.bincount(c, minlength=N_CELLS)
        self.sum_r += np.bincount(c, weights=r, minlength=N_CELLS)
        self.sum_r2 += np.bincount(c, weights=r * r, minlength=N_CELLS)
        self.sum_hard += np.bincount(c, weights=(r > 0.5).astype(np.float64), minlength=N_CELLS)
        self._tables = None

    def __add__(self, other: "Calibration") -> "Calibration":
        return Calibration(
            self.counts + other.counts, self.sum_r + other.sum_r,
            self.sum_r2 + other.sum_r2, self.sum_hard + other.sum_hard,
            min_calls=self.min_calls,
        )

    @property
    def n_calls(self) -> int:
        return int(self.counts.sum())

    # ------------------------------------------------------------------ tables
    def _fallback(self, sums: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Per-cell (sum, count) with sparse cells replaced by their p-bin, then global."""
        counts = self.counts.reshape(N_P_BINS, N_COV)
        sums = sums.reshape(N_P_BINS, N_COV)
        bin_counts = counts.sum(axis=1, keepdims=True)
        bin_sums = sums.sum(axis=1, keepdims=True)
        use_cell = counts >= self.min_calls
        use_bin = ~use_cell & (bin_counts >= self.min_calls)
        out_sum = np.where(use_cell, sums, np.where(use_bin, bin_sums, sums.sum()))
        out_cnt = np.where(use_cell, counts, np.where(use_bin, bin_counts, counts.sum()))
        return out_sum.ravel(), out_cnt.ravel()

    def tables(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Per-cell ``(e_r, v_r, p_hat)``."""
        if self._tables is None:
            if self.n_calls == 0:
                raise ValueError("calibration table is empty")
            s1, c = self._fallback(self.sum_r)
            s2, _ = self._fallback(self.sum_r2)
            sh, _ = self._fallback(self.sum_hard)
            e_r = s1 / c
            v_r = np.maximum(s2 / c - e_r**2, 0.0)
            p_hat = np.clip(sh / c, P_HAT_CLIP, 1 - P_HAT_CLIP)
            self._tables = (e_r, v_r, p_hat)
        return self._tables

    def expected(self, m: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(p_hat, e_r, v_r)`` for every call."""
        e_r, v_r, p_hat = self.tables()
        c = cell_index(m, n)
        return p_hat[c], e_r[c], v_r[c]

    # ------------------------------------------------------------------ io
    def to_dict(self) -> dict:
        e_r, v_r, p_hat = self.tables()
        return {
            "n_p_bins": N_P_BINS,
            "cov_edges": COV_EDGES.tolist(),
            "min_calls": self.min_calls,
            "n_calls": self.n_calls,
            "counts": self.counts.tolist(),
            "sum_r": self.sum_r.tolist(),
            "sum_r2": self.sum_r2.tolist(),
            "sum_hard": self.sum_hard.tolist(),
            "p_hat": np.round(p_hat, 4).reshape(N_P_BINS, N_COV).tolist(),
            "e_r": np.round(e_r, 4).reshape(N_P_BINS, N_COV).tolist(),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Calibration":
        if d.get("n_p_bins") != N_P_BINS or list(d.get("cov_edges", [])) != COV_EDGES.tolist():
            raise ValueError("calibration file was written with different bin edges")
        return cls(d["counts"], d["sum_r"], d["sum_r2"], d["sum_hard"], min_calls=d["min_calls"])

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=1))

    @classmethod
    def load(cls, path: str | Path) -> "Calibration":
        return cls.from_dict(json.loads(Path(path).read_text()))
