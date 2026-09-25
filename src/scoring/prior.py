"""Site model without a calibration table: a Beta prior plus caller moments.

The atlas gives pooled counts (m of n) per CpG. Shallow sites are read with a
credibility weight through a Beta prior,

    p = (m + alpha) / (n + alpha + beta)

which pools the real reads with alpha methylated and beta unmethylated
pseudo-reads. With ``alpha = beta = 1`` the pull goes toward 0.5, but the
genome is mostly methylated, so the prior should pull toward the typical site.
``fit_beta_prior`` estimates alpha and beta from the high-coverage sites of the
atlas itself (method of moments for a beta-binomial), so no reads are needed
and nothing leaks from the samples being scored.

To score a read we also need what the caller does at a truly methylated or
unmethylated CpG: the class moments ``mu1, v1`` (mean and variance of the soft
call r when methylated) and ``mu0, v0`` (when unmethylated). Under the model
``E[r | p]`` and ``E[r^2 | p]`` are linear in p, so two least-squares fits of
r and r^2 on [1, p] over a sample's calls recover them. The normal equations
are additive over chunks and samples.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.utils.log import logger


@dataclass
class BetaPrior:
    alpha: float
    beta: float

    def p(self, m: np.ndarray, n: np.ndarray) -> np.ndarray:
        m = np.asarray(m, dtype=np.float64)
        n = np.asarray(n, dtype=np.float64)
        return (m + self.alpha) / (n + self.alpha + self.beta)

    @property
    def mean(self) -> float:
        return self.alpha / (self.alpha + self.beta)


def fit_beta_prior(methylated: np.ndarray, total: np.ndarray, *, min_cov: int = 100,
                   min_sites: int = 1000) -> BetaPrior:
    """Beta prior matching the mean and spread of site methylation at well-covered sites.

    The observed fractions f = m / n scatter more than the true site
    probabilities because of binomial sampling; the method-of-moments
    estimate removes that part: Var(p) = (Var(f) - mu(1-mu) E[1/n]) / (1 - E[1/n]).
    """
    m = np.asarray(methylated, dtype=np.float64)
    n = np.asarray(total, dtype=np.float64)
    hi = n >= min_cov
    while hi.sum() < min_sites and min_cov > 1:
        min_cov = max(1, min_cov // 2)
        hi = n >= min_cov
    if hi.sum() < 2:
        raise ValueError("not enough atlas sites to fit a Beta prior")
    f = m[hi] / n[hi]
    mu = float(f.mean())
    e_inv_n = float((1.0 / n[hi]).mean())
    var_p = (float(f.var()) - mu * (1 - mu) * e_inv_n) / max(1.0 - e_inv_n, 1e-9)
    var_p = float(np.clip(var_p, 1e-4, mu * (1 - mu) * 0.999))
    k = mu * (1 - mu) / var_p - 1.0
    prior = BetaPrior(mu * k, (1 - mu) * k)
    logger.info(
        f"Fitted Beta prior on {int(hi.sum()):,} atlas sites with coverage >= {min_cov}: "
        f"alpha={prior.alpha:.3f}, beta={prior.beta:.3f} (pulls toward {prior.mean:.2f})"
    )
    return prior


@dataclass
class Moments:
    mu1: float
    v1: float
    mu0: float
    v0: float
    n_calls: int


TAIL_HI, TAIL_LO = 0.9, 0.1   # "near-certain" atlas sites used for the variances
V_FLOOR = 0.002               # never claim a caller more precise than this


class MomentStats:
    """Additive sufficient statistics for the caller moments.

    The means come from the regression of r on [1, p] (``E[r | p]`` is linear
    in p under the model). The variances from the matching regression of r^2
    are fragile, because p is only an estimate of the true site probability;
    instead they are measured directly at near-certain sites (p >= 0.9 for
    v1, p <= 0.1 for v0), minus the small residual mixing term
    ``p (1 - p) (mu1 - mu0)^2``, and floored at ``V_FLOOR``. Fit the
    statistics only at atlas sites with high pooled coverage when the sample
    is part of the atlas, so that its own reads do not dominate p.
    """

    def __init__(self) -> None:
        self.A = np.zeros((2, 2))
        self.b1 = np.zeros(2)
        self.b2 = np.zeros(2)
        self.n = 0
        self.hi = np.zeros(4)  # count, sum r, sum r^2, sum p(1-p)  at p >= TAIL_HI
        self.lo = np.zeros(4)  # same at p <= TAIL_LO

    def add(self, p: np.ndarray, r: np.ndarray) -> None:
        p = np.asarray(p, dtype=np.float64)
        r = np.asarray(r, dtype=np.float64)
        X = np.column_stack([np.ones(len(p)), p])
        self.A += X.T @ X
        self.b1 += X.T @ r
        self.b2 += X.T @ (r * r)
        self.n += len(p)
        for acc, sel in ((self.hi, p >= TAIL_HI), (self.lo, p <= TAIL_LO)):
            rs, ps = r[sel], p[sel]
            acc += (len(rs), rs.sum(), (rs * rs).sum(), (ps * (1 - ps)).sum())

    def __add__(self, other: "MomentStats") -> "MomentStats":
        out = MomentStats()
        out.A, out.b1, out.b2, out.n = self.A + other.A, self.b1 + other.b1, self.b2 + other.b2, self.n + other.n
        out.hi, out.lo = self.hi + other.hi, self.lo + other.lo
        return out

    def fit(self, *, min_calls: int = 100, min_tail: int = 1000) -> Moments:
        if self.n < min_calls:
            raise ValueError(f"only {self.n} calls available for the moment fit")
        mu0, slope1 = np.linalg.solve(self.A, self.b1)
        M0, slope2 = np.linalg.solve(self.A, self.b2)
        mu1, M1 = mu0 + slope1, M0 + slope2

        def tail_var(acc: np.ndarray, fallback: float, label: str) -> float:
            c, s1, s2, spq = acc
            if c < min_tail:
                logger.warning(f"only {int(c)} calls at {label} atlas sites; variance from the regression")
                return fallback
            return s2 / c - (s1 / c) ** 2 - spq / c * (mu1 - mu0) ** 2

        v1 = tail_var(self.hi, M1 - mu1**2, f"p >= {TAIL_HI}")
        v0 = tail_var(self.lo, M0 - mu0**2, f"p <= {TAIL_LO}")
        return Moments(float(mu1), float(max(v1, V_FLOOR)), float(mu0), float(max(v0, V_FLOOR)), int(self.n))


class PriorModel:
    """``expected(m, n) -> (p, e_r, v_r)`` from a Beta prior and caller moments."""

    def __init__(self, prior: BetaPrior, moments: Moments, info: dict | None = None) -> None:
        self.prior = prior
        self.moments = moments
        self.info = info or {}

    def expected(self, m: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        p = self.prior.p(m, n)
        mo = self.moments
        e_r = p * mo.mu1 + (1 - p) * mo.mu0
        v_r = p * (mo.v1 + mo.mu1**2) + (1 - p) * (mo.v0 + mo.mu0**2) - e_r**2
        return p, e_r, np.maximum(v_r, 0.0)

    def to_dict(self) -> dict:
        return {"prior": asdict(self.prior), "moments": asdict(self.moments), "info": self.info}

    @classmethod
    def from_dict(cls, d: dict) -> "PriorModel":
        return cls(BetaPrior(**d["prior"]), Moments(**d["moments"]), d.get("info"))

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def load(cls, path: str | Path) -> "PriorModel":
        return cls.from_dict(json.loads(Path(path).read_text()))
