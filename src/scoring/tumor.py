"""Tumour model: the alternative to the control atlas in the read likelihood ratio.

Two variants give the tumour methylation probability ``p_tumour`` per call.

``TumorPriorModel`` treats the tumour atlas exactly like the control atlas: a
Beta prior fitted on (or given for) the pooled tumour counts,

    p_tumour = (m_t + alpha_t) / (n_t + alpha_t + beta_t)

A site without tumour reads gets the prior mean. This is the model to use once
the tumour atlas is deep (many high-fraction samples).

``TumorModel`` (calibration variant) exists for sparse tumour atlases. The dense
layer is a calibration table fitted on tumour reads (``Calibration``): for a
site where the control atlas says p (in a coverage stratum), it gives the
probability that a tumour read is methylated there, ``f``. The sparse layer is
the pooled tumour count (m_t of n_t) at the site. They combine as a Beta
posterior with the dense value as prior mean and ``strength`` pseudo-reads:

    p_tumour = (m_t + strength * f) / (n_t + strength)

Either way the read scores built from p_tumour are

    llr_tumor_control = sum over CpGs of  log(p_tumour e^L + 1 - p_tumour) - log(p_control e^L + 1 - p_control)
    llr_tumor_random  = sum over CpGs of  log(p_tumour e^L + 1 - p_tumour) - log(0.5 e^L + 0.5)
    z_tumor           = the z score with p_tumour in place of p_control (same caller moments)

all positive / high when the read is better explained by the tumour model.
"""

from __future__ import annotations

import numpy as np

from src.scoring.atlas import Atlas
from src.scoring.calibration import Calibration
from src.scoring.prior import BetaPrior


class TumorPriorModel:
    """Beta prior on the pooled tumour atlas, mirroring the control-side site model."""

    def __init__(self, prior: BetaPrior, atlas: Atlas, info: dict | None = None) -> None:
        self.prior = prior
        self.atlas = atlas
        self.info = info or {}

    def p(self, m_c: np.ndarray, n_c: np.ndarray, chrom, start, end) -> np.ndarray:
        """Tumour methylation probability per call; the control counts are unused."""
        m_t, n_t, count = self.atlas.lookup(chrom, start, end)
        m_t = np.where(count > 0, m_t, 0.0)
        n_t = np.where(count > 0, n_t, 0.0)
        return self.prior.p(m_t, n_t)

    def describe(self) -> dict:
        return {
            "variant": "beta prior",
            "prior": {"alpha": self.prior.alpha, "beta": self.prior.beta, **self.info},
            "atlas_sites": self.atlas.n_sites,
        }


class TumorModel:
    def __init__(self, calibration: Calibration, atlas: Atlas | None = None, *, strength: float = 2.0) -> None:
        if strength <= 0:
            raise ValueError("strength must be positive")
        self.calibration = calibration
        self.atlas = atlas
        self.strength = strength

    def p(self, m_c: np.ndarray, n_c: np.ndarray, chrom, start, end) -> np.ndarray:
        """Tumour methylation probability per call, given the control counts and the site."""
        f, _, _ = self.calibration.expected(m_c, n_c)
        if self.atlas is None:
            return f
        m_t, n_t, count = self.atlas.lookup(chrom, start, end)
        m_t = np.where(count > 0, m_t, 0.0)
        n_t = np.where(count > 0, n_t, 0.0)
        return (m_t + self.strength * f) / (n_t + self.strength)

    def describe(self) -> dict:
        return {
            "variant": "calibration table",
            "calibration_calls": self.calibration.n_calls,
            "atlas_sites": None if self.atlas is None else self.atlas.n_sites,
            "strength": self.strength,
        }
