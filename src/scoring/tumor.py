"""Tumour model: the alternative to the control atlas in the read likelihood ratio.

Per-site tumour counts from a few samples are sparse, so the model has two
layers. The dense layer is a calibration table fitted on tumour reads
(``Calibration``): for a site where the control atlas says p (in a coverage
stratum), it gives the probability that a tumour read is methylated there,
``f``. The sparse layer is the pooled tumour count (m_t of n_t) at the site,
where the samples have coverage. They combine as a Beta posterior with the
dense value as prior mean and ``strength`` pseudo-reads:

    p_tumour = (m_t + strength * f) / (n_t + strength)

so a site without tumour reads gets f, and a well-covered site its own count.
The read score is then

    llr_tumor = sum over CpGs of  log(p_tumour e^L + 1 - p_tumour) - log(p_control e^L + 1 - p_control)

positive when the read is better explained by the tumour model.
"""

from __future__ import annotations

import numpy as np

from src.scoring.atlas import Atlas
from src.scoring.calibration import Calibration


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
            "calibration_calls": self.calibration.n_calls,
            "atlas_sites": None if self.atlas is None else self.atlas.n_sites,
            "strength": self.strength,
        }
