"""Pooled CpG reference atlas: per-site methylated and total read counts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.feather as feather

from src.utils.log import logger

DEFAULT_COLUMNS = ("chrom_x", "start_genomic", "methylated_pooled", "total_pooled")


class Atlas:
    """Sorted lookup of pooled counts by (chromosome, 1-based plus-strand C position).

    Sites are keyed as ``(chromosome code << 32) | position`` and kept sorted, so a
    single CpG and a run of consecutive CpGs (a grouped nanopolish call) use the
    same range lookup: counts are averaged over every atlas site inside the range.
    """

    def __init__(
        self,
        chrom: np.ndarray,
        pos: np.ndarray,
        methylated: np.ndarray,
        total: np.ndarray,
        *,
        min_cov: int = 1,
    ) -> None:
        pos = np.asarray(pos, dtype=np.int64)
        methylated = np.asarray(methylated, dtype=np.int64)
        total = np.asarray(total, dtype=np.int64)
        if not (len(chrom) == len(pos) == len(methylated) == len(total)):
            raise ValueError("atlas columns have different lengths")

        self.n_sites_total = len(total)
        keep = total >= min_cov
        self.min_cov = min_cov

        codes, names = pd.factorize(np.asarray(chrom)[keep])
        self.chrom_codes: dict[str, int] = {str(c): i for i, c in enumerate(names)}
        key = (codes.astype(np.int64) << 32) | pos[keep]
        order = np.argsort(key, kind="stable")
        self.key = key[order]
        if len(self.key) > 1 and np.any(np.diff(self.key) == 0):
            raise ValueError("atlas contains duplicate (chromosome, position) sites")
        self.methylated = methylated[keep][order]
        self.total = total[keep][order]
        self._cum_m = np.concatenate([[0], np.cumsum(self.methylated)])
        self._cum_n = np.concatenate([[0], np.cumsum(self.total)])

    @property
    def n_sites(self) -> int:
        return len(self.key)

    @classmethod
    def from_feather(
        cls,
        path: str | Path,
        *,
        columns: tuple[str, str, str, str] = DEFAULT_COLUMNS,
        min_cov: int = 1,
    ) -> "Atlas":
        chrom_col, pos_col, meth_col, total_col = columns
        table = feather.read_table(path, columns=list(columns))
        # dictionary-encode the chromosome column so 30M strings never become python objects
        chrom = pc.dictionary_encode(table[chrom_col]).combine_chunks()
        names = np.asarray(chrom.dictionary.to_pylist(), dtype=object)
        chrom_arr = names[chrom.indices.to_numpy()]
        atlas = cls(
            chrom_arr,
            table[pos_col].to_numpy(),
            table[meth_col].to_numpy(),
            table[total_col].to_numpy(),
            min_cov=min_cov,
        )
        logger.info(
            f"Loaded atlas with {atlas.n_sites:,} sites "
            f"({atlas.n_sites_total - atlas.n_sites:,} below coverage {min_cov} dropped)"
        )
        return atlas

    def _codes(self, chrom: np.ndarray | pd.Series) -> np.ndarray:
        cat = pd.Categorical(np.asarray(chrom))
        per_cat = np.asarray([self.chrom_codes.get(str(c), -1) for c in cat.categories], dtype=np.int64)
        codes = per_cat[cat.codes] if len(per_cat) else np.full(len(cat), -1, dtype=np.int64)
        codes[cat.codes < 0] = -1
        return codes

    def lookup(
        self, chrom, start, end
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Mean methylated and total counts over atlas sites in [start, end] (1-based, inclusive).

        Returns ``(m, n, count)``; ``count`` is the number of atlas sites in the
        range and is 0 where nothing matched (then m and n are 0).
        """
        start = np.asarray(start, dtype=np.int64)
        end = np.asarray(end, dtype=np.int64)
        codes = self._codes(chrom)
        known = codes >= 0
        base = np.where(known, codes, 0).astype(np.int64) << 32
        i0 = np.searchsorted(self.key, base | start, side="left")
        i1 = np.searchsorted(self.key, base | end, side="right")
        count = np.where(known, i1 - i0, 0)
        safe = np.maximum(count, 1)
        m = (self._cum_m[i1] - self._cum_m[i0]) / safe
        n = (self._cum_n[i1] - self._cum_n[i0]) / safe
        m[count == 0] = 0.0
        n[count == 0] = 0.0
        return m, n, count

    def describe(self) -> dict:
        return {
            "n_sites": self.n_sites,
            "n_sites_total": self.n_sites_total,
            "min_cov": self.min_cov,
            "n_chromosomes": len(self.chrom_codes),
        }
