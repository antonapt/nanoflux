"""Adapters turning per-read methylation call tables into one canonical layout.

Every adapter yields ``(calls, stats)`` per chunk, where ``calls`` has the columns

    read_id  read name
    chrom    chromosome
    start    1-based position of the first CpG's plus-strand C
    end      1-based position of the last CpG's plus-strand C (== start for single calls)
    r        the read's probability that the CpG(s) are methylated

so the scoring code never needs to know which caller produced the file. A caller
that reports a log-likelihood ratio is converted with ``r = sigmoid(tau * LLR)``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd

CANONICAL_COLUMNS = ["read_id", "chrom", "start", "end", "r"]

# SAM flag bits: unmapped, secondary, supplementary
_NON_PRIMARY = 0x4 | 0x100 | 0x800


def _expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def iter_modkit_calls(
    path: str | Path, *, chunksize: int, tau: float = 1.0
) -> Iterator[tuple[pd.DataFrame, dict]]:
    """``modkit extract --read-calls-path`` table (one row per read and CpG).

    ``ref_position`` is 0-based and refers to the modified base itself. On a
    plus-strand read that is the C, so the 1-based plus-strand C is
    ``ref_position + 1``. On a minus-strand read the modified base is the C of
    the reverse strand, which sits on the plus-strand G one base to the right,
    so ``ref_position`` already equals the 1-based position of the plus-strand C.

    Only primary alignments and the calls ``m`` (5mC) / ``-`` (canonical) are
    kept; ``h`` (5hmC) has no counterpart in the atlas. ``call_prob`` is the
    probability of the *called* class, so the probability of methylation is
    ``call_prob`` for ``m`` and ``1 - call_prob`` for ``-``. ``tau`` is unused.
    """
    usecols = ["read_id", "ref_position", "chrom", "ref_strand", "call_prob", "call_code", "flag"]
    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=chunksize,
        dtype={
            "read_id": "string",
            "ref_position": np.int64,
            "chrom": "category",
            "ref_strand": "category",
            "call_prob": np.float64,
            "call_code": "category",
            "flag": np.int64,
        },
    )
    for chunk in reader:
        primary = (chunk["flag"].to_numpy() & _NON_PRIMARY) == 0
        code = chunk["call_code"].astype(str).to_numpy()
        is_m = code == "m"
        is_canonical = code == "-"
        keep = primary & (is_m | is_canonical)

        sub = chunk.loc[keep]
        plus = sub["ref_strand"].astype(str).to_numpy() == "+"
        pos = sub["ref_position"].to_numpy(np.int64) + plus.astype(np.int64)
        prob = sub["call_prob"].to_numpy(np.float64)
        r = np.where(is_m[keep], prob, 1.0 - prob)
        calls = pd.DataFrame(
            {
                "read_id": sub["read_id"].to_numpy(),
                "chrom": sub["chrom"].astype(str).to_numpy(),
                "start": pos,
                "end": pos,
                "r": r,
            }
        )
        stats = {
            "n_rows": int(len(chunk)),
            "n_non_primary": int((~primary).sum()),
            "n_other_mod": int((primary & ~(is_m | is_canonical)).sum()),
            "n_calls": int(len(calls)),
        }
        yield calls, stats


def iter_nanopolish_calls(
    path: str | Path, *, chunksize: int, tau: float = 1.0
) -> Iterator[tuple[pd.DataFrame, dict]]:
    """nanopolish / f5c ``call-methylation`` table.

    ``start`` and ``end`` are 0-based plus-strand C positions of the first and
    last CpG of the call for both read strands, so ``start + 1`` / ``end + 1``
    are the 1-based atlas positions. A call with ``num_motifs > 1`` is one
    measurement over consecutive CpGs; the atlas lookup averages over them.
    """
    usecols = ["read_name", "chromosome", "start", "end", "log_lik_ratio"]
    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=chunksize,
        dtype={"read_name": "string", "chromosome": "category", "start": np.int64,
               "end": np.int64, "log_lik_ratio": np.float64},
    )
    for chunk in reader:
        calls = pd.DataFrame(
            {
                "read_id": chunk["read_name"].to_numpy(),
                "chrom": chunk["chromosome"].astype(str).to_numpy(),
                "start": chunk["start"].to_numpy(np.int64) + 1,
                "end": chunk["end"].to_numpy(np.int64) + 1,
                "r": _expit(tau * chunk["log_lik_ratio"].to_numpy(np.float64)),
            }
        )
        stats = {"n_rows": int(len(chunk)), "n_non_primary": 0, "n_other_mod": 0,
                 "n_calls": int(len(calls))}
        yield calls, stats


FORMATS = {"modkit": iter_modkit_calls, "nanopolish": iter_nanopolish_calls}
