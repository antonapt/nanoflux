"""Weighted pileup: counting, thresholds, strand combining and the infer-compatible layout."""

import numpy as np
import pandas as pd
import pytest

from src.scoring.pileup import _SiteIndex, load_annotation, load_weights, weighted_pileup
from src.tests.test_scoring import modkit_row, write_modkit
from src.utils.data_utils import load_pipeline_data


def write_weights(path, weights: dict):
    pd.DataFrame({"read_id": list(weights), "weight": list(weights.values())}).to_parquet(path)


@pytest.fixture
def annotation(tmp_path):
    path = tmp_path / "anno.bed"
    pd.DataFrame([["chr1", 100, 102, "cgA"], ["chr1", 200, 202, "cgB"], ["chr2", 100, 102, "cgC"]]).to_csv(
        path, sep="\t", header=False, index=False)
    return path


def test_site_index_overlap_rule(annotation):
    sites = _SiteIndex(load_annotation(annotation))
    # plus-strand C at 100 -> CpG [100, 102) overlaps probe [100, 102); 99 and 101 overlap too, 98/102 do not
    assert sites.lookup(np.array(["chr1"] * 5), np.array([98, 99, 100, 101, 102])).tolist() == [-1, 0, 0, 0, -1]
    assert sites.lookup(np.array(["chr2", "chr3"]), np.array([100, 100])).tolist() == [2, -1]


def test_weighted_pileup_counts_and_layout(tmp_path, annotation):
    rows = [
        # cgA (C at 0-based 100): three reads
        modkit_row("r1", 100, "+", 0.95, "m"),   # weight 1.0, methylated
        modkit_row("r2", 101, "-", 0.90, "-"),   # minus strand: G at 101 -> same CpG; weight 0.25, canonical
        modkit_row("r3", 100, "+", 0.60, "m"),   # below threshold 0.8 -> fail
        # cgB (C at 200): one unscored read
        modkit_row("r4", 200, "+", 0.99, "-"),
        # not on a probe
        modkit_row("r1", 150, "+", 0.99, "m"),
        # secondary alignment and 5hmC are ignored
        modkit_row("r5", 100, "+", 0.99, "m", flag=256),
        modkit_row("r5", 100, "+", 0.99, "h"),
    ]
    write_modkit(tmp_path / "reads.tsv", rows)
    write_weights(tmp_path / "w.parquet", {"r1": 1.0, "r2": 0.25, "r3": 1.0})

    out = tmp_path / "methylation.bed"
    stats = weighted_pileup(tmp_path / "reads.tsv", annotation, load_weights(tmp_path / "w.parquet"),
                            out, filter_threshold=0.8, unscored_weight=0.25)
    assert stats["n_calls_on_probes"] == 4 and stats["n_probes_covered"] == 2
    assert stats["n_reads_unscored"] == 1 and stats["filter_threshold"] == 0.8

    bed = pd.read_csv(out, sep="\t", header=None)
    assert bed.shape == (2, 14)
    a = bed[bed[13] == "cgA"].iloc[0]
    valid, pct, nmod, ncan, _, _, nfail, _, _ = a[9].split(" ")
    assert float(valid) == 1.25 and float(nmod) == 1.0 and float(ncan) == 0.25 and float(nfail) == 1.0
    assert float(pct) == pytest.approx(80.0, abs=0.01)
    assert (a[0], a[1], a[2]) == ("chr1", 100, 101)
    b = bed[bed[13] == "cgB"].iloc[0]
    assert b[9].split(" ")[0] == "0.25" and b[9].split(" ")[1] == "0.00"

    # infer's loader reads it: percent -> -1 / 1 / NaN per probe
    X = load_pipeline_data(out)
    assert X.loc["cgA", "mod"] == pytest.approx(1.0) and X.loc["cgB", "mod"] == -1.0


def test_unit_weights_equal_plain_counts_and_auto_threshold(tmp_path, annotation):
    rng = np.random.default_rng(0)
    rows, truth = [], {"cgA": [0, 0], "cgB": [0, 0]}
    probs = rng.uniform(0.5, 1.0, 400)
    thr = np.quantile(probs, 0.1)
    for i, p in enumerate(probs):
        probe, pos = ("cgA", 100) if i % 2 == 0 else ("cgB", 200)
        meth = rng.random() < 0.7
        strand = "+" if rng.random() < 0.5 else "-"
        rows.append(modkit_row(f"r{i}", pos if strand == "+" else pos + 1, strand, p, "m" if meth else "-"))
        if p >= thr:
            truth[probe][0] += 1
            truth[probe][1] += meth
    write_modkit(tmp_path / "reads.tsv", rows)
    write_weights(tmp_path / "w.parquet", {f"r{i}": 1.0 for i in range(400)})
    stats = weighted_pileup(tmp_path / "reads.tsv", annotation, load_weights(tmp_path / "w.parquet"),
                            tmp_path / "out.bed")
    assert stats["filter_threshold"] == pytest.approx(thr, abs=2 / 1024)  # histogram resolution
    bed = pd.read_csv(tmp_path / "out.bed", sep="\t", header=None)
    for probe, (n, k) in truth.items():
        valid, pct, *_ = bed[bed[13] == probe].iloc[0][9].split(" ")
        assert float(valid) == pytest.approx(n, abs=1) and float(pct) == pytest.approx(100 * k / n, abs=1.0)


def test_zero_weight_reads_drop_sites(tmp_path, annotation):
    write_modkit(tmp_path / "reads.tsv", [modkit_row("r1", 100, "+", 0.99, "m")])
    write_weights(tmp_path / "w.parquet", {"r1": 0.0})
    stats = weighted_pileup(tmp_path / "reads.tsv", annotation, load_weights(tmp_path / "w.parquet"),
                            tmp_path / "out.bed", filter_threshold=0.5)
    assert stats["n_probes_covered"] == 0
    assert (tmp_path / "out.bed").stat().st_size == 0
