"""Unit and end-to-end tests for nanoflux score on synthetic data."""

from argparse import Namespace

import numpy as np
import pandas as pd
import pytest

from src.scoring.atlas import Atlas
from src.scoring.calibration import N_CELLS, Calibration, cell_index
from src.scoring.calls import iter_modkit_calls
from src.scoring.prior import BetaPrior, MomentStats, PriorModel, fit_beta_prior
from src.scoring.scores import finalize_reads, score_calls
from src.scoring.weights import N_BIN_LABELS, bin_quantiles, compute_weights, parse_thresholds
from src.tasks import score as score_task  # noqa: F401  (lazy wrapper must import cleanly)
from src.tasks.score import OUTPUT_COLUMNS, main as score_main, parse_prior

MODKIT_COLUMNS = [
    "read_id", "forward_read_position", "ref_position", "chrom", "mod_strand", "ref_strand",
    "ref_mod_strand", "fw_soft_clipped_start", "fw_soft_clipped_end", "alignment_start",
    "alignment_end", "read_length", "call_prob", "call_code", "base_qual", "ref_kmer",
    "query_kmer", "canonical_base", "modified_primary_base", "fail", "inferred",
    "within_alignment", "flag",
]


def modkit_row(read, pos, strand, prob, code, flag=None):
    flag = (0 if strand == "+" else 16) if flag is None else flag
    return [read, 0, pos, "chr1", "+", strand, strand, 0, 0, 0, 0, 100, prob, code, 30,
            "ccgcc", "CCGCC", "C", "C", "false", "false", "true", flag]


def write_modkit(path, rows):
    pd.DataFrame(rows, columns=MODKIT_COLUMNS).to_csv(path, sep="\t", index=False)


# ------------------------------------------------------------------ atlas
def test_atlas_lookup_single_and_range():
    atlas = Atlas(np.array(["chr1"] * 4 + ["chr2"]), [10, 20, 30, 40, 10], [1, 2, 3, 4, 5], [10, 10, 10, 10, 10])
    m, n, count = atlas.lookup(["chr1", "chr1", "chr2", "chr3", "chr1"], [20, 20, 10, 10, 15], [20, 40, 10, 10, 35])
    assert count.tolist() == [1, 3, 1, 0, 2]
    assert m.tolist() == [2.0, 3.0, 5.0, 0.0, 2.5]
    assert n.tolist() == [10.0, 10.0, 10.0, 0.0, 10.0]


def test_atlas_min_cov_and_duplicates():
    atlas = Atlas(np.array(["chr1", "chr1"]), [10, 20], [1, 1], [1, 5], min_cov=3)
    assert atlas.n_sites == 1 and atlas.n_sites_total == 2
    with pytest.raises(ValueError, match="duplicate"):
        Atlas(np.array(["chr1", "chr1"]), [10, 10], [1, 1], [5, 5])


# ------------------------------------------------------------------ modkit adapter
def test_modkit_strand_mapping_and_filters(tmp_path):
    path = tmp_path / "calls.tsv"
    write_modkit(path, [
        modkit_row("a", 10, "+", 0.9, "m"),          # plus: C at 0-based 10 -> 1-based 11
        modkit_row("a", 11, "-", 0.8, "-"),          # minus: G at 0-based 11 == 1-based C 11
        modkit_row("b", 10, "+", 0.7, "h"),          # 5hmC dropped
        modkit_row("b", 10, "+", 0.7, "m", flag=256),  # secondary dropped
    ])
    (calls, stats), = list(iter_modkit_calls(path, chunksize=10))
    assert calls["start"].tolist() == [11, 11]
    assert calls["end"].tolist() == [11, 11]
    assert calls["r"].tolist() == pytest.approx([0.9, 0.2])
    assert stats == {"n_rows": 4, "n_non_primary": 1, "n_other_mod": 1, "n_calls": 2}


# ------------------------------------------------------------------ calibration
def test_cell_index_strata():
    m = np.array([0, 5, 9, 50, 100, 200])
    n = np.array([1, 5, 10, 100, 101, 1000])
    cells = cell_index(m, n)
    p_bins, cov_bins = cells // 6, cells % 6
    assert cov_bins.tolist() == [0, 0, 1, 4, 5, 5]
    assert p_bins.tolist() == [6, 17, 16, 10, 19, 4]  # floor(20 * (m+1)/(n+2))


def test_calibration_additive_and_fallback():
    rng = np.random.default_rng(0)
    n = np.full(1000, 10)
    m = rng.integers(0, 11, 1000)
    r = rng.random(1000)
    a, b = Calibration(min_calls=1), Calibration(min_calls=1)
    a.add(m[:500], n[:500], r[:500])
    b.add(m[500:], n[500:], r[500:])
    both = Calibration(min_calls=1)
    both.add(m, n, r)
    assert np.allclose((a + b).counts, both.counts) and np.allclose((a + b).sum_r2, both.sum_r2)

    # a sparse cell falls back to its p-bin, an empty p-bin to the global mean
    c = Calibration(min_calls=10)
    c.add(np.array([5] * 20), np.array([10] * 20), np.array([0.9] * 20))  # p-bin 10, cov 1
    c.add(np.array([50] * 2), np.array([100] * 2), np.array([0.1] * 2))   # same p-bin, cov 4: sparse
    _, e_r, _ = c.expected(np.array([50, 0]), np.array([100, 1]))
    assert e_r[0] == pytest.approx((20 * 0.9 + 2 * 0.1) / 22)  # p-bin pooled
    assert e_r[1] == pytest.approx((20 * 0.9 + 2 * 0.1) / 22)  # global
    assert len(c.tables()[0]) == N_CELLS


def test_calibration_roundtrip(tmp_path):
    c = Calibration()
    c.add(np.array([3, 8]), np.array([10, 10]), np.array([0.2, 0.95]))
    c.save(tmp_path / "c.json")
    d = Calibration.load(tmp_path / "c.json")
    assert d.n_calls == 2 and np.allclose(d.sum_r, c.sum_r)


# ------------------------------------------------------------------ scores
def test_reads_following_the_atlas_score_near_zero():
    """Reads drawn from the atlas get z centred at 0 with unit spread; flipped reads go negative."""
    rng = np.random.default_rng(1)
    n_reads, n_sites = 2000, 20
    p_true = rng.beta(0.5, 0.5, n_sites)
    n = np.full(n_sites, 200)
    m = rng.binomial(n, p_true)
    state = rng.random((n_reads, n_sites)) < p_true            # true methylation per read and site
    r = np.where(state, rng.uniform(0.8, 0.99, state.shape), rng.uniform(0.01, 0.2, state.shape))
    flipped = np.arange(n_reads) >= n_reads - 200                 # last 200 reads are discordant
    r[flipped] = 1 - r[flipped]

    calls = pd.DataFrame({"read_id": np.repeat(np.arange(n_reads).astype(str), n_sites),
                          "r": r.ravel()})
    mm, nn = np.tile(m, n_reads), np.tile(n, n_reads)
    calib = Calibration(min_calls=1)
    calib.add(mm[~np.repeat(flipped, n_sites)], nn[~np.repeat(flipped, n_sites)],
              calls["r"].to_numpy()[~np.repeat(flipped, n_sites)])
    reads = finalize_reads([score_calls(calls, mm, nn, np.ones(len(calls), dtype=int), calib)])

    z_ok = reads["z"].to_numpy()[~flipped]
    z_bad = reads["z"].to_numpy()[flipped]
    assert abs(z_ok.mean()) < 0.15 and 0.8 < z_ok.std() < 1.25
    assert np.median(z_bad) < -4
    assert (reads["n_cpg"] == n_sites).all() and (reads["n_bin"].astype(str) == "10+").all()


# ------------------------------------------------------------------ prior
def test_fit_beta_prior_recovers_parameters():
    rng = np.random.default_rng(3)
    p_true = rng.beta(1.5, 0.7, 20000)
    n = rng.integers(100, 400, 20000)
    m = rng.binomial(n, p_true)
    prior = fit_beta_prior(m, n, min_cov=100)
    assert prior.alpha == pytest.approx(1.5, rel=0.15)
    assert prior.beta == pytest.approx(0.7, rel=0.15)
    assert prior.p(np.array([3]), np.array([10]))[0] == pytest.approx((3 + prior.alpha) / (10 + prior.alpha + prior.beta))


def test_fit_beta_prior_lowers_coverage_when_sparse():
    n = np.full(3000, 20)
    m = np.random.default_rng(4).binomial(n, 0.7)
    prior = fit_beta_prior(m, n, min_cov=100, min_sites=1000)  # no site reaches 100
    assert prior.alpha > 0 and prior.beta > 0 and 0.6 < prior.mean < 0.8


def test_moment_fit_recovers_class_moments(tmp_path):
    rng = np.random.default_rng(5)
    p = rng.random(50000)
    state = rng.random(50000) < p
    r = np.where(state, rng.normal(0.85, 0.05, 50000), rng.normal(0.15, 0.05, 50000))
    stats = MomentStats()
    stats.add(p[:25000], r[:25000])
    stats.add(p[25000:], r[25000:])
    mo = stats.fit()
    assert mo.mu1 == pytest.approx(0.85, abs=0.02) and mo.mu0 == pytest.approx(0.15, abs=0.02)
    # the variances are second-moment differences and inherit the error of mu, so only loosely
    assert 0 <= mo.v1 < 0.02 and 0 <= mo.v0 < 0.02 and mo.n_calls == 50000

    model = PriorModel(BetaPrior(1.0, 1.0), mo)
    pp, e_r, v_r = model.expected(np.array([9, 0]), np.array([10, 10]))
    assert pp.tolist() == pytest.approx([10 / 12, 1 / 12])
    assert e_r[0] > e_r[1] and (v_r >= 0).all()
    model.save(tmp_path / "m.json")
    assert PriorModel.load(tmp_path / "m.json").to_dict() == model.to_dict()


def test_parse_prior():
    assert parse_prior(["auto"]) is None
    assert parse_prior(["1", "0.5"]) == BetaPrior(1.0, 0.5)
    for bad in (["1"], ["a", "b"], ["0", "1"], ["1", "2", "3"]):
        with pytest.raises(ValueError):
            parse_prior(bad)


# ------------------------------------------------------------------ weights
def test_parse_thresholds():
    t = parse_thresholds("1:-3, 2:-3.8,3-4:-4.6,5-9:-5.8,10+:-8")
    assert t == {"1": -3.0, "2": -3.8, "3-4": -4.6, "5-9": -5.8, "10+": -8.0}
    with pytest.raises(ValueError, match="missing"):
        parse_thresholds("1:-3")
    with pytest.raises(ValueError, match="unknown bin"):
        parse_thresholds("1:-3,2:-3,3-4:-3,5-9:-3,11:-3")


def test_weights_soft_and_hard():
    z = np.array([-10.0, -3.0, 0.5, 2.0, np.nan])
    n_bin = ["1"] * 5
    cutoffs, medians = {b: -3.0 for b in N_BIN_LABELS}, {b: 0.5 for b in N_BIN_LABELS}
    soft = compute_weights(z, n_bin, cutoffs, medians, mode="soft", temperature=1.0, w_min=0.1)
    scale = 0.5 - (-3.0)                                       # median - cutoff
    assert soft[0] == pytest.approx(0.1 + 0.9 / (1 + np.exp(-7.0 / scale)))  # two scales past cutoff
    assert soft[1] == pytest.approx(0.55)                      # cutoff: midway between floor and 1
    assert soft[2] == pytest.approx(0.1 + 0.9 / (1 + np.e))   # median: sigmoid(-1)
    assert soft[0] > soft[1] > soft[2] > soft[3] and np.isnan(soft[4])
    hard = compute_weights(z, n_bin, cutoffs, medians, mode="hard", w_min=0.0)
    assert hard[:4].tolist() == [1.0, 0.0, 0.0, 0.0] and np.isnan(hard[4])
    zero_temp = compute_weights(z, n_bin, cutoffs, medians, mode="soft", temperature=0.0)
    assert zero_temp[:4].tolist() == hard[:4].tolist()


def test_bin_quantiles():
    q = bin_quantiles(np.array([1.0, 2.0, 3.0, 10.0]), ["1", "1", "1", "2"], 0.5)
    assert q["1"] == 2.0 and q["2"] == 10.0 and np.isnan(q["10+"])


# ------------------------------------------------------------------ end to end
def test_score_command_end_to_end(tmp_path):
    rng = np.random.default_rng(2)
    positions = np.arange(100, 100 + 40 * 50, 50)               # 40 CpGs, 1-based C positions
    p_true = rng.beta(0.5, 0.5, len(positions))
    total = np.full(len(positions), 30)
    meth = rng.binomial(total, p_true)
    pd.DataFrame({"chrom_x": "chr1", "start_genomic": positions, "methylated_pooled": meth,
                  "total_pooled": total}).to_feather(tmp_path / "atlas.feather")

    rows = []
    for i in range(300):
        sites = rng.choice(len(positions), size=rng.integers(1, 15), replace=False)
        flip = i >= 270
        for s in sites:
            state = rng.random() < p_true[s]
            if flip:
                state = not state
            prob = rng.uniform(0.85, 0.99)
            strand = "+" if rng.random() < 0.5 else "-"
            pos0 = positions[s] - 1 if strand == "+" else positions[s]
            rows.append(modkit_row(f"read{i}", pos0, strand, prob, "m" if state else "-"))
    rows.append(modkit_row("unmatched", 5, "+", 0.9, "m"))     # no atlas CpG here
    write_modkit(tmp_path / "calls.tsv", rows)

    args = Namespace(
        input=tmp_path / "calls.tsv", output=tmp_path / "out", create_dir=True, debug=False,
        format="modkit", atlas=tmp_path / "atlas.feather",
        atlas_columns=["chrom_x", "start_genomic", "methylated_pooled", "total_pooled"],
        min_atlas_cov=1, tau=1.0, prior=["auto"], prior_min_cov=20, moments=None,
        calibration=None, fit_calibration=False, weight_mode="soft", weight_quantile=0.1,
        weight_thresholds=None, weight_temperature=1.0, weight_min=0.0, chunksize=1000,
    )
    score_main(args)

    out = pd.read_parquet(tmp_path / "out" / "read_scores.parquet")
    assert list(out.columns) == OUTPUT_COLUMNS
    assert len(out) == 300 and "unmatched" not in set(out["read_id"])
    assert out["weight"].between(0, 1).all()
    flipped = out["read_id"].str.slice(4).astype(int) >= 270
    assert out.loc[flipped, "weight"].mean() > out.loc[~flipped, "weight"].mean()
    assert (tmp_path / "out" / "moments.json").exists()
    assert not (tmp_path / "out" / "calibration.json").exists()
    summary = pd.read_json(tmp_path / "out" / "score_summary.json", typ="series")
    assert summary["site_model"] == "beta prior + moments"

    # reusing the moments file reproduces the scores
    args2 = Namespace(**{**vars(args), "output": tmp_path / "out2",
                         "moments": tmp_path / "out" / "moments.json"})
    score_main(args2)
    out2 = pd.read_parquet(tmp_path / "out2" / "read_scores.parquet")
    assert np.allclose(out["z"].to_numpy(), out2["z"].to_numpy(), equal_nan=True)

    # an explicit prior and explicit thresholds
    args3 = Namespace(**{**vars(args), "output": tmp_path / "out3", "prior": ["1", "1"],
                         "weight_thresholds": "1:-2,2:-2.5,3-4:-3,5-9:-3.5,10+:-4",
                         "weight_mode": "hard"})
    score_main(args3)
    out3 = pd.read_parquet(tmp_path / "out3" / "read_scores.parquet")
    assert set(out3["weight"].unique()) <= {0.0, 1.0}

    # calibration table: fit on this (out-of-atlas) sample, then reuse
    args4 = Namespace(**{**vars(args), "output": tmp_path / "out4", "fit_calibration": True})
    score_main(args4)
    assert (tmp_path / "out4" / "calibration.json").exists()
    args5 = Namespace(**{**vars(args), "output": tmp_path / "out5",
                         "calibration": tmp_path / "out4" / "calibration.json"})
    score_main(args5)
    out4 = pd.read_parquet(tmp_path / "out4" / "read_scores.parquet")
    out5 = pd.read_parquet(tmp_path / "out5" / "read_scores.parquet")
    assert np.allclose(out4["z"].to_numpy(), out5["z"].to_numpy(), equal_nan=True)
