"""Tests for the pure position-resolved reductions.

Every case here is answerable BY CONSTRUCTION on a synthetic tensor -- no checkpoint, no D:, no
forward pass -- which is the whole reason the reductions were split out of the scripts.

The two load-bearing properties are the parity inverse (`global_mean_from_bin_means` must reproduce
`z.mean(axis=2)`, because the scoring run's parity gate is exactly that arithmetic over STORED
features) and the non-finite refusals (a NaN that reached the scorer would be mapped to 0.0 and read
as a clean null rather than as a broken reduction).
"""

from __future__ import annotations

import numpy as np
import pytest

from dna_decode.glm.positional_readout import (
    PositionalReadoutError,
    bin_maxes,
    bin_means,
    channel_autocorrelation,
    flatten_bins,
    global_mean_from_bin_means,
    orient,
)

# --------------------------------------------------------------------------------------------------
# bin_means / bin_maxes -- answers known by construction
# --------------------------------------------------------------------------------------------------


def _staircase(n_b=2, n_c=3, per_bin=2, n_bins=3):
    """z whose per-bin means are exact integers: bin j of row (b, c) holds `offset + j + 1` twice."""
    z = np.zeros((n_b, n_c, per_bin * n_bins))
    expect = np.zeros((n_b, n_c, n_bins))
    for b in range(n_b):
        for c in range(n_c):
            off = 10 * b + c
            for j in range(n_bins):
                v = off + j + 1
                # straddle the mean so the bin mean is exactly v but no element equals it
                z[b, c, j * per_bin:(j + 1) * per_bin] = [v - 1, v + 1]
                expect[b, c, j] = v
    return z, expect


def test_bin_means_returns_the_constructed_bin_means_exactly():
    z, expect = _staircase()
    got = bin_means(z, 3)
    assert got.shape == expect.shape
    assert np.array_equal(got, expect), "bin means are exact integers here; no tolerance is needed"


def test_bin_maxes_returns_the_constructed_bin_maxes_exactly():
    z, expect = _staircase()
    # each bin holds [v-1, v+1], so the max is v+1
    assert np.array_equal(bin_maxes(z, 3), expect + 1)


def test_one_bin_is_the_global_mean_and_L_bins_is_the_identity():
    rng = np.random.default_rng(11)
    z = rng.normal(size=(4, 5, 12))
    assert np.allclose(bin_means(z, 1)[:, :, 0], z.mean(axis=2))
    assert np.allclose(bin_means(z, 12), z)  # bin width 1


# --------------------------------------------------------------------------------------------------
# the parity inverse -- this IS the scoring run's Gate B arithmetic
# --------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("n_bins", [1, 2, 3, 5, 6, 10, 15, 25, 30, 50, 75, 150])
def test_global_mean_from_bin_means_reproduces_the_mean_pool_at_every_dividing_rung(n_bins):
    """The real geometry: 150 positions, 128 channels. Every divisor must recover the mean-pool."""
    rng = np.random.default_rng(2026)
    z = rng.normal(size=(7, 128, 150))
    recovered = global_mean_from_bin_means(bin_means(z, n_bins))
    assert recovered.shape == (7, 128)
    assert np.allclose(recovered, z.mean(axis=2), rtol=0, atol=1e-12)


def test_the_parity_inverse_is_not_vacuous_a_perturbed_block_fails_it():
    """If a perturbed block still passed, Gate B would be theatre."""
    rng = np.random.default_rng(3)
    z = rng.normal(size=(5, 8, 10))
    bm = bin_means(z, 5)
    bm[2, 4, 1] += 0.5
    assert not np.allclose(global_mean_from_bin_means(bm), z.mean(axis=2), atol=1e-9)


# --------------------------------------------------------------------------------------------------
# the ragged refusal -- an unequal final bin would break the inverse quietly
# --------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("n_bins", [4, 7, 8, 9, 11, 100])
@pytest.mark.parametrize("fn", [bin_means, bin_maxes])
def test_non_dividing_n_bins_is_refused_rather_than_truncated_or_padded(fn, n_bins):
    z = np.zeros((1, 2, 150))
    with pytest.raises(PositionalReadoutError, match="does not divide"):
        fn(z, n_bins)


@pytest.mark.parametrize("fn", [bin_means, bin_maxes])
def test_n_bins_out_of_range_is_refused(fn):
    z = np.zeros((1, 2, 10))
    with pytest.raises(PositionalReadoutError, match="exceeds the position axis"):
        fn(z, 11)
    with pytest.raises(PositionalReadoutError, match=r"n_bins must be >= 1"):
        fn(z, 0)


# --------------------------------------------------------------------------------------------------
# channel_autocorrelation -- peak at the period, trough at half it
# --------------------------------------------------------------------------------------------------

LAGS = list(range(1, 21))


def test_a_pure_sinusoid_peaks_at_its_period_and_troughs_at_half_of_it():
    period = 10
    t = np.arange(150)
    z = np.sin(2 * np.pi * t / period)[None, None, :]
    r = channel_autocorrelation(z, LAGS)[0]
    assert r.shape == (len(LAGS),)
    assert r[LAGS.index(period)] > 0.999, "r at the period must saturate"
    assert r[LAGS.index(period // 2)] < -0.999, "r at half the period must anti-saturate"
    # the per-lag form has no length attenuation, so harmonics saturate too: the argmax lag is a
    # MULTIPLE of the period, not necessarily the period itself. This is a feature block, not a
    # period estimator.
    assert LAGS[int(np.argmax(r))] % period == 0
    assert LAGS[int(np.argmin(r))] % period == period // 2


@pytest.mark.parametrize("period", [5, 6, 10, 15, 20])
def test_the_period_peak_generalises_across_periods(period):
    """The property the biased sample ACF FAILS at long periods: see the module docstring."""
    t = np.arange(150)
    z = np.cos(2 * np.pi * t / period)[None, None, :]
    r = channel_autocorrelation(z, LAGS)[0]
    assert r[LAGS.index(period)] > 0.999, "no length attenuation may eat the period peak"
    assert LAGS[int(np.argmax(r))] % period == 0, "argmax must be a multiple of the period"
    assert r.min() < -0.75, "a periodic channel must anti-correlate somewhere in range"


def test_the_estimator_carries_no_length_attenuation():
    """Pin the exact defect that the biased sample ACF had, so a revert cannot pass unnoticed.

    Under the biased form a period-15 cosine on L=150 scored r_15 = 0.9000 < r_1 = 0.9014 -- the
    peak moved to lag 1 purely because of the (L-k)/L shrinkage.
    """
    t = np.arange(150)
    z = np.cos(2 * np.pi * t / 15)[None, None, :]
    r = channel_autocorrelation(z, LAGS)[0]
    assert r[LAGS.index(15)] > r[LAGS.index(1)]
    # and a constant-amplitude signal must score the same at a long lag as at a short one
    z2 = np.cos(2 * np.pi * t / 2)[None, None, :]  # period 2: every even lag is a full cycle
    r2 = channel_autocorrelation(z2, [2, 20])[0]
    assert abs(r2[0] - r2[1]) < 1e-9, "lag 20 must not be shrunk relative to lag 2"


def test_minimum_overlap_is_enforced_because_a_two_point_pearson_is_always_saturated():
    z = np.arange(10, dtype=float).reshape(1, 1, 10)
    with pytest.raises(PositionalReadoutError, match="MIN_OVERLAP"):
        channel_autocorrelation(z, [8])  # overlap 2
    assert channel_autocorrelation(z, [7]).shape == (1, 1)  # overlap 3 is allowed


def test_white_noise_gives_near_zero_at_every_non_zero_lag():
    """Averaged over 64 independent channels the per-lag ACF concentrates hard on 0."""
    rng = np.random.default_rng(7)
    z = rng.normal(size=(1, 64, 150))
    r = channel_autocorrelation(z, LAGS).reshape(64, len(LAGS))
    per_lag_mean = r.mean(axis=0)
    assert np.abs(per_lag_mean).max() < 0.06, "noise must not fake periodic structure"
    assert np.abs(r).mean() < 0.12


def test_autocorrelation_is_bounded_and_translation_invariant():
    """Anchor-free is the property that justifies this readout on an un-registered window."""
    rng = np.random.default_rng(5)
    base = rng.normal(size=200)
    a = base[:150][None, None, :]
    b = base[7:157][None, None, :]  # same series, shifted start
    ra = channel_autocorrelation(a, LAGS)[0]
    rb = channel_autocorrelation(b, LAGS)[0]
    assert np.abs(ra).max() <= 1.0 + 1e-12
    # not bit-identical (different windows of the series), but close -- a bin readout would differ
    # far more under the same shift
    assert np.abs(ra - rb).max() < 0.5


def test_a_constant_channel_yields_zero_not_nan():
    """Zero variance means the ACF is undefined; 0.0 is the honest encoding, NaN is not."""
    z = np.zeros((1, 2, 50))
    z[0, 1, :] = np.arange(50)  # one live channel beside the constant one
    r = channel_autocorrelation(z, [1, 2, 3]).reshape(2, 3)
    assert np.all(np.isfinite(r))
    assert np.array_equal(r[0], np.zeros(3))
    assert np.abs(r[1]).sum() > 0, "the live channel must not have been zeroed too"


def test_lag_zero_is_refused_because_it_is_a_constant_column():
    z = np.zeros((1, 1, 10))
    with pytest.raises(PositionalReadoutError, match="non-positive lag"):
        channel_autocorrelation(z, [0, 1, 2])


def test_negative_duplicate_empty_and_oversized_lag_sets_are_refused():
    z = np.zeros((1, 1, 10))
    with pytest.raises(PositionalReadoutError, match="non-positive lag"):
        channel_autocorrelation(z, [-1, 2])
    with pytest.raises(PositionalReadoutError, match="duplicate lag"):
        channel_autocorrelation(z, [3, 3, 4])
    with pytest.raises(PositionalReadoutError, match="is empty"):
        channel_autocorrelation(z, [])
    with pytest.raises(PositionalReadoutError, match=r">= the position axis"):
        channel_autocorrelation(z, [9, 10])


def test_autocorrelation_column_order_is_channel_major():
    """Must match `flatten_bins` so the two feature blocks concatenate predictably."""
    rng = np.random.default_rng(13)
    z = rng.normal(size=(3, 4, 60))
    lags = [1, 5, 9]
    flat = channel_autocorrelation(z, lags)
    assert flat.shape == (3, 4 * len(lags))
    for c in range(4):
        single = channel_autocorrelation(z[:, c:c + 1, :], lags)
        assert np.allclose(flat[:, c * len(lags):(c + 1) * len(lags)], single)


# --------------------------------------------------------------------------------------------------
# orient -- an involution, and a no-op on an all-plus cohort
# --------------------------------------------------------------------------------------------------


def test_orient_applied_twice_is_the_identity():
    rng = np.random.default_rng(17)
    bm = rng.normal(size=(6, 5, 10))
    strands = ["+", "-", "-", "+", "-", "+"]
    assert np.array_equal(orient(orient(bm, strands), strands), bm)


def test_orient_is_a_no_op_on_an_all_plus_cohort():
    rng = np.random.default_rng(19)
    bm = rng.normal(size=(4, 3, 8))
    assert np.array_equal(orient(bm, ["+"] * 4), bm)


def test_orient_actually_reverses_the_minus_rows_and_leaves_plus_rows_alone():
    """Non-vacuity: an `orient` that did nothing would pass the involution test too."""
    bm = np.arange(2 * 1 * 4, dtype=float).reshape(2, 1, 4)
    out = orient(bm, ["+", "-"])
    assert np.array_equal(out[0], bm[0])
    assert np.array_equal(out[1, 0], bm[1, 0][::-1])
    assert not np.array_equal(out, bm)


def test_orient_reverses_only_the_bin_axis_never_the_channel_axis():
    """The channels are learned latents; there is no known complement permutation over them."""
    bm = np.zeros((1, 3, 4))
    bm[0, 0, 0] = 1.0
    out = orient(bm, ["-"])
    assert out[0, 0, 3] == 1.0, "the 1 should move along the BIN axis"
    assert out[0, 2, 3] == 0.0, "and must not have moved to another channel"


def test_orient_refuses_a_strand_length_mismatch():
    bm = np.zeros((3, 2, 4))
    with pytest.raises(PositionalReadoutError, match="strands for"):
        orient(bm, ["+", "-"])


@pytest.mark.parametrize("token", ["", ".", "plus", "+1", "*", None])
def test_orient_refuses_an_unrecognised_strand_token(token):
    bm = np.zeros((2, 2, 4))
    with pytest.raises(PositionalReadoutError, match="unrecognised strand"):
        orient(bm, ["+", token])


# --------------------------------------------------------------------------------------------------
# flatten_bins
# --------------------------------------------------------------------------------------------------


def test_flatten_bins_is_channel_major_and_round_trips():
    rng = np.random.default_rng(23)
    bm = rng.normal(size=(5, 7, 10))
    flat = flatten_bins(bm)
    assert flat.shape == (5, 70)
    assert np.array_equal(flat.reshape(5, 7, 10), bm)


def test_flatten_bins_feature_count_matches_the_bar_rungs():
    """128 channels x 10 bins = 1,280, the figure the plan and the bar both name."""
    z = np.zeros((2, 128, 150))
    assert flatten_bins(bin_means(z, 10)).shape[1] == 1280
    assert flatten_bins(bin_maxes(z, 10)).shape[1] == 1280
    assert flatten_bins(bin_means(z, 2)).shape[1] == 256
    assert flatten_bins(bin_means(z, 5)).shape[1] == 640


# --------------------------------------------------------------------------------------------------
# non-finite refusals at EVERY entry point -- the clean-null trap
# --------------------------------------------------------------------------------------------------

ENTRY_POINTS = [
    ("bin_means", lambda z: bin_means(z, 5)),
    ("bin_maxes", lambda z: bin_maxes(z, 5)),
    ("global_mean_from_bin_means", lambda z: global_mean_from_bin_means(z)),
    ("channel_autocorrelation", lambda z: channel_autocorrelation(z, [1, 2])),
    ("orient", lambda z: orient(z, ["+"] * z.shape[0])),
    ("flatten_bins", lambda z: flatten_bins(z)),
]


@pytest.mark.parametrize("name,fn", ENTRY_POINTS, ids=[n for n, _ in ENTRY_POINTS])
@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_every_entry_point_refuses_a_non_finite_input(name, fn, bad):
    z = np.zeros((2, 3, 10))
    z[1, 2, 7] = bad
    with pytest.raises(PositionalReadoutError, match="non-finite"):
        fn(z)


@pytest.mark.parametrize("name,fn", ENTRY_POINTS, ids=[n for n, _ in ENTRY_POINTS])
def test_every_entry_point_refuses_a_wrong_rank_input(name, fn):
    with pytest.raises(PositionalReadoutError, match="must be 3-D"):
        fn(np.zeros((4, 10)))


@pytest.mark.parametrize("name,fn", ENTRY_POINTS, ids=[n for n, _ in ENTRY_POINTS])
def test_every_entry_point_refuses_an_empty_axis(name, fn):
    with pytest.raises(PositionalReadoutError, match="empty axis"):
        fn(np.zeros((0, 3, 10)))
