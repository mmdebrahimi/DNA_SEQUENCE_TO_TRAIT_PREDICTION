"""Guards for the oracle-selection measurement.

The expensive failure here is a FLATTERING metric, and this file was written after one near-miss and one
real defect: the smoke run printed a confident substantive verdict on a single axis against a bar needing
three, and a 0.10 margin bar turned out to be the same order as the measured noise.
"""
from __future__ import annotations

import random

import pytest

from dna_decode.glm.selection import (
    SelectionResult,
    assert_k1_is_degenerate,
    midrank_percentiles,
    permutation_null,
    simulate_selection,
)


def _pool(n=300, seed=0):
    r = random.Random(seed)
    return [r.random() for _ in range(n)]


# ---------------------------------------------------------------------------------------------------
# mid-ranks -- ties must not be resolved by index order
# ---------------------------------------------------------------------------------------------------
def test_midrank_ties_share_one_percentile():
    """Ordinal ranks hand a tie-block's whole spread to whichever element sorted first. This repo has
    already measured a case where that silently shifted a published median."""
    assert midrank_percentiles([5, 1, 5, 3]) == [pytest.approx(5 / 6), 0.0, pytest.approx(5 / 6),
                                                 pytest.approx(1 / 3)]


def test_midrank_all_tied_is_the_midpoint_not_an_ordering():
    assert midrank_percentiles([2, 2, 2]) == [0.5, 0.5, 0.5]


def test_midrank_degenerate_sizes():
    assert midrank_percentiles([]) == []
    assert midrank_percentiles([7.0]) == [0.5]


def test_midrank_is_index_order_invariant():
    """Reversing the input must permute the output identically, never change the values."""
    vals = [3.0, 1.0, 3.0, 2.0, 1.0]
    assert midrank_percentiles(vals)[::-1] == midrank_percentiles(vals[::-1])


# ---------------------------------------------------------------------------------------------------
# the three arms
# ---------------------------------------------------------------------------------------------------
def test_a_perfect_scorer_captures_exactly_all_the_headroom():
    """The upper anchor. Scoring with the truth itself must give headroom 1.0 exactly -- if it does not,
    the ceiling arm is not the ceiling."""
    m = _pool()
    r = simulate_selection(m, m, k=10, trials=500)
    assert r.headroom_captured == pytest.approx(1.0)
    assert r.mean_pct_picked == pytest.approx(r.mean_pct_ceiling)


def test_an_anticorrelated_scorer_captures_NEGATIVE_headroom():
    """The lower anchor. Without it, a metric that cannot go below zero would look reassuring on junk."""
    m = _pool()
    r = simulate_selection([-v for v in m], m, k=10, trials=500)
    assert r.headroom_captured < -0.8
    assert r.win_rate_vs_random == 0.0


def test_a_useless_scorer_is_centred_on_ZERO_across_realisations():
    """THE load-bearing null check, and the measurement that motivated `permutation_null`.

    A single useless scorer does NOT give headroom ~0 -- it gives one draw from a distribution whose spread
    `trials` cannot reduce (measured stdev ~0.1 at n=300, k=10). Only averaging over SCORER REALISATIONS
    recovers the zero. Getting this backwards is how a +0.13 reading gets mistaken for a harness bias.
    """
    m = _pool()
    hs = []
    for s in range(40):
        rr = random.Random(5000 + s)
        junk = [rr.random() for _ in range(len(m))]
        hs.append(simulate_selection(junk, m, k=10, trials=300, seed=s).headroom_captured)
    mean = sum(hs) / len(hs)
    assert abs(mean) < 0.05, f"useless scorers must centre on 0, got {mean}"
    assert any(h < 0 for h in hs), "a correct null must go negative sometimes"
    assert any(h > 0 for h in hs)


def test_k1_is_degenerate_by_construction_and_headroom_refuses_to_be_a_number():
    """With one candidate the pick IS the random draw IS the best. A lift at k=1 localises a bug to the
    harness, so this is the simulation's free self-test."""
    m = _pool()
    r = simulate_selection(_pool(seed=9), m, k=1, trials=400)
    assert r.headroom_captured is None, "0/0 must be None, never a fabricated number"
    assert_k1_is_degenerate(r)


def test_assert_k1_catches_a_broken_harness():
    """Non-vacuity: the self-test must FAIL on arms that differ."""
    bad = SelectionResult(scorer="x", k=1, n_trials=1, n_pool=10,
                          mean_pct_picked=0.9, mean_pct_random=0.5, mean_pct_ceiling=0.9,
                          median_measured_picked=1, median_measured_random=1, median_measured_ceiling=1,
                          fold_over_random=1.0, headroom_captured=None,
                          win_rate_vs_random=0.0, tie_rate_vs_random=0.0)
    with pytest.raises(AssertionError, match="must be degenerate"):
        assert_k1_is_degenerate(bad)
    with pytest.raises(ValueError, match="only meaningful at k=1"):
        assert_k1_is_degenerate(SelectionResult(
            scorer="x", k=2, n_trials=1, n_pool=10, mean_pct_picked=0.5, mean_pct_random=0.5,
            mean_pct_ceiling=0.5, median_measured_picked=1, median_measured_random=1,
            median_measured_ceiling=1, fold_over_random=1.0, headroom_captured=None,
            win_rate_vs_random=0.0, tie_rate_vs_random=0.0))


def test_headroom_rises_with_k_for_a_real_scorer():
    """More candidates means more available headroom AND more for a working scorer to capture. This is the
    propose-k / score-k / keep-best shape the inverse cell also measured."""
    m = _pool()
    r = random.Random(3)
    noisy = [v + r.gauss(0, 0.25) for v in m]        # a decent but imperfect scorer
    hs = [simulate_selection(noisy, m, k=k, trials=1500, seed=1).headroom_captured for k in (2, 5, 20)]
    assert hs[0] < hs[2], f"headroom should rise with k, got {hs}"


# ---------------------------------------------------------------------------------------------------
# ties in the SCORER must break randomly
# ---------------------------------------------------------------------------------------------------
def test_an_all_tied_scorer_is_indistinguishable_from_random():
    """A scorer emitting one value for everything must capture ~0, not whatever index order gave it.

    GC content over 150 bp takes few distinct values, so heavy ties are real here, and `sorted()`-order
    tie-breaking has already shifted a median elsewhere in this repo.
    """
    m = _pool()
    r = simulate_selection([1.0] * len(m), m, k=10, trials=3000)
    assert abs(r.headroom_captured) < 0.06, f"tied scorer must not gain, got {r.headroom_captured}"


# ---------------------------------------------------------------------------------------------------
# the null band
# ---------------------------------------------------------------------------------------------------
def test_permutation_null_is_centred_on_zero_and_reports_a_real_spread():
    nb = permutation_null(_pool(), k=10, trials=200, n_perms=30)
    assert nb["n_perms"] == 30
    assert abs(nb["mean"]) < 0.06
    assert nb["stdev"] > 0.0, "a null with no spread would make the band decorative"
    assert nb["min"] < 0.0 < nb["max"]


def test_permutation_null_refuses_to_invent_a_band_at_k1():
    """At k=1 every headroom is None, so there is no band. It must say so rather than return zeros."""
    nb = permutation_null(_pool(), k=1, trials=50, n_perms=5)
    assert nb["n_perms"] == 0
    assert "no null band" in nb["note"]


def test_permutation_null_is_scorer_independent_by_construction():
    """The null is a property of (pool, k) alone -- a shuffled scorer is a uniform permutation whatever it
    was shuffled from. That is WHY the runner computes it once per cell instead of once per scorer."""
    m = _pool()
    a = permutation_null(m, k=5, trials=150, n_perms=20, seed=4)
    b = permutation_null(m, k=5, trials=150, n_perms=20, seed=4)
    assert a == b


# ---------------------------------------------------------------------------------------------------
# refusals
# ---------------------------------------------------------------------------------------------------
def test_refuses_misaligned_inputs():
    with pytest.raises(ValueError, match="must align"):
        simulate_selection([1.0, 2.0], [1.0, 2.0, 3.0], k=2)


def test_refuses_k_larger_than_the_pool():
    """Drawing k without replacement from fewer than k is not a smaller sample, it is a different
    experiment -- so it raises rather than silently degrading."""
    with pytest.raises(ValueError, match="exceeds pool size"):
        simulate_selection([1.0, 2.0], [1.0, 2.0], k=3)


def test_refuses_nonsense_k_and_trials():
    with pytest.raises(ValueError, match="k must be >= 1"):
        simulate_selection([1.0], [1.0], k=0)
    with pytest.raises(ValueError, match="trials must be >= 1"):
        simulate_selection([1.0], [1.0], k=1, trials=0)


def test_refuses_percentiles_that_do_not_align():
    with pytest.raises(ValueError, match="percentiles must align"):
        simulate_selection([1.0, 2.0], [1.0, 2.0], k=1, percentiles=[0.5])


def test_seed_is_deterministic():
    m = _pool()
    s = _pool(seed=11)
    a = simulate_selection(s, m, k=5, trials=200, seed=42).as_dict()
    b = simulate_selection(s, m, k=5, trials=200, seed=42).as_dict()
    assert a == b


def test_as_dict_keeps_none_headroom_as_none():
    """A rounding helper that turned None into 0.0 would convert 'undefined' into 'no effect'."""
    m = _pool()
    d = simulate_selection(_pool(seed=2), m, k=1, trials=50).as_dict()
    assert d["headroom_captured"] is None
