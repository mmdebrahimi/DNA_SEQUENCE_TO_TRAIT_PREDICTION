"""Tests for the pool-rule diagnosis probe's pure logic (`scripts/fba_pool_rule_probe.py`).

Offline and substrate-free: both functions under test are pure. The model-dependent arms are not
tested here — they need iML1515 and are exercised by the probe's own run.

The properties worth pinning are the three that already produced wrong conclusions during the probe's
own development, each of which looked like a finding rather than a bug:

  * a NaN deleted-max means the deletion made the LP INFEASIBLE, which is the STRONGEST form of
    load-bearing — and `nan <= tol` is False, so a two-state test files it as irrelevant;
  * `collapsed` alone is too strict on a genome-scale model, so a graded `reduced` state must exist
    without counting as pool-excluding;
  * the set comparison must report the DIRECTION of a disagreement, because "the rule excluded things
    sign detection never flagged" and "sign detection flagged things the rule kept" are different
    defects with different fixes.
"""
from __future__ import annotations

import pytest

from scripts.fba_pool_rule_probe import (
    LOAD_BEARING_STATES,
    classify_flux_role,
    compare_sets,
)

NAN = float("nan")
TOL = 1e-6


# --------------------------------------------------------------- classify_flux_role


def test_nan_deleted_max_is_infeasible_not_unaffected():
    """The bug that made the probe report an empty load-bearing set for all five targets."""
    assert classify_flux_role(10.0, NAN, TOL) == "infeasible"


def test_none_deleted_max_is_treated_like_nan():
    assert classify_flux_role(10.0, None, TOL) == "infeasible"


def test_infeasible_counts_as_load_bearing():
    """An infeasible deletion solve is genuine essentiality in this repo, not a solver failure."""
    assert "infeasible" in LOAD_BEARING_STATES


def test_a_naive_threshold_would_misfile_nan_and_that_is_why_the_state_exists():
    """Pins the language trap directly: the comparison a two-state test would have used is False."""
    assert (NAN <= TOL) is False
    assert classify_flux_role(10.0, NAN, TOL) in LOAD_BEARING_STATES


@pytest.mark.parametrize("deleted", [0.0, 1e-9, 1e-7])
def test_at_or_under_tolerance_is_collapsed(deleted):
    assert classify_flux_role(10.0, deleted, TOL) == "collapsed"


def test_collapsed_counts_as_load_bearing():
    assert "collapsed" in LOAD_BEARING_STATES


def test_a_large_reduction_is_reduced_not_collapsed():
    """The measured LDH_D shape: lactate max 7.256 -> 0.806, a 9x drop that is not zero."""
    assert classify_flux_role(7.256, 0.806, TOL) == "reduced"


def test_reduced_does_NOT_count_as_load_bearing():
    """Graded evidence is reported, but it must not silently exclude a route from the pool."""
    assert "reduced" not in LOAD_BEARING_STATES


def test_a_small_change_is_unaffected():
    assert classify_flux_role(10.0, 9.9, TOL) == "unaffected"


def test_a_deletion_that_RAISES_the_ceiling_is_unaffected():
    """The measured PFL shape: deleting it lifts the formate/lactate ceiling rather than lowering it."""
    assert classify_flux_role(7.256, 25.66, TOL) == "unaffected"


def test_an_unknown_intact_max_cannot_manufacture_a_reduced_verdict():
    """Without a baseline there is nothing to be reduced relative to."""
    assert classify_flux_role(NAN, 5.0, TOL) == "unaffected"
    assert classify_flux_role(None, 5.0, TOL) == "unaffected"


def test_a_zero_intact_max_cannot_manufacture_a_reduced_verdict():
    """intact_max <= tol means the product was never achievable; a ratio would divide into noise."""
    assert classify_flux_role(0.0, 5.0, TOL) == "unaffected"


def test_every_state_is_one_of_the_four_declared_values():
    four = {"infeasible", "collapsed", "reduced", "unaffected"}
    for intact in (None, NAN, 0.0, 10.0):
        for deleted in (None, NAN, 0.0, 1e-9, 0.8, 9.9, 25.0):
            assert classify_flux_role(intact, deleted, TOL) in four


# -------------------------------------------------------------------- compare_sets


def test_identical_sets_agree_and_have_no_difference():
    r = compare_sets(["A", "B"], ["B", "A"])
    assert r["sets_agree"] is True
    assert r["symmetric_difference"] == []


def test_the_measured_acetate_shape_reports_the_excluded_not_detected_direction():
    """The rule excluded ACKr/PTAr/PTA2; sign detection flags only POX for ac_c."""
    r = compare_sets(["POX"], ["ACKr", "PTAr", "PTA2", "POX"])
    assert r["sets_agree"] is False
    assert r["excluded_but_NOT_sign_detected"] == ["ACKr", "PTA2", "PTAr"]
    assert r["sign_detected_but_NOT_excluded"] == []
    assert "DISAGREE" in r["reading"]


def test_the_opposite_direction_is_reported_separately():
    """Two different defects with two different fixes must not collapse into one field."""
    r = compare_sets(["X", "Y"], ["X"])
    assert r["sign_detected_but_NOT_excluded"] == ["Y"]
    assert r["excluded_but_NOT_sign_detected"] == []


def test_the_measured_formate_shape_agrees_so_breadth_is_the_only_charge():
    r = compare_sets(["PFL"], ["PFL"])
    assert r["sets_agree"] is True
    assert "not the problem" in r["reading"]


def test_empty_inputs_agree_rather_than_raising():
    r = compare_sets([], [])
    assert r["sets_agree"] is True
    assert r["n_sign_detected"] == 0 and r["n_rule_excluded"] == 0


def test_duplicates_do_not_inflate_the_counts():
    r = compare_sets(["A", "A", "B"], ["B", "B"])
    assert r["n_sign_detected"] == 2
    assert r["n_rule_excluded"] == 1


def test_symmetric_difference_is_the_union_of_both_directions():
    r = compare_sets(["A", "B"], ["B", "C"])
    assert set(r["symmetric_difference"]) == {"A", "C"}
    assert set(r["excluded_but_NOT_sign_detected"]) | set(r["sign_detected_but_NOT_excluded"]) == {"A", "C"}
