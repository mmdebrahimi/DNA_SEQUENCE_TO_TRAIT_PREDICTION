"""Guards for kill-count reporting (plan Step 3, 2026-10-01).

Encodes the two historical retractions as named regression cases: a filter that removed nothing must
never be reportable as an active control.
"""

from __future__ import annotations

import pytest

from dna_decode.constraints import report as RP
from dna_decode.constraints.registry import INAPPLICABLE, REFUTED, SATISFIED, ConstraintVerdict


def _v(kind):
    return ConstraintVerdict(kind)


def test_the_three_states_are_distinct_and_derived():
    r = RP.ConstraintReport()
    r.record("refutes", _v(REFUTED))
    for _ in range(5):
        r.record("removes_nothing", _v(SATISFIED))
    r.note_not_evaluated("never_ran")

    assert r.status("refutes") == RP.ACTIVE
    assert r.status("removes_nothing") == RP.INERT
    assert r.status("never_ran") == RP.NOT_EVALUATED
    assert len({RP.ACTIVE, RP.INERT, RP.NOT_EVALUATED}) == 3


def test_inert_is_NOT_conflated_with_satisfied_or_with_active():
    """The whole point: 5 satisfied cases and 0 refutations is 'inert', not 'working'."""
    r = RP.ConstraintReport()
    for _ in range(5):
        r.record("k", _v(SATISFIED))
    assert r.counts["k"].n_satisfied == 5
    assert r.counts["k"].n_refuted == 0
    assert r.status("k") == RP.INERT
    assert r.status("k") != RP.ACTIVE


def test_inapplicable_does_not_count_as_evaluated():
    """"could not check" must not become "checked and found nothing"."""
    r = RP.ConstraintReport()
    for _ in range(3):
        r.record("k", _v(INAPPLICABLE))
    assert r.counts["k"].n_inapplicable == 3
    assert r.counts["k"].n_evaluated == 0
    assert r.status("k") == RP.NOT_EVALUATED


def test_a_silent_inert_constraint_RAISES_and_the_guard_is_non_vacuous():
    r = RP.ConstraintReport()
    r.record("inert_one", _v(SATISFIED))
    with pytest.raises(RP.SilentInertConstraintError, match="not a control"):
        r.assert_no_silent_inert()
    r.acknowledge_inert("inert_one", "curated input already satisfies this law")
    r.assert_no_silent_inert()                 # passes once acknowledged -- both directions proven


def test_an_acknowledgement_requires_a_reason():
    r = RP.ConstraintReport()
    r.record("k", _v(SATISFIED))
    with pytest.raises(ValueError):
        r.acknowledge_inert("k", "")


def test_an_active_constraint_needs_no_acknowledgement():
    r = RP.ConstraintReport()
    r.record("k", _v(REFUTED))
    r.assert_no_silent_inert()


def test_regression_the_hiv_censoring_guard_shape():
    """HISTORICAL CASE 1. `hiv_nnrti_mutant_catalog` filtered operator-prefixed folds (`<`/`>`) when
    there are zero such values in any Stanford column. Under this report that filter is `inert` and
    cannot be presented as a control without an explicit acknowledgement."""
    r = RP.ConstraintReport()
    for _ in range(2168):                      # the real EFV cohort size; zero censored values
        r.record("fold_operator_filter", _v(SATISFIED))
    assert r.status("fold_operator_filter") == RP.INERT
    with pytest.raises(RP.SilentInertConstraintError):
        r.assert_no_silent_inert()


def test_regression_the_resfinder_point_row_shape():
    """HISTORICAL CASE 2. The ResFinder comparison excluded AMRFinder `POINT` rows as a control when
    there are zero POINT rows across all 1,818 runs. Same shape, same refusal."""
    r = RP.ConstraintReport()
    for _ in range(1818):
        r.record("point_row_exclusion", _v(SATISFIED))
    assert r.counts["point_row_exclusion"].n_refuted == 0
    assert r.status("point_row_exclusion") == RP.INERT


def test_the_artifact_makes_an_inert_constraint_VISIBLE():
    """A block carried only in memory is not a disclosure -- inert constraints must appear in the dict
    with their status and reason, not be absent from it."""
    r = RP.ConstraintReport()
    r.record("a", _v(REFUTED))
    r.record("b", _v(SATISFIED))
    r.note_not_evaluated("c")
    r.acknowledge_inert("b", "reference data is already clean")
    d = r.as_dict()
    assert set(d["constraints"]) == {"a", "b", "c"}
    assert d["constraints"]["b"]["status"] == RP.INERT
    assert d["constraints"]["b"]["acknowledged_reason"] == "reference data is already clean"
    assert d["constraints"]["c"]["status"] == RP.NOT_EVALUATED
    assert d["summary"] == {"n_active": 1, "n_inert": 1, "n_not_evaluated": 1}
    assert "not evidence of protection" in d["reading"].lower().replace("--", "")


def test_counts_accumulate_across_many_records():
    r = RP.ConstraintReport()
    for kind, n in ((REFUTED, 3), (SATISFIED, 7), (INAPPLICABLE, 2)):
        for _ in range(n):
            r.record("k", _v(kind))
    c = r.counts["k"]
    assert (c.n_refuted, c.n_satisfied, c.n_inapplicable, c.n_evaluated) == (3, 7, 2, 10)
    assert r.status("k") == RP.ACTIVE
