"""A falsified hypothesis and a mis-scoped claim must both stay visible in the record."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "wiki" / "salmserovar_o7_allele_length_probe_2026-09-04.json"


@pytest.fixture(scope="module")
def d():
    if not ART.exists():
        pytest.skip("O7 probe artifact absent")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_length_hypothesis_is_recorded_as_falsified(d):
    assert d["verdict"] == "ALLELE_LENGTH_HYPOTHESIS_FALSIFIED"
    lo, hi = d["o_allele_length_range"]
    assert lo < d["o7_reference_length"] < hi, "O7 must sit inside the ordinary length range"


def test_the_concentration_claim_is_corrected_not_deleted(d):
    """The narrow claim was true of its subset; the general one was not. Both scopes stay visible."""
    by = d["what_the_pattern_actually_is"]["by_o_group"]
    assert by.get("3,10", 0) >= by.get("7", 0), (
        "the corrected finding is that another O group has at least as many partial hits")
    note = d["what_the_pattern_actually_is"]["note"]
    assert "not o7-exclusive" in note.lower() and "abstention" in note.lower()


def test_the_revised_reading_explains_why_coverage_was_the_right_lever(d):
    r = d["revised_reading"]
    assert "partial homology" in r
    assert "COVERAGE" in r and "identity" in r


def test_it_does_not_overclaim_a_replacement_cause(d):
    """Four causes about this cell were asserted and measured wrong; this must not add a fifth."""
    assert any("does not establish" in s.lower() or "not a closed question" in s.lower()
               for s in d["honest_limits"]), "the remaining unknown must be stated"
