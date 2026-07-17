"""Tests for the multi-edit verdict (the handoff's second experiment)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.forward_inverse_multiedit import DMS_ID, _spearman

ART = Path(__file__).resolve().parent.parent / "wiki" / "forward_inverse_multiedit_2026-07-17.json"


def test_spearman_is_a_rank_correlation():
    assert _spearman([1, 2, 3], [10, 20, 30]) == pytest.approx(1.0)
    assert _spearman([1, 2, 3], [30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_ignores_monotone_rescaling():
    """It must not be fooled by a scale change -- the whole point of using rank here."""
    assert _spearman([1, 2, 3], [1, 100, 10000]) == pytest.approx(1.0)


def test_spearman_of_a_constant_is_zero_not_nan():
    assert _spearman([1, 2, 3], [5, 5, 5]) == 0.0


@pytest.mark.skipif(not ART.exists(), reason="multiedit artifact not present")
def test_the_extension_is_recorded_as_not_warranted_on_its_own_motivation():
    """The handoff's motivation was 'for targets NO SINGLE EDIT REACHES'. That is a factual claim about the
    data, and it is false where checkable: 0/694 real doubles land outside the single-edit range."""
    r = json.loads(ART.read_text(encoding="utf-8"))
    assert r["Q1_reach"]["verdict"] == "NOT_MOTIVATED"
    assert r["Q1_reach"]["n_multi_outside_single_range"] == 0
    assert r["verdict"] == "MULTIEDIT_EXTENSION_NOT_WARRANTED"


@pytest.mark.skipif(not ART.exists(), reason="multiedit artifact not present")
def test_additivity_ranks_well_but_does_not_pin_the_magnitude():
    """The third independent sighting of tonight's theme. Additive composition ranks doubles strongly
    (Spearman ~0.80) yet its absolute error is WORSE than predicting the assay median -- summing two
    effects overshoots a bounded assay's floor. Rank and magnitude are separate capabilities."""
    q2 = json.loads(ART.read_text(encoding="utf-8"))["Q2_additivity"]
    assert q2["spearman_additive_vs_measured"] > 0.7          # ranks well...
    assert q2["additive_beats_null"] is False                  # ...but does not pin the dose


@pytest.mark.skipif(not ART.exists(), reason="multiedit artifact not present")
def test_the_artifact_states_its_n_equals_1_substrate_limit():
    """One assay cannot prove multi-edits NEVER extend reach; the honest claim is narrower."""
    r = json.loads(ART.read_text(encoding="utf-8"))
    assert "n=1 substrate" in r["honest_scope"]
    assert r["substrate"]["dms_id"] == DMS_ID
    assert r["substrate"]["orders_present"] == [2]             # doubles only -- no triples to speak for
