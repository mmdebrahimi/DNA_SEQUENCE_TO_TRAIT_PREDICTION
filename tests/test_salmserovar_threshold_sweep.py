"""The sweep's job was to check the deployed cut; it also caught its own selection rule misfiring.

These pin the thing that is easy to lose later: the deployed value survived a held-out check, and the
pre-registered rule that would have changed it was rejected on evidence rather than on preference.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "wiki" / "salmserovar_threshold_sweep_2026-09-04.json"


@pytest.fixture(scope="module")
def d():
    if not ART.exists():
        pytest.skip("sweep artifact absent")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_deployed_cut_did_not_change(d):
    """The sweep's outcome was KEEP, so the runner must still carry 40."""
    from dna_decode.salmserovar.runner import SEROVAR_COVERAGE_THRESHOLD
    assert SEROVAR_COVERAGE_THRESHOLD == 40.0
    assert d["verdict"] == "SELECTION_RULE_MISFIRED_KEEP_DEPLOYED_40"


def test_the_choice_was_made_on_a_half_that_did_not_report_the_number(d):
    p = d["preregistered"]
    assert p["registered_before_sweep"] is True
    assert "SELECT" in p["split"] and "CONFIRM" in p["reported_number_comes_from"]
    assert d["split"]["select_n"] > 0 and d["split"]["confirm_n"] > 0


def test_the_rule_was_not_argmax(d):
    """Argmax picks the luckiest cut; the rule deliberately did not."""
    assert "NOT argmax" in d["preregistered"]["selection_rule"]


def test_the_misfire_is_recorded_with_its_cause_and_not_quietly_dropped(d):
    m = d["selection_rule_MISFIRED"]
    assert m["what_the_rule_picked"] != 40.0, "the rule genuinely disagreed with the deployed value"
    assert "PREMISE" in m["why_the_rule_was_wrong"]
    assert "KEEP" in m["action"]
    # the reason it is not rule-breaking must be stated, or this reads as convenient
    assert "CHOOSE on SELECT, CHECK on CONFIRM" in m["why_this_is_NOT_rule_breaking"]
    assert "already deployed" in m["why_this_is_NOT_rule_breaking"]


def test_the_heldout_evidence_actually_dominates_on_both_axes(d):
    """'Dominates' is a strong word -- it must mean more correct AND no more wrong."""
    c = d["confirm_half"]
    assert c["deployed_40"]["correct"] > c["chosen"]["correct"]
    assert c["deployed_40"]["wrong"] <= c["chosen"]["wrong"]


def test_the_dose_response_that_falsified_the_premise_is_visible(d):
    """The premise died because the response is graded, not noisy -- that must be checkable."""
    sweep = {s["coverage_cut"]: s["correct"] for s in d["select_half"]["sweep"]}
    base = d["select_half"]["baseline_80"]["correct"]
    assert sweep[70.0] > base and sweep[60.0] > sweep[70.0] and sweep[40.0] > sweep[50.0]
    assert sweep[30.0] < sweep[40.0], "the turning point is what makes 40 defensible"


def test_optimality_is_not_claimed(d):
    """The debt this run closed is 'defensible', NOT 'optimal' -- the artifact must not blur them."""
    disclaimers = [s for s in d["honest_limits"] if "optimal" in s.lower()]
    assert disclaimers, "no honest_limit mentions optimality at all"
    assert any(("not establish" in s.lower() or "not shown" in s.lower() or "might find" in s.lower())
               for s in disclaimers), f"optimality is mentioned but not disclaimed: {disclaimers}"
    assert any("defensible" in s.lower() for s in d["honest_limits"]), (
        "the weaker claim actually supported must be stated too")
