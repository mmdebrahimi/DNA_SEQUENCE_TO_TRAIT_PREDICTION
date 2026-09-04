"""The adopted threshold must stay tied to the evidence that justified it.

The coverage cut was lowered on a pre-registered rule, not a hunch. These pin the rule, the deployed
value, and the asymmetry that makes it safe — so a later "tidy-up" cannot quietly revert or loosen it.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "wiki" / "salmserovar_threshold_tradeoff_2026-09-04.json"


def test_deployed_thresholds_match_the_adopted_setting():
    from dna_decode.salmserovar.runner import (
        SEROVAR_COVERAGE_THRESHOLD, SEROVAR_IDENTITY_THRESHOLD)
    assert SEROVAR_COVERAGE_THRESHOLD == 40.0, "the adopted coverage cut"
    assert SEROVAR_IDENTITY_THRESHOLD == 90.0, (
        "identity must NOT move -- it was never the failing axis; relaxing it would admit "
        "genuinely different alleles")


@pytest.mark.skipif(not ART.exists(), reason="tradeoff artifact absent")
def test_the_rule_was_registered_before_the_run(d=None):
    d = json.loads(ART.read_text(encoding="utf-8"))
    p = d["preregistered"]
    assert p["registered_before_run"] is True
    assert p["adopt_if"]["newly_wrong_max"] < p["adopt_if"]["rescued_correct_min"], (
        "the bar must be ASYMMETRIC -- a confident wrong serovar is not recoverable, an abstention is")
    assert "reject_if" in p, "a rule with no reject branch is not a rule"


@pytest.mark.skipif(not ART.exists(), reason="tradeoff artifact absent")
def test_the_adopted_result_actually_cleared_its_own_bar():
    d = json.loads(ART.read_text(encoding="utf-8"))
    p = d["preregistered"]["adopt_if"]
    assert d["verdict"] == "ADOPT"
    assert d["rescued_correct"] >= p["rescued_correct_min"]
    assert d["newly_wrong"] <= p["newly_wrong_max"]
    assert d["net"] >= p["net_gain_min"]


@pytest.mark.skipif(not ART.exists(), reason="tradeoff artifact absent")
def test_the_trade_was_measured_on_BOTH_selective_axes():
    """Coverage alone is not enough -- a caller that abstains less can simply be wrong more."""
    sc = json.loads(ART.read_text(encoding="utf-8"))["selective_classification"]
    assert sc["relaxed"]["coverage"] > sc["deployed"]["coverage"]
    assert sc["relaxed"]["accuracy_on_covered"] >= sc["deployed"]["accuracy_on_covered"], (
        "accuracy-on-covered must not have been sacrificed to buy coverage")


@pytest.mark.skipif(not ART.exists(), reason="tradeoff artifact absent")
def test_the_loss_side_is_recorded_not_just_the_gain():
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert "correct->wrong" in d["transition_matrix"] or d["newly_wrong"] == 0
    assert isinstance(d["newly_wrong_detail"], list)
    assert any("does NOT mean the O-antigen defect is unfixable" in s for s in d["honest_limits"])
