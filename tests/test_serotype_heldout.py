"""The held-out confirmation must stay honest about shrinkage and about the O-axis cost.

A replication that quietly reports the discovery number, or hides a small regression on the other
axis, would be worse than no replication. Both are pinned.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "wiki" / "serotype_heldout_confirm_2026-09-04.json"


@pytest.fixture(scope="module")
def d():
    if not ART.exists():
        pytest.skip("held-out artifact absent")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_prediction_was_registered_before_the_run(d):
    assert d["preregistered"]["registered_before_run"] is True
    assert d["preregistered"]["quantitative_bar"] == 0.05
    assert "falsified_if" in d["preregistered"], "a prediction with no failure condition is not one"


def test_held_out_gain_is_reported_and_is_SMALLER_than_discovery(d):
    """The unblinded discovery number must not be the one carried forward."""
    r = d["regression_from_discovery"]
    assert d["H_gain"] < r["discovery_H_gain"], "held-out gain should shrink vs an unblinded discovery"
    assert r["shrinkage"] > 0
    assert "quote the held-out" in r["interpretation"].lower()


def test_the_O_axis_regression_is_recorded_not_buried(d):
    """identity-primary is very slightly worse on O; that trade-off must stay visible."""
    assert d["O_gain"] < 0, "the artifact should record the real O-axis cost"
    ne = d["net_effect"]
    assert ne["O_misses_introduced"] > 0
    assert ne["net_calls_corrected"] == ne["H_misses_fixed"] - ne["O_misses_introduced"]


def test_the_mechanism_prediction_was_checked_not_just_the_effect(d):
    """Resolution must be unchanged -- the rule picks WHICH allele, not WHETHER one is picked."""
    s = d["secondary_prediction_outcome"]
    assert s["met"] is True and s["observed_no_call_delta"] == 0.0


def test_disjointness_was_measured_not_assumed(d):
    h = d["heldout"]
    assert h["n_excluded_as_discovery_overlap"] > 0, "overlap must be measured and removed"
    assert h["seed"] != 23, "must use a different seed from the discovery cohort"
    assert h["cohort_meta"]["passes_source_diversity_bar"] is True


def test_lineage_limit_is_stated(d):
    assert any("not by lineage" in s.lower() or "lineage" in s.lower() for s in d["honest_limits"])
