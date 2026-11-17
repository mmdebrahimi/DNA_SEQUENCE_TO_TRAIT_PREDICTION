"""Tests for the active-frame re-grade (v2).

The whole point of the re-grade is that it is auditable as a ONE-FIELD change to a frozen
pre-registration. These tests pin that: the original verdict survives, the inputs are the recorded
ones, exactly one bar field moves and it is the null gate, the effect-size threshold is untouched,
and a differently-shaped artifact is refused rather than silently graded.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.active_frame_regrade import RegradeRefused, load_arms, recorded_ci, regrade
from scripts.glm_active_frame_gate import PREREGISTERED_ACTIVE_FRAME

GATE = Path("wiki/glm_active_frame_gate_2026-11-16.json")

pytestmark = pytest.mark.skipif(
    not GATE.exists(), reason="the recorded gate artifact is absent"
)


@pytest.fixture(scope="module")
def artifact() -> dict:
    return json.loads(GATE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def result(artifact: dict) -> dict:
    return regrade(artifact)


# --------------------------------------------------------------------------------------
# the one-field claim, which is the whole argument
# --------------------------------------------------------------------------------------
def test_exactly_one_field_of_the_frozen_bar_changed_and_it_is_the_null_gate(
    result: dict,
) -> None:
    diff = result["bar_diff_against_the_frozen_preregistration"]
    assert result["bar_diff_field_count"] == 1, diff
    assert set(diff) == {"null_max_abs"}
    assert diff["null_max_abs"]["frozen"] == PREREGISTERED_ACTIVE_FRAME["null_max_abs"]


def test_the_effect_size_threshold_is_untouched(result: dict) -> None:
    """A noise-derived bar cannot evaluate an effect-size threshold. v1 of this re-grade replaced
    +0.15 anyway; this pins that it no longer does."""
    diff = result["bar_diff_against_the_frozen_preregistration"]
    assert "explains_margin" not in diff
    assert "explains_min_seeds" not in diff
    assert "contributes_margin" not in diff
    assert result["why_only_the_null_gate_needed_changing"]["explains_margin_bar"] == 0.15


def test_the_withdrawn_claim_is_recorded_not_quietly_dropped(result: dict) -> None:
    """The 'correction runs in both directions' argument was invalid and reached the user. An
    artifact that silently stopped making it would leave the retraction unfindable."""
    w = result["a_claim_WITHDRAWN"]
    assert "WITHDRAWN" in w
    assert "both" in w.lower() and "direction" in w.lower()
    assert "EFFECT-SIZE" in w or "effect-size" in w


def test_the_effect_size_criterion_was_already_met_on_the_recorded_data(result: dict) -> None:
    """This is why only one gate needed changing."""
    w = result["why_only_the_null_gate_needed_changing"]
    assert w["effect_size_criterion_was_already_met"] is True
    assert w["median_margin"] >= w["explains_margin_bar"]
    assert w["seeds_clearing_the_UNCHANGED_effect_size_bar"] >= w["explains_min_seeds"]


# --------------------------------------------------------------------------------------
# provenance: nothing retrained, nothing recomputed
# --------------------------------------------------------------------------------------
def test_the_original_verdict_is_preserved_verbatim(artifact: dict, result: dict) -> None:
    assert result["original_verdict_PRESERVED"]["verdict"] == artifact["verdict"]
    assert result["original_verdict_PRESERVED"]["verdict"] == "INDETERMINATE_NULL_NOT_CLEAN"
    assert result["original_verdict_PRESERVED"]["detail"] == artifact["verdict_detail"]


def test_the_recorded_ci_is_read_and_not_recomputed(artifact: dict) -> None:
    """A recomputed bootstrap CI is a different number from the one the frozen ladder was applied
    to, so the re-grade must read it."""
    assert recorded_ci(artifact) == tuple(artifact["phase2_active"]["margin_ci"])


def test_a_missing_ci_is_refused_rather_than_recomputed(artifact: dict) -> None:
    stripped = copy.deepcopy(artifact)
    stripped["phase2_active"].pop("margin_ci", None)
    with pytest.raises(RegradeRefused, match="paired-bootstrap CI"):
        regrade(stripped)


# --------------------------------------------------------------------------------------
# the null arm
# --------------------------------------------------------------------------------------
def test_the_headline_is_that_a_clean_pipeline_fails_the_frozen_null_bar(result: dict) -> None:
    cfr = result["null_arm"]["frozen_bar_clean_failure_rate"]
    assert cfr > 0.85, f"the whole argument rests on this being near-certain; got {cfr}"


def test_every_frozen_bar_breach_is_in_the_harmless_direction(result: dict) -> None:
    na = result["null_arm"]
    assert na["n_breaches_of_the_frozen_bar"] == na["n_breaches_that_are_negative"]
    assert na["n_positive_nulls_above_the_frozen_bar"] == 0


def test_the_replacement_null_bar_is_the_CONSERVATIVE_one(result: dict) -> None:
    """deff=1 must be TIGHTER than the measured-ICC bar, else choosing it proves nothing."""
    na = result["null_arm"]
    assert na["replacement_bar_used"] < na["looser_bar_at_the_measured_icc_NOT_used"]
    assert result["independence_unit"]["used_for_the_bar"] is False


def test_three_named_mis_specifications_are_recorded(result: dict) -> None:
    reasons = result["null_arm"]["three_independent_mis_specifications"]
    assert len(reasons) == 3
    joined = " ".join(reasons).lower()
    assert "max" in joined and "two-sided" in joined and "independent units" in joined


# --------------------------------------------------------------------------------------
# the paired empirical null
# --------------------------------------------------------------------------------------
def test_the_paired_null_is_arithmetically_right_and_all_positive(
    artifact: dict, result: dict
) -> None:
    ps = artifact["phase2_active"]["per_seed"]
    seeds = sorted(ps, key=int)
    expected = [
        round(float(ps[s]["encoder_multi"]) - float(ps[s]["encoder_shuffled"]), 4)
        for s in seeds
    ]
    pe = result["paired_empirical_null"]
    assert pe["per_seed"] == expected
    assert pe["all_positive"] is True
    assert pe["min"] > 0


def test_the_pairing_rationale_names_the_v1_error(result: dict) -> None:
    """v1 compared each margin to the UNPAIRED null maximum and reported a misleading 1.07x."""
    pe = result["paired_empirical_null"]
    assert "unpaired" in pe["why_paired"].lower()
    assert "1.07" in pe["why_paired"]


def test_no_p_value_is_claimed_for_the_paired_null(result: dict) -> None:
    """The seeds share peaks, so a sign test would overstate independence."""
    pe = result["paired_empirical_null"]
    assert "overstate" in pe["no_p_value_on_purpose"].lower()
    blob = json.dumps(result).lower()
    assert "p_value" not in blob.replace("no_p_value_on_purpose", "")
    assert "p =" not in blob and "p<" not in blob


# --------------------------------------------------------------------------------------
# non-vacuity: the re-grade must still be able to refuse
# --------------------------------------------------------------------------------------
def test_a_weak_run_does_not_reach_the_effect_size_verdict(artifact: dict) -> None:
    weak = copy.deepcopy(artifact)
    for s in weak["phase2_active"]["per_seed"]:
        weak["phase2_active"]["per_seed"][s]["margin"] = 0.01
    weak["phase2_active"]["margin_ci"] = [-0.01, 0.02]
    assert regrade(weak)["verdict"] != "FRAME_EXPLAINS_THE_GAP"


def test_a_dirty_null_still_blocks_the_grade(artifact: dict) -> None:
    """The null arm must retain its power to block, in the direction that matters -- UPWARD."""
    dirty = copy.deepcopy(artifact)
    for s in dirty["phase2_active"]["per_seed"]:
        dirty["phase2_active"]["per_seed"][s]["encoder_shuffled"] = 0.9
    assert regrade(dirty)["verdict"] == "INDETERMINATE_NULL_NOT_CLEAN"


def test_mismatched_held_out_peak_counts_are_refused(artifact: dict) -> None:
    odd = copy.deepcopy(artifact)
    seeds = sorted(odd["phase2_active"]["per_seed"], key=int)
    odd["phase2_active"]["per_seed"][seeds[0]]["held_out_peaks"] = 1
    with pytest.raises(RegradeRefused, match="held-out peak counts differ"):
        regrade(odd)


def test_an_artifact_without_per_seed_records_is_refused() -> None:
    with pytest.raises(RegradeRefused, match="nothing to re-grade"):
        load_arms({"phase2_active": {"per_seed": {}}})


def test_a_wrong_shaped_artifact_is_refused_rather_than_defaulted() -> None:
    with pytest.raises(RegradeRefused):
        load_arms({"not_the_phase": {}})


# --------------------------------------------------------------------------------------
# honesty surface
# --------------------------------------------------------------------------------------
def test_honest_limits_name_the_user_decision_and_the_weak_seed(result: dict) -> None:
    limits = " ".join(result["honest_limits"])
    assert "USER decision" in limits
    assert "weakest seed" in limits
    assert "0.941" in limits


def test_the_framing_recommendation_comes_from_the_frozen_map(result: dict) -> None:
    """It must be the pre-committed string for this verdict, not a fresh sentence."""
    assert result["framing_recommendation"] == PREREGISTERED_ACTIVE_FRAME[
        "framing_recommendation"
    ][result["verdict"]]
