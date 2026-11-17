"""Tests for the active-frame re-grade.

The point of these is that the re-grade must be auditable as a BAR change and nothing else: the
original verdict survives, the inputs are the recorded ones, and a differently-shaped artifact is
refused rather than silently graded.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.active_frame_regrade import RegradeRefused, load_arms, regrade

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


def test_the_original_verdict_is_preserved_verbatim(artifact: dict, result: dict) -> None:
    """A re-grade that overwrote the first verdict would be indistinguishable from having only
    ever run the second."""
    assert result["original_verdict_PRESERVED"]["verdict"] == artifact["verdict"]
    assert result["original_verdict_PRESERVED"]["verdict"] == "INDETERMINATE_NULL_NOT_CLEAN"
    assert result["original_verdict_PRESERVED"]["detail"] == artifact["verdict_detail"]


def test_the_margins_are_read_and_not_recomputed(artifact: dict, result: dict) -> None:
    recorded = [
        float(artifact["phase2_active"]["per_seed"][s]["margin"])
        for s in sorted(artifact["phase2_active"]["per_seed"], key=int)
    ]
    assert result["margin_arm"]["per_seed_margins"] == recorded


def test_the_headline_is_that_a_clean_pipeline_fails_the_frozen_null_bar(result: dict) -> None:
    cfr = result["null_arm"]["frozen_bar_validation"]["clean_failure_rate"]
    assert cfr > 0.85, f"the whole argument rests on this being near-certain; got {cfr}"


def test_every_frozen_bar_breach_is_in_the_harmless_direction(result: dict) -> None:
    """Leakage can only push a null upward, so a negative breach is not evidence of the hazard."""
    na = result["null_arm"]
    assert na["n_nulls_breaching_frozen_bar"] == na["n_breaches_that_are_negative"]
    assert na["n_positive_nulls_above_frozen_bar"] == 0


def test_the_correction_runs_in_both_directions(result: dict) -> None:
    """The null bar was too TIGHT and the margin bar too CONSERVATIVE. A bar-shopping exercise
    would not surrender a finding that one of its own thresholds was needlessly strict."""
    assert result["null_arm"]["frozen_bar"] < result["null_arm"][
        "derived_bar_one_sided_hazard_matched"
    ]
    assert result["margin_arm"]["frozen_bar"] > result["margin_arm"][
        "derived_bar_5pct_clean_failure_PROXY"
    ]
    # negligible, NOT exactly zero: this is a Monte Carlo rate over 40,000 trials, so asserting
    # equality to 0.0 pins sampling luck rather than the finding.
    assert result["margin_arm"]["frozen_bar_validation"]["clean_failure_rate"] < 0.005


def test_the_design_effect_is_well_above_one(result: dict) -> None:
    u = result["independence_unit"]
    assert u["design_effect"] > 5.0
    assert u["n_effective"] < u["n_rows"] / 5.0
    assert u["measured_from"]


def test_the_verdict_is_graded_on_this_data(result: dict) -> None:
    assert result["verdict"] == "FRAME_EXPLAINS_THE_GAP"
    assert result["margin_arm"]["seeds_clearing_derived_bar"] >= 8


def test_the_verdict_is_not_vacuous_a_weak_run_stays_indeterminate(artifact: dict) -> None:
    """Proof the re-grade can still refuse: shrink every margin below the derived bar and the
    verdict must change. Without this, FRAME_EXPLAINS_THE_GAP could be unconditional."""
    weak = copy.deepcopy(artifact)
    for s in weak["phase2_active"]["per_seed"]:
        weak["phase2_active"]["per_seed"][s]["margin"] = 0.01
    assert regrade(weak)["verdict"] == "INDETERMINATE_UNDER_DERIVED_BAR"


def test_a_dirty_null_still_blocks_the_grade(artifact: dict) -> None:
    """The null arm must retain its power to block, in the direction that matters."""
    dirty = copy.deepcopy(artifact)
    for s in dirty["phase2_active"]["per_seed"]:
        dirty["phase2_active"]["per_seed"][s]["encoder_shuffled"] = 0.9
    assert regrade(dirty)["verdict"] == "INDETERMINATE_UNDER_DERIVED_BAR"


def test_mismatched_held_out_peak_counts_are_refused(artifact: dict) -> None:
    """One independence unit cannot describe seeds that held out different numbers of peaks."""
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


def test_the_primary_margin_null_is_empirical_and_arithmetically_right(
    artifact: dict, result: dict
) -> None:
    """The primary null must be the recorded shuffled-label difference, not a simulated one, so
    it carries no variance model and no design effect. Re-derive it from the artifact."""
    ps = artifact["phase2_active"]["per_seed"]
    seeds = sorted(ps, key=int)
    expected = [
        round(float(ps[s]["encoder_shuffled"]) - float(ps[s]["gc"]), 4) for s in seeds
    ]
    emp = result["margin_arm"]["PRIMARY_empirical_null"]
    assert emp["per_seed_margin_nulls"] == expected
    assert emp["margin_null_max"] == pytest.approx(max(expected), abs=1e-4)
    assert emp["seeds_exceeding_the_empirical_null_max"] >= 8


def test_the_thinness_of_the_weakest_seed_is_reported_not_smoothed(result: dict) -> None:
    """The weakest seed clears the empirical null by a hair. If that ever stops being stated the
    artifact has started flattering itself."""
    emp = result["margin_arm"]["PRIMARY_empirical_null"]
    assert emp["weakest_margin_over_null_max"] < 1.5
    assert "thin" in emp["honest_reading"].lower()
    assert any("weakest seed clears" in lim for lim in result["honest_limits"])


def test_the_verdict_does_not_depend_on_the_icc_modelling_choice(result: dict) -> None:
    """The deff=1 probe is the tightest bar the design can produce. If the verdict survives it,
    the most contestable input to the re-grade is not load-bearing."""
    rb = result["robustness_does_the_verdict_need_the_ICC_at_all"]
    assert rb["verdict_holds_at_deff_1"] is True
    # the deff=1 bars must genuinely be TIGHTER, else the probe is not conservative
    assert rb["deff_1_null_bar_one_sided"] < result["null_arm"][
        "derived_bar_one_sided_hazard_matched"
    ]


def test_the_simulated_margin_bars_are_labelled_a_proxy(result: dict) -> None:
    """`difference-of-spearman` has no variance model, so the simulated figures must not read as
    primary."""
    assert "PROXY" in result["margin_arm"]["simulated_bars_are_a_PROXY_not_primary"].upper()
    assert "derived_bar_5pct_clean_failure_PROXY" in result["margin_arm"]


def test_honest_limits_are_present_and_name_the_user_decision(result: dict) -> None:
    limits = " ".join(result["honest_limits"]).lower()
    assert "user decision" in limits
    assert "approximation" in limits or "asymptotic" in limits
    assert "0.941" in " ".join(result["honest_limits"])
