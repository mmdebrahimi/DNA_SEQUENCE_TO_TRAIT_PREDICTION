"""Tests for the frozen frame-conditioning acceptance bar.

The bar's job is to make two specific post-hoc moves structurally unavailable: choosing the frame
granularity after seeing which one answers better, and grading absorption without the controls that
separate it from "the baseline simply got better". So the tests check that the primary granularity
is NAMED, that both controls are named and wired to the outcome they feed, and that the bar's
outcome set is EXACTLY the grader's check order -- a bar listing an outcome the grader cannot reach,
or vice versa, is how a dead branch survives review.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from dna_decode.glm.frame_verdict import CHECK_ORDER

REPO = pathlib.Path(__file__).resolve().parents[1]
BAR = REPO / "wiki/promoteratlas_frame_acceptance_bar.json"
PARENT = REPO / "wiki/promoteratlas_arm_acceptance_bar.json"

pytestmark = pytest.mark.skipif(not BAR.is_file(), reason="frame bar not frozen yet")


@pytest.fixture(scope="module")
def bar() -> dict:
    return json.loads(BAR.read_text(encoding="utf-8"))


# ------------------------------------------------------------------------------------- schema


def test_the_bar_declares_its_schema_and_freeze_date(bar):
    assert bar["schema"] == "promoteratlas-frame-acceptance-bar-v1"
    assert bar["frozen"]


def test_every_outcome_carries_a_meaning_and_a_CHECKABLE_predicate(bar):
    """Checkable means it relates MEASURED QUANTITIES, which is not the same as carrying a digit.

    Two legitimate predicate shapes exist here and each defeated a stricter version of this test.
    INDETERMINATE_BASELINE_NOT_STRENGTHENED compares two measured training fits to each other and
    correctly carries NO threshold constant -- demanding a digit would have forced a meaningless
    `- 0.0` into the bar. INDETERMINATE_HARNESS_DRIFT is a tolerance ("within 0.02") and carries no
    relational SYMBOL; the parent exposure bar words its drift gate the same way, so requiring an
    operator would have imposed a form the codebase does not use. The check is therefore a
    disjunction, and `test_the_checkability_check_is_non_vacuous` bounds the loosening -- a
    predicate that relates nothing and quantifies nothing still fails.
    """
    for name, body in bar["outcomes"].items():
        assert body.get("meaning"), name
        pred = body.get("predicate", "")
        assert pred, name
        assert _is_checkable(pred), (
            "%s's predicate neither relates nor quantifies anything, so it is a sentiment rather "
            "than a bar: %r" % (name, pred)
        )


def _is_checkable(pred: str) -> bool:
    """A predicate must either relate quantities or state a numeric tolerance."""
    return any(op in pred for op in ("<", ">")) or any(c.isdigit() for c in pred)


def test_the_checkability_check_is_non_vacuous():
    """Proof the loosened check can still fail: a sentiment passes neither limb."""
    assert not _is_checkable("the encoder looks exposed")
    assert not _is_checkable("the margin is substantially absorbed")
    assert _is_checkable("a < b")
    assert _is_checkable("within 0.02")


def test_every_margin_bearing_predicate_states_the_inherited_number(bar):
    """The thresholded outcomes must carry 0.02 explicitly rather than referring to it vaguely."""
    thresholded = [
        n for n in bar["outcomes"]
        if n in ("FRAME_EXPLAINS", "FRAME_PARTIAL", "FRAME_EXCLUDED_CONDITIONAL",
                 "INDETERMINATE_FRAME_UNINFORMATIVE",
                 "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE")
    ]
    assert len(thresholded) == 5
    for name in thresholded:
        assert "0.02" in bar["outcomes"][name]["predicate"], name


def test_the_check_order_lists_every_outcome_exactly_once(bar):
    order = bar["verdict_check_order"]
    assert len(order) == len(set(order)), "a duplicated outcome in the check order"
    assert set(order) == set(bar["outcomes"]), (
        "order/outcomes mismatch: %r" % (set(order) ^ set(bar["outcomes"]))
    )


def test_the_check_order_MATCHES_the_graders_check_order_exactly(bar):
    """Not just the same set -- the same sequence, because first-match-wins makes order semantic."""
    assert tuple(bar["verdict_check_order"]) == CHECK_ORDER


def test_every_indeterminate_is_ordered_before_every_substantive_outcome(bar):
    order = bar["verdict_check_order"]
    last_indet = max(i for i, n in enumerate(order) if n.startswith("INDETERMINATE_"))
    first_subst = min(i for i, n in enumerate(order) if not n.startswith("INDETERMINATE_"))
    assert last_indet < first_subst, "validity must be settled before any substantive claim"


def test_the_plumbing_gate_is_ordered_before_the_informativeness_gate(bar):
    """A baseline that never received the feature can itself make the frame look uninformative."""
    order = bar["verdict_check_order"]
    assert order.index("INDETERMINATE_BASELINE_NOT_STRENGTHENED") < order.index(
        "INDETERMINATE_FRAME_UNINFORMATIVE"
    )


# ------------------------------------------------------------------------- the frozen choices


def test_the_bar_NAMES_its_primary_granularity_so_the_choice_cannot_be_made_later(bar):
    assert bar["primary_frame_granularity"] == "in_window_any"
    assert bar["reported_also"] == "in_window_full"
    assert bar["primary_frame_granularity"] != bar["reported_also"]
    assert bar["in_window_any_threshold_bases"] == 1, "the threshold is named, not implied"


def test_the_bar_records_which_definition_of_FULL_it_means(bar):
    """Union coverage and inside-a-single-window are different predicates."""
    assert "UNION" in bar["granularity_rationale"]


def test_the_margin_is_inherited_from_the_parent_and_still_agrees_with_it(bar):
    assert bar["inherited_margin"] == 0.02
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    assert float(parent["margin"]["value"]) == bar["inherited_margin"], (
        "the frame bar has decoupled from the parent arm's margin"
    )
    assert "promoteratlas_arm_acceptance_bar.json" in bar["inherited_margin_source"]


def test_the_bar_refuses_to_be_a_new_arm(bar):
    """A positional or fine-tuned variant is candidate 8 under R6; this measures a LIMIT."""
    text = bar["not_a_new_arm"]
    assert "NEW ARM" in text and "R6" in text
    assert "must not move arms.native150" in text


# --------------------------------------------------------------------- the mandatory controls


def test_both_missing_controls_are_named_and_wired_to_the_outcome_they_feed(bar):
    ctrl = bar["controls_required_before_grading"]
    assert "matched_strength_series" in ctrl
    assert "size_matched_random_stratification" in ctrl
    for key in ("matched_strength_series", "size_matched_random_stratification"):
        assert ctrl[key]["feeds_outcome"] == "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE", (
            "%s names no outcome, so it could be reported and ignored" % key
        )
    assert ctrl["matched_strength_series"]["feeds_outcome"] in bar["outcomes"]


def test_the_matched_series_names_the_actual_comparators(bar):
    what = bar["controls_required_before_grading"]["matched_strength_series"]["what"]
    for comparator in ("gc", "dinucleotide", "kmer4", "frame"):
        assert comparator in what, comparator


def test_the_per_stratum_spread_is_required_so_range_restriction_is_distinguishable(bar):
    ctrl = bar["controls_required_before_grading"]["per_stratum_expression_spread"]
    assert "range restriction" in ctrl["why"]


# ------------------------------------------------------------------------- feasibility first


def test_feasibility_was_established_before_the_bar_was_frozen(bar):
    feas = bar["feasibility_established_first"]
    assert feas["windows_sorted"] == 3496
    assert feas["bases_covered_sorted"] == 625255
    assert feas["locate"]["n_unique"] == 43061
    assert feas["tiles_in_window_any"] == 28622
    assert feas["tiles_in_window_full"] == 8239
    assert feas["n_tiles_scored_by_the_arm"] == 44106


def test_no_frame_feature_was_constant_at_freeze_time(bar):
    """A constant feature could absorb nothing, so the experiment would be unanswerable."""
    var = bar["feasibility_established_first"]["frame_feature_variance"]
    assert var, "no variance recorded"
    for name, body in var.items():
        assert body["is_constant"] is False, name
        assert body["n_distinct_values"] >= 2, name


def test_the_bar_was_frozen_before_any_conditioning_number_existed(bar):
    before = " ".join(bar["frozen_before"]).lower()
    assert "ridge" in before and "margin" in before


def test_the_locate_floor_is_recorded_with_headroom_over_the_measured_rate(bar):
    feas = bar["feasibility_established_first"]
    assert feas["locate_floor"] == 0.95
    assert feas["locate"]["unique_rate"] > feas["locate_floor"]


# ------------------------------------------------------------------------------- forecast


def test_the_forecast_covers_every_outcome_and_sums_to_one(bar):
    probs = {k[2:]: v for k, v in bar["forecast"].items() if k.startswith("p_")}
    assert set(probs) == set(bar["outcomes"]), (
        "forecast/outcome mismatch: %r" % (set(probs) ^ set(bar["outcomes"]))
    )
    assert sum(probs.values()) == pytest.approx(1.0, abs=1e-6)
    assert all(0.0 <= v <= 1.0 for v in probs.values())


def test_the_forecast_states_the_standing_record_and_its_direction(bar):
    note = bar["forecast"]["note"]
    assert "2 of 7" in note
    assert "PESSIMISTIC" in note
    assert "not to be deferred to" in note or "never" in note


# ---------------------------------------------------------------- pre-commitments and limits


def test_the_bounded_verdict_pre_commitment_is_present_and_explicit(bar):
    joined = " ".join(bar["pre_commitments"])
    assert "may NEVER be written as 'position excluded'" in joined
    assert "ON THE SAME STRATUM" in joined
    assert "arms.native150 keeps its measured value" in joined
    assert "Neither control may be dropped" in joined


def test_the_pre_commitment_against_weakening_the_test_is_present(bar):
    joined = " ".join(bar["pre_commitments"])
    assert "never a weaker test" in joined


def test_the_honest_limits_name_the_scope_of_the_replica(bar):
    joined = " ".join(bar["honest_limits_fixed_in_advance"])
    assert "AS SPECIFIED" in joined
    assert "as the authors ran it" in joined
    assert "689" in joined, "the duplicate-sequence property of the substrate"
    assert "325" in joined and "31" in joined, "the excluded counts"


def test_the_deviation_from_the_plan_is_recorded_rather_than_taken_quietly(bar):
    dev = bar["deviation_from_the_technical_plan"]
    assert "SEVEN" in dev["what"] and "EIGHT" in dev["what"]
    assert "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE" in dev["what"]
    assert "ratify" in dev["status"], "the deviation must be surfaced for a decision"
    assert dev["what"].count("INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE") >= 1
