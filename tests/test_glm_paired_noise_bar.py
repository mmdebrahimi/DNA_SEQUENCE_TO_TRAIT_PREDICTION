"""The paired-noise pre-registration must stay a PRE-registration, and stay SCALE-FREE.

Two defects these tests exist to prevent, both recorded:

  1. A BAND DERIVED FROM THE WRONG QUANTITY. The v1 draft of this plan judged a two-model
     difference against the seed spread of ONE model. That is the documented error in
     `feedback_tolerance_band_from_the_tested_quantity`, whose first recorded instance was this
     same generator. The fix is a per-seed PAIRED delta, and these tests pin that the file says so.

  2. A NUMERIC THRESHOLD FROZEN BEFORE THE SCALE WAS MEASURED. The magnitude of the sampling
     variation is what the run measures, so a constant band frozen now would be a guess. The
     readable outcomes must therefore compare the delta to its OWN measured spread.

A NOTE ON HOW NOT TO TEST THIS. A substring search cannot distinguish an assertion from a denial:
the file contains the word "clearable" only inside "is NOT restated here as clearable", so grepping
for it returns a misleading True. Nor can "the conditions contain no decimals" work, because
`+0.0271` appears legitimately to fix the SIGN of the expected effect and the reproduction checks
legitimately cite 0.6679 / 0.5865. So these tests assert on STRUCTURE -- which key compares against
what -- rather than on the presence or absence of tokens.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BAR = ROOT / "wiki" / "glm_generator_paired_noise_bar_2026-11-20.json"

#: The two committed arms this pre-registration is about. Mode-both, lower is better.
COMMITTED_ZEROSHOT = 0.7787
COMMITTED_FINETUNED = 0.8058


@pytest.fixture(scope="module")
def bar() -> dict:
    assert BAR.is_file(), "the pre-registration is missing: %s" % BAR
    return json.loads(BAR.read_text(encoding="utf-8"))


def test_the_bar_parses_and_declares_itself_frozen_before_any_number(bar):
    assert bar["record"] == "glm-generator-paired-noise-bar-v1"
    assert bar["frozen"] == "2026-11-20"
    fb = bar["frozen_before"].lower()
    # It must say the run has NOT happened yet, not merely carry a date.
    assert "no multi-seed run has been executed" in fb
    assert "not the post-hoc reading" in fb or "not a post-hoc reading" in fb


def test_the_decision_quantity_is_a_PER_SEED_DIFFERENCE_not_one_arms_spread(bar):
    """Defect 1. The whole point of v2 over v1."""
    dq = bar["decision_quantity"]
    assert "delta_s" in dq["definition"]
    assert "PER SEED" in dq["definition"] or "per seed" in dq["definition"]
    # The subtraction must be present in the definition, not merely implied in prose.
    assert "-" in dq["definition"] and "score_zs" in dq["definition"]
    # And the file must state WHY pairing is the fix, naming the variance identity.
    assert "Var(A - B)" in dq["why_paired"]
    assert "NOT the uncertainty" in dq["why_paired"]


def test_the_bar_names_the_anti_pattern_it_was_built_to_avoid(bar):
    ap = bar["decision_quantity"]["anti_pattern_this_avoids"]
    assert "feedback_tolerance_band_from_the_tested_quantity" in ap
    assert "v1" in ap, "the file should record that its own first draft made this error"


@pytest.mark.parametrize("outcome", ["DELTA_ROBUST", "DELTA_NOT_ESTABLISHED"])
def test_every_readable_delta_outcome_compares_against_the_MEASURED_spread(bar, outcome):
    """Defect 2, in its positive form.

    This is the load-bearing assertion. A readable verdict must be reached by comparing the delta
    to the spread the run itself measures -- never to a constant frozen before the scale was known.
    Asserting "no decimals appear" would be the wrong test (see the module docstring); asserting
    "the comparison is against the spread" is the property actually wanted.
    """
    cond = bar["delta_outcomes"][outcome]["condition"]
    assert "spread of delta_s" in cond, (
        "%s must be decided against the measured spread, not a constant: %r" % (outcome, cond))


def test_a_numeric_band_is_explicitly_refused(bar):
    key = "no_numeric_band_is_frozen_here"
    assert key in bar, "the refusal must be explicit, not merely implied by the conditions"
    txt = bar[key].lower()
    assert "deliberate" in txt
    assert "mis-specified" in txt or "guessing at the scale" in txt


def test_instrument_invalid_covers_a_missing_cell_and_is_not_averaged_over(bar):
    inv = bar["delta_outcomes"]["INSTRUMENT_INVALID"]
    assert "missing" in inv["condition"]
    assert "not averaged over" in inv["reading"]


def test_all_four_trajectory_outcomes_exist_and_are_paired(bar):
    traj = bar["trajectory_outcomes"]
    for name in ("OVERFIT_MONOTONE", "FLAT_NOT_EPOCHS",
                 "HELPS_THEN_DEGRADES", "TRAJECTORY_UNRESOLVED"):
        assert name in traj, "missing trajectory outcome: %s" % name
    assert "paired" in traj["scored_on"]
    assert "same seed set" in traj["scored_on"]


def test_the_unresolved_trajectory_verdict_prescribes_more_draws_not_a_conclusion(bar):
    reading = bar["trajectory_outcomes"]["TRAJECTORY_UNRESOLVED"]["reading"]
    assert "NOT a conclusion" in reading


def test_the_seed_set_contains_the_committed_arms_seed(bar):
    seeds = bar["seed_set"]
    assert len(seeds) == len(set(seeds)), "duplicate seeds: %r" % seeds
    assert len(seeds) >= 3, "fewer than 3 seeds cannot bound a spread"
    assert 7 in seeds, (
        "7 is the falsifier CLI default and therefore the seed both committed arms ran at; "
        "without it the anchors are not checkable")


def test_the_anchor_rule_names_both_committed_numbers(bar):
    anchor = bar["anchor_rule"]["ANCHOR_DRIFT"]
    assert str(COMMITTED_ZEROSHOT) in anchor
    assert str(COMMITTED_FINETUNED) in anchor
    # It must be checked against the measured spread, not for exact equality.
    assert "rather than for exact equality" in anchor


def test_the_variant_is_branded_and_does_not_claim_the_frozen_arm(bar):
    assert bar["claim_status"] == "DIAGNOSTIC_METHOD_VARIANT"
    nfa = bar["not_the_frozen_arm"]
    assert "method_frozen.epochs = 3" in nfa
    assert "Epochs 1 and 2" in nfa
    # Epoch 3 IS the frozen method and the file must say so, else the one legitimate
    # comparison to the committed arm gets disclaimed away along with the other two.
    assert "Epoch 3 IS the frozen method" in nfa


def test_the_bar_states_what_it_does_NOT_bound(bar):
    """Seed pairing varies GENERATION only. One training trajectory means training noise
    is unmeasured, and a DELTA_ROBUST verdict must not be read as fully bounded."""
    txt = bar["what_this_design_does_NOT_bound"]
    assert "TRAINING variance" in txt
    assert "ONE training trajectory" in txt
    assert "out of scope" in txt


def test_a_forecast_is_recorded_with_reasoning(bar):
    f = bar["forecast"]
    assert 0.0 < f["P_DELTA_ROBUST"] < 1.0
    # Reasoning must contain BOTH directions, or it is not a forecast, it is a prediction.
    assert "For:" in f["reasoning"] and "Against:" in f["reasoning"]
    assert "retro-fitted" in f["reasoning"]


def test_the_long_range_blindness_limit_survives_into_this_bar(bar):
    """Inherited from the arm artifacts: a k-mer discriminator sees LOCAL composition only, so a
    low score is necessary-not-sufficient. Dropping it here would let a good number over-claim."""
    joined = " ".join(bar["honest_limits"])
    assert "LONG-RANGE" in joined
    assert "necessary, NOT sufficient" in joined


def test_this_bar_does_not_supersede_the_frozen_acceptance_bar(bar):
    txt = bar["supersedes_nothing"]
    assert "glm_generator_finetune_acceptance_bar.json" in txt
    assert "0.6679" in txt, "it must name the bar it is NOT replacing"
