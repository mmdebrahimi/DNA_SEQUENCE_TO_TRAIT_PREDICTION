"""The human-falsifier bar must stay a PRE-registration, and must keep forbidding the easy comparison.

THE DEFECT THESE TESTS EXIST TO PREVENT is a cross-substrate comparison that looks right and is
wrong. Putting the human distinguishability level next to the bacterial 0.6679 is one subtraction
away at all times. The predecessor memo (wiki/eukaryotic_arm_result_K562_2026-11-01.md) states that
only a WITHIN-substrate contrast crosses a substrate boundary, and these two substrates also differ
in length (230 bp vs 150 bp), which changes a k-mer discriminator's job directly.

A SECOND defect: freezing a numeric tolerance before the scale is known. The tolerance here is the
human self-control's own seed spread, derived in the run. That is the same error recorded in
feedback_tolerance_band_from_the_tested_quantity, whose first instance was this generator family.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BAR = ROOT / "wiki" / "human_falsifier_bar_2026-11-21.json"
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"

#: The committed bacterial values the bar's anchors are derived from, mode-both.
BACT_NULL_BOTH = 0.6679
BACT_SELF_BOTH = 0.5865
BACT_NULL_KMER = 0.5990
BACT_SELF_KMER = 0.5447


@pytest.fixture(scope="module")
def bar() -> dict:
    assert BAR.is_file(), "the pre-registration is missing: %s" % BAR
    return json.loads(BAR.read_text(encoding="utf-8"))


def test_the_bar_parses_and_says_nothing_was_measured_yet(bar):
    assert bar["record"] == "human-falsifier-margin-bar-v1"
    assert bar["frozen"] == "2026-11-21"
    fb = bar["frozen_before"]
    # It must record the VERIFIED fact that no human measurement exists, not merely a date.
    assert "ZERO references to human data" in fb
    assert "not a post-hoc reading" in fb or "not the post-hoc reading" in fb


def test_the_decision_quantity_is_a_WITHIN_SUBSTRATE_MARGIN(bar):
    dq = bar["decision_quantity"]
    assert "null_distinguishability - self_control_ceiling" in dq["definition"]
    assert "SEPARATELY per substrate" in dq["definition"]
    # It must cite the predecessor's rule rather than asserting the methodology itself.
    why = dq["why_a_substrate_boundary"] if "why_a_substrate_boundary" in dq else dq["why_a_margin_and_not_a_level"]
    assert "eukaryotic_arm_result_K562_2026-11-01" in why
    assert "within-substrate" in why
    assert "230 bp vs 150 bp" in why, "the length difference must be named as a reason"


def test_the_bacterial_anchors_are_arithmetically_correct(bar):
    """A margin quoted wrong would silently shift every verdict."""
    a = bar["decision_quantity"]["bacterial_anchors"]
    both, kmer = a["mode_both"], a["mode_kmer"]
    assert both["null"] == BACT_NULL_BOTH and both["self_control"] == BACT_SELF_BOTH
    assert kmer["null"] == BACT_NULL_KMER and kmer["self_control"] == BACT_SELF_KMER
    assert abs(both["margin"] - (BACT_NULL_BOTH - BACT_SELF_BOTH)) < 1e-9
    assert abs(kmer["margin"] - (BACT_NULL_KMER - BACT_SELF_KMER)) < 1e-9


def test_the_bacterial_anchors_match_the_committed_grid(bar):
    """The anchors must come from the committed artifact, not from prose."""
    if not GRID.is_file():
        pytest.skip("the bacterial grid artifact is absent")
    g = json.loads(GRID.read_text(encoding="utf-8"))["grid"]
    nulls = {round(c["null"]["both"], 4) for c in g if c.get("null") and c["null"].get("both")}
    assert BACT_NULL_BOTH in nulls, (
        "the bar's bacterial null %.4f is not among the grid's observed nulls %s"
        % (BACT_NULL_BOTH, sorted(nulls)))
    a = bar["decision_quantity"]["bacterial_anchors"]
    assert "glm_generator_paired_sweep_2026-11-20.json" in a["source"]


def test_the_invalid_cross_substrate_comparison_is_EXPLICITLY_PROHIBITED(bar):
    """The load-bearing assertion. The easy comparison must be forbidden in the bar itself."""
    p = bar["PROHIBITED_COMPARISON"]
    assert "0.6679" in p, "the prohibition must name the number it forbids comparing to"
    assert "may NOT be compared" in p
    # It must say WHY, or a later reader will treat it as arbitrary and drop it.
    assert "length" in p and "assay" in p


def test_no_numeric_tolerance_is_frozen(bar):
    key = "no_numeric_tolerance_is_frozen_here"
    assert key in bar, "the refusal to freeze a tolerance must be explicit"
    txt = bar[key]
    assert "seed spread" in txt
    assert "derived in the run" in txt


@pytest.mark.parametrize("outcome", ["INSTRUMENT_INVALID", "CHAIN_IS_WEAKER_ON_HUMAN",
                                     "CHAIN_IS_COMPARABLE", "CHAIN_IS_STRONGER_ON_HUMAN"])
def test_all_four_outcomes_exist_with_a_condition_and_a_reading(bar, outcome):
    o = bar["outcomes_ordered_first_match_wins"][outcome]
    assert o["condition"].strip()
    assert o["reading"].strip()


def test_every_readable_outcome_is_decided_against_the_MEASURED_spread(bar):
    """Scale-free by construction: the comparison is margin-vs-margin against a spread the run
    measures, never against a constant invented here."""
    outs = bar["outcomes_ordered_first_match_wins"]
    for name in ("CHAIN_IS_WEAKER_ON_HUMAN", "CHAIN_IS_COMPARABLE", "CHAIN_IS_STRONGER_ON_HUMAN"):
        cond = outs[name]["condition"]
        assert "spread" in cond, "%s must be decided against the measured spread: %r" % (name, cond)


def test_the_favourable_outcome_disclaims_what_it_does_not_show(bar):
    """CHAIN_IS_WEAKER_ON_HUMAN is the result most likely to be over-read as a green light."""
    o = bar["outcomes_ordered_first_match_wins"]["CHAIN_IS_WEAKER_ON_HUMAN"]
    assert "explicitly_does_not_show" in o
    assert "would SUCCEED" in o["explicitly_does_not_show"]


def test_the_uniform_control_is_the_instrument_check(bar):
    inv = bar["outcomes_ordered_first_match_wins"]["INSTRUMENT_INVALID"]
    assert "uniform control" in inv["condition"]
    assert "0.9641" in inv["reading"], "the bacterial uniform value must be named as the expectation"


def test_the_geometry_decisions_are_made_IN_ADVANCE(bar):
    """Adaptor handling and cell precedence are both choosable after seeing results, so the bar
    must fix them first."""
    g = bar["geometry_decisions_the_plan_delegated_to_this_bar"]
    primary = g["adaptors_primary_is_FULL_230bp"]
    assert "PRIMARY is the full 230 bp" in primary
    assert "SECONDARY" in primary, "the 200 nt element must still be measured"
    assert "cannot be chosen after seeing" in primary
    cells = g["both_cells_reported_and_disagreement_resolves_CONSERVATIVELY"]
    assert "WEAKER + STRONGER reports CHAIN_IS_COMPARABLE" in cells
    assert "fixed in advance" in cells


def test_the_bar_refuses_to_borrow_the_bacterial_ceiling(bar):
    """There is no human replicate ceiling; the bacterial 0.941 must not be substituted."""
    joined = " ".join(bar["what_this_design_does_NOT_bound"])
    assert "0.941" in joined
    assert "must NOT be borrowed" in joined
    assert "human_frame_shape_2026-11-16" in joined


def test_the_no_generator_scope_choice_is_recorded(bar):
    joined = " ".join(bar["what_this_design_does_NOT_bound"])
    assert "NO generator is run" in joined
    assert "GENERator-v2-prokaryote-1.2b-base" in joined, (
        "the reason must be concrete -- the current generator is prokaryote-locked")


def test_the_long_range_blindness_limit_survives(bar):
    joined = " ".join(bar["what_this_design_does_NOT_bound"])
    assert "long-range" in joined.lower()
    assert "necessary-NOT-sufficient" in joined or "necessary, NOT sufficient" in joined


def test_a_two_sided_forecast_is_recorded(bar):
    f = bar["forecast"]
    assert 0.0 < f["P_CHAIN_IS_WEAKER_ON_HUMAN"] < 1.0
    assert "FOR:" in f["reasoning"] and "AGAINST:" in f["reasoning"]
    assert "retro-fitted" in f["reasoning"]
    # The for-case must cite the measured encoder numbers, not intuition alone.
    assert "0.4853" in f["reasoning"] and "0.2362" in f["reasoning"]


def test_the_published_incumbents_are_recorded_as_unavailable(bar):
    joined = " ".join(bar["honest_limits"])
    assert "UNRESOLVED_TABLE_IDENTITY_UNVERIFIED" in joined
