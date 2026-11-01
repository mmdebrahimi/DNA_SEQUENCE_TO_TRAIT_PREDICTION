"""Reachability + ordering guards for the exposure grader.

**The load-bearing test is `test_every_outcome_is_reachable`.** Run 1 of the parent PromoterAtlas
arm emitted its strongest verdict while a validity gate was failing, because that gate's predicate
could never fire: the clause was present, ordered first, and dead. Presence is not enough, so every
outcome in the frozen check order must be reached by a constructed input, and the set of outcomes
the grader can actually produce must EQUAL the frozen order -- not merely be contained in it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.glm.exposure_verdict import (
    CHECK_ORDER,
    CHEAP_SEPARABILITY_FLOOR,
    ExposureInput,
    grade,
)

BAR = Path("wiki/promoteratlas_exposure_acceptance_bar.json")

#: One constructed input per outcome. Kept as data so the reachability test cannot silently cover
#: fewer outcomes than the order declares.
CASES: dict[str, ExposureInput] = {
    "INDETERMINATE_HARNESS_DRIFT": ExposureInput(
        harness_reconciled=False, sep_auroc_gc_only=0.90, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.5000, atlas_minus_random_dir=0.5070, gc_baseline=0.3112,
    ),
    "INDETERMINATE_CONTROLS_UNSEPARABLE": ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.52, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.5000, atlas_minus_random_dir=0.5070, gc_baseline=0.3112,
    ),
    "INDETERMINATE_REMOVAL_COST_UNBOUNDED": ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.90, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.4000, atlas_minus_random_dir=0.4000, gc_baseline=0.3112,
    ),
    "EXPOSURE_EXPLAINS": ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.90, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.3200, atlas_minus_random_dir=0.5070, gc_baseline=0.3112,
    ),
    "EXPOSURE_PARTIAL": ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.90, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.4300, atlas_minus_random_dir=0.5070, gc_baseline=0.3112,
    ),
    "EXPOSURE_EXCLUDED_LINEAR": ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.90, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.5000, atlas_minus_random_dir=0.5070, gc_baseline=0.3112,
    ),
}


@pytest.mark.parametrize("want", list(CASES))
def test_every_outcome_is_reachable(want):
    """No branch may be dead code. This is the defect that let run 1 of the parent arm pass."""
    assert grade(CASES[want]).verdict == want


def test_the_reachable_set_EQUALS_the_frozen_check_order():
    """Containment is not enough: an outcome declared in the bar but unreachable by the grader is
    exactly the run-1 defect, and an outcome the grader can emit that the bar never declared is an
    ungraded verdict."""
    assert set(CASES) == set(CHECK_ORDER)
    assert {grade(v).verdict for v in CASES.values()} == set(CHECK_ORDER)


def test_validity_outcomes_are_marked_invalid_and_substantive_ones_valid():
    for name, x in CASES.items():
        assert grade(x).is_valid is (not name.startswith("INDETERMINATE_"))


def test_harness_drift_wins_over_every_substantive_outcome():
    """A failing harness must not be overridden by a strong-looking number -- the precise shape of
    the run-1 failure, where the strongest verdict was emitted over a failing gate."""
    x = ExposureInput(
        harness_reconciled=False, sep_auroc_gc_only=0.99, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.5076, atlas_minus_random_dir=0.5076, gc_baseline=0.3112,
    )
    assert grade(x).verdict == "INDETERMINATE_HARNESS_DRIFT"


def test_unseparable_controls_win_over_exposure_explains():
    """A total collapse with unseparable controls says nothing about the encoder."""
    x = ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.50, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.3112, atlas_minus_random_dir=0.5070, gc_baseline=0.3112,
    )
    assert grade(x).verdict == "INDETERMINATE_CONTROLS_UNSEPARABLE"


def test_random_direction_cost_gates_the_substantive_branches():
    """If removing ANY direction costs more than the margin, a drop is not attributable."""
    x = ExposureInput(
        harness_reconciled=True, sep_auroc_gc_only=0.90, atlas_baseline=0.5076,
        atlas_minus_separation_dir=0.3100, atlas_minus_random_dir=0.4000, gc_baseline=0.3112,
    )
    assert grade(x).verdict == "INDETERMINATE_REMOVAL_COST_UNBOUNDED"


def test_excluded_linear_detail_always_carries_its_own_bound():
    """The pre-commitment is that this outcome may never be written as 'exposure excluded'. The
    grader's own detail string must carry the bound, so a downstream quoter cannot drop it."""
    d = grade(CASES["EXPOSURE_EXCLUDED_LINEAR"]).detail
    assert "nonlinear" in d and "never" in d.lower()


def test_cheap_separability_floor_is_above_chance():
    assert 0.5 < CHEAP_SEPARABILITY_FLOOR < 1.0


@pytest.mark.skipif(not BAR.is_file(), reason="bar not frozen yet")
def test_grader_order_matches_the_frozen_bar_exactly():
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert tuple(bar["verdict_check_order"]) == CHECK_ORDER
    assert set(bar["outcomes"]) == set(CHECK_ORDER)


@pytest.mark.skipif(not BAR.is_file(), reason="bar not frozen yet")
def test_frozen_bar_margin_is_inherited_not_invented():
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert bar["inherited_margin"] == 0.02
    assert "arm_acceptance_bar" in bar["inherited_margin_source"]


@pytest.mark.skipif(not BAR.is_file(), reason="bar not frozen yet")
def test_frozen_bar_preserves_the_parent_arms_measured_metric():
    """A measured metric is not retro-edited to match a later finding."""
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    names = {p["name"] for p in bar["pre_commitments"]}
    assert "a_measured_metric_is_not_retro_edited" in names
    assert "exclusion_is_linear_only" in names
    assert "non_vacuity_of_the_control" in names
