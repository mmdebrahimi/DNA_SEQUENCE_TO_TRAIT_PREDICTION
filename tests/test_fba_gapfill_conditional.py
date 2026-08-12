"""Gap-fill-vs-conditional-metric experiment -- pure logic (no cobra, no network)."""
from __future__ import annotations

from scripts.fba_gapfill_conditional_test import count_flips


def test_count_flips_is_zero_when_nothing_changed():
    """THE headline number of the experiment: every arm flipped zero calls, so this must be exact."""
    base = {"c1": {"a": True, "b": False}, "c2": {"a": False, "b": False}}
    assert count_flips(base, base) == 0


def test_count_flips_counts_each_changed_cell_in_both_directions():
    base = {"c1": {"a": True, "b": False}}
    arm = {"c1": {"a": False, "b": True}}          # essential->dispensable AND dispensable->essential
    assert count_flips(base, arm) == 2


def test_count_flips_ignores_genes_absent_from_the_arm():
    """An arm that could not score a gene must not be counted as having flipped it -- that would
    manufacture a change out of a missing measurement."""
    base = {"c1": {"a": True, "b": False}}
    assert count_flips(base, {"c1": {"a": True}}) == 0


def test_count_flips_ignores_conditions_absent_from_the_arm():
    base = {"c1": {"a": True}, "c2": {"a": True}}
    assert count_flips(base, {"c1": {"a": True}}) == 0


def test_a_ratio_change_that_does_not_cross_the_threshold_is_not_a_flip():
    """The mechanism the experiment uncovered: 78 of 268 ratios moved and zero calls flipped. Encoded
    here so the distinction between 'the model changed' and 'the metric changed' stays explicit."""
    frac = 0.01
    base_ratio, arm_ratio = 0.8860, 1.0000        # a real observed pair (b2277)
    base = {"c1": {"g": base_ratio <= frac}}
    arm = {"c1": {"g": arm_ratio <= frac}}
    assert base_ratio != arm_ratio                 # the model DID change
    assert count_flips(base, arm) == 0             # the metric did NOT
