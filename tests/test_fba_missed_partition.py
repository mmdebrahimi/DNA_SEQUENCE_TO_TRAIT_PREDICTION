"""Tests for the missed-gene partition (`scripts/fba_missed_gene_partition.py`).

The whole value of this script is that its three classes imply DIFFERENT remedial work, so a
misclassification would point the next experiment at the wrong layer. `classify` is pure; these pin the
boundaries, including the one the real run turned out to be empty (MIS_CONDITIONED).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.fba_missed_gene_partition import classify  # noqa: E402


def test_no_prediction_anywhere_is_never_fires():
    assert classify(set(), {"glucose"}) == "NEVER_FIRES"


def test_firing_only_in_the_wrong_conditions_is_mis_conditioned():
    """The class the real run found EMPTY (0 of 131). Kept tested so the distinction stays real rather
    than becoming an unreachable branch nobody notices."""
    assert classify({"acetate"}, {"glucose"}) == "MIS_CONDITIONED"


def test_any_correct_condition_is_partial_overlap():
    assert classify({"glucose", "acetate"}, {"glucose"}) == "PARTIAL_OVERLAP"


def test_exact_match_is_partial_overlap_not_a_separate_class():
    """A perfect match is still 'overlaps'; the script deliberately does not carve out a PERFECT class,
    because the per-cell catch rate is reported separately."""
    assert classify({"glucose"}, {"glucose"}) == "PARTIAL_OVERLAP"


def test_never_fires_beats_the_empty_truth_edge_case():
    """A gene with no true conditions cannot overlap; an empty prediction is still NEVER_FIRES."""
    assert classify(set(), set()) == "NEVER_FIRES"
    assert classify({"glucose"}, set()) == "MIS_CONDITIONED"
