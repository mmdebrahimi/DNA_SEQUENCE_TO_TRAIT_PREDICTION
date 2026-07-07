"""Pins the compound-allele resolver logic (`compound_caller`) independent of any gene.

The load-bearing rule: a haplotype's component-SNP SET resolves most-specific-first (the largest matching
compound wins), so {A,B} -> the compound, {A} -> the single, {} -> reference.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx.compound_caller import CompoundAllele, _haplotype_star  # noqa: E402

_RULES = [
    CompoundAllele("*3A", frozenset({"A", "B"})),
    CompoundAllele("*3B", frozenset({"A"})),
    CompoundAllele("*3C", frozenset({"B"})),
]


def test_both_components_pick_compound():
    assert _haplotype_star(frozenset({"A", "B"}), _RULES, "*1") == "*3A"


def test_single_component_picks_single():
    assert _haplotype_star(frozenset({"A"}), _RULES, "*1") == "*3B"
    assert _haplotype_star(frozenset({"B"}), _RULES, "*1") == "*3C"


def test_no_component_is_reference():
    assert _haplotype_star(frozenset(), _RULES, "*1") == "*1"


def test_most_specific_wins_over_subset():
    # {A,B} must resolve to *3A even though {A}->*3B and {B}->*3C both also match subsets
    assert _haplotype_star(frozenset({"A", "B"}), _RULES, "*1") == "*3A"


def test_superset_still_matches_compound():
    # an extra unmodelled component tag present alongside {A,B} still yields the compound
    assert _haplotype_star(frozenset({"A", "B", "X"}), _RULES, "*1") == "*3A"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
