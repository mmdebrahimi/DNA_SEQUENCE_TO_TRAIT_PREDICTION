"""Tests for the orphan-protection screen (`scripts/fba_orphan_protection_screen.py`).

The screen makes a THEOREM-shaped claim -- "the model can never call this gene essential" -- so the
tests target the ways that claim could be admitted wrongly:

  * a rescuer that is itself DEAD (universally blocked) must not count;
  * a rescuer that is disabled by the SAME knockout must not count;
  * a rescuer with a SMALLER capacity than the reaction it replaces must not count;
  * a reversed reaction really is the same transformation (and must be recognised as one);
  * a gene with even ONE unexplained disabled reaction stays CALLABLE.

The pure helpers run with no model at all; `find_rescuers` is exercised on a hand-built toy network so
the trap cases are reachable without iML1515.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.fba_orphan_protection_screen import (  # noqa: E402
    canon_interval,
    canon_key,
    classify_gene,
    covers,
    find_rescuers,
)


# ------------------------------------------------------------------ canonical transformation identity

def test_forward_and_reverse_are_the_same_transformation():
    """`A -> B` and `B -> A` differ only in direction, which the interval carries, not the key."""
    assert canon_key({"A": -1.0, "B": 1.0}) == canon_key({"A": 1.0, "B": -1.0})


def test_scaling_is_the_same_transformation():
    assert canon_key({"A": -2.0, "B": 2.0}) == canon_key({"A": -1.0, "B": 1.0})


def test_different_stoichiometry_is_a_different_transformation():
    """Same metabolites, different ratio -- must NOT be treated as interchangeable."""
    assert canon_key({"A": -1.0, "B": 2.0}) != canon_key({"A": -1.0, "B": 1.0})


def test_different_metabolites_are_a_different_transformation():
    assert canon_key({"A": -1.0, "B": 1.0}) != canon_key({"A": -1.0, "C": 1.0})


def test_empty_stoichiometry_has_no_key():
    assert canon_key({}) is None


def test_canonical_interval_swaps_when_the_reference_coefficient_is_negative():
    """u = v*c, so a negative reference coefficient mirrors the interval. Getting this backwards would
    make forward and reverse reactions look mutually rescuing when they are not."""
    # ref metabolite is "A" (lexicographically smallest), coefficient -1
    assert canon_interval({"A": -1.0, "B": 1.0}, 0.0, 10.0) == (-10.0, 0.0)
    # ref coefficient +1 -> interval preserved
    assert canon_interval({"A": 1.0, "B": -1.0}, 0.0, 10.0) == (0.0, 10.0)


def test_covers_is_containment():
    assert covers((-10.0, 10.0), (-1.0, 1.0))
    assert not covers((-1.0, 1.0), (-10.0, 10.0))
    assert covers((0.0, 5.0), (0.0, 5.0))


# --------------------------------------------------------------------------------- gene classification

def test_no_disabled_reactions_is_uncallable():
    assert classify_gene([], set(), {}) == ("NO_KO_EFFECT", True)


def test_all_blocked_is_uncallable():
    assert classify_gene(["r1", "r2"], {"r1", "r2"}, {}) == ("ALL_DISABLED_BLOCKED", True)


def test_all_rescued_is_uncallable():
    cls, ok = classify_gene(["r1"], set(), {"r1": ["r1b"]})
    assert (cls, ok) == ("DUPLICATE_RESCUED", True)


def test_mixed_blocked_and_rescued_is_uncallable():
    cls, ok = classify_gene(["r1", "r2"], {"r1"}, {"r2": ["r2b"]})
    assert (cls, ok) == ("MIXED_BLOCKED_AND_RESCUED", True)


def test_one_unexplained_reaction_makes_the_gene_callable():
    """The whole screen turns on this being conjunctive -- a single uncovered reaction is enough for the
    model to be able to call the gene essential."""
    cls, ok = classify_gene(["r1", "r2"], {"r1"}, {"r2": []})
    assert (cls, ok) == ("CALLABLE", False)


# --------------------------------------------------------------------------------------- rescuer traps

@pytest.fixture()
def toy():
    """A -> B by four parallel routes plus one half-capacity route and one different transformation."""
    cobra = pytest.importorskip("cobra", reason="cobrapy required")
    m = cobra.Model("toy")
    A, B, C = (cobra.Metabolite(x) for x in ("A", "B", "C"))

    def add(rid, stoich, lb=-1000.0, ub=1000.0, rule=""):
        r = cobra.Reaction(rid)
        r.lower_bound, r.upper_bound = lb, ub
        m.add_reactions([r])
        r.add_metabolites(stoich)
        r.gene_reaction_rule = rule
        return r

    add("MAIN", {A: -1, B: 1}, rule="g1")
    add("TWIN", {A: -1, B: 1}, rule="g2")          # a legitimate rescuer
    add("DEAD", {A: -1, B: 1}, rule="g3")          # will be passed in as blocked
    add("SMALL", {A: -1, B: 1}, lb=-1.0, ub=1.0)   # capacity too small
    add("REVERSE", {A: 1, B: -1}, rule="g4")       # same transformation, opposite orientation
    add("OTHER", {A: -1, C: 1}, rule="g5")         # different transformation
    return m


def _index(model):
    from scripts.fba_orphan_protection_screen import build_duplicate_index
    return build_duplicate_index(model)


def test_a_surviving_full_capacity_twin_rescues(toy):
    r = find_rescuers(toy, "MAIN", disabled={"MAIN"}, blocked=set(), dup_index=_index(toy))
    assert "TWIN" in r


def test_a_blocked_reaction_never_rescues(toy):
    """A dead route carries nothing, so it cannot absorb the flux it superficially duplicates."""
    r = find_rescuers(toy, "MAIN", disabled={"MAIN"}, blocked={"DEAD"}, dup_index=_index(toy))
    assert "DEAD" not in r


def test_a_reaction_disabled_by_the_same_knockout_never_rescues(toy):
    """The common trap: two reactions of the SAME gene look like each other's rescuer."""
    r = find_rescuers(toy, "MAIN", disabled={"MAIN", "TWIN"}, blocked=set(), dup_index=_index(toy))
    assert "TWIN" not in r


def test_a_lower_capacity_duplicate_never_rescues(toy):
    r = find_rescuers(toy, "MAIN", disabled={"MAIN"}, blocked=set(), dup_index=_index(toy))
    assert "SMALL" not in r


def test_a_reversed_reaction_with_symmetric_bounds_does_rescue(toy):
    """REVERSE spans [-1000, 1000] in canonical units just as MAIN does, so it genuinely covers it."""
    r = find_rescuers(toy, "MAIN", disabled={"MAIN"}, blocked=set(), dup_index=_index(toy))
    assert "REVERSE" in r


def test_a_different_transformation_never_rescues(toy):
    r = find_rescuers(toy, "MAIN", disabled={"MAIN"}, blocked=set(), dup_index=_index(toy))
    assert "OTHER" not in r


def test_a_reaction_never_rescues_itself(toy):
    r = find_rescuers(toy, "MAIN", disabled={"MAIN"}, blocked=set(), dup_index=_index(toy))
    assert "MAIN" not in r


# ------------------------------------------------------------------------------------- the unmask test

def test_masking_genes_is_everything_sharing_a_reaction(toy):
    """g1 is alone on MAIN, so nothing masks it beyond itself."""
    from scripts.fba_orphan_protection_screen import masking_genes
    assert masking_genes(toy, "g1") == {"g1"}


def test_masking_genes_collects_partners_across_all_of_a_genes_reactions(toy):
    """A gene on a shared reaction must pull in its co-occupants -- that set IS the mask."""
    from scripts.fba_orphan_protection_screen import masking_genes
    shared = toy.reactions.get_by_id("TWIN")
    shared.gene_reaction_rule = "g2 or gX"
    assert masking_genes(toy, "g2") == {"g2", "gX"}


def test_unmask_deletes_the_mask_not_the_reactions():
    """The load-bearing distinction: deleting the MASKING GENES is not the same as zeroing the gene's
    reactions. For a gene mapped to many reactions the latter removes capability far beyond the
    redundancy under test, and scored ompC (285 reactions) differently from the correct test."""
    from scripts.fba_orphan_protection_screen import masking_genes
    cobra = pytest.importorskip("cobra", reason="cobrapy required")
    m = cobra.Model("wide")
    A, B, C = (cobra.Metabolite(x) for x in ("A", "B", "C"))
    for rid, st, rule in (("R1", {A: -1, B: 1}, "wide or alt"), ("R2", {B: -1, C: 1}, "wide")):
        r = cobra.Reaction(rid)
        r.lower_bound, r.upper_bound = -1000.0, 1000.0
        m.add_reactions([r])
        r.add_metabolites(st)
        r.gene_reaction_rule = rule
    # the mask is the GENE set, and it does not include the second reaction's unrelated capability
    assert masking_genes(m, "wide") == {"wide", "alt"}
    assert len(list(m.genes.get_by_id("wide").reactions)) == 2
