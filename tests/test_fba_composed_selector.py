"""Tests for the composed expression selector (`scripts/fba_composed_selector.py`).

The whole point of this script is `compose_safe_set` — the joint-verification step both prior gating
attempts lacked. The trap it exists to catch is real and was measured: `gpr_disabled_reactions` is a
SINGLE-GENE property, so gating one member of an isozyme pair is harmless while gating BOTH kills the
reaction. Per-element safety does not compose over a set.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

pytest.importorskip("cobra", reason="cobrapy required")

from scripts.fba_composed_selector import compose_safe_set, disabled_by_set  # noqa: E402
from scripts.fba_orphan_protection_screen import gpr_disabled_reactions  # noqa: E402


@pytest.fixture()
def isozyme_model():
    """R1 is carried by an isozyme pair (gA or gB); R2 is carried by gC alone."""
    import cobra

    m = cobra.Model("iso")
    m.compartments = {"c": "cytosol"}
    A, B, C = (cobra.Metabolite(x, compartment="c") for x in ("A", "B", "C"))

    def add(rid, st, rule):
        r = cobra.Reaction(rid)
        r.lower_bound, r.upper_bound = 0.0, 1000.0
        m.add_reactions([r])
        r.add_metabolites(st)
        r.gene_reaction_rule = rule
        return r

    add("R1", {A: -1, B: 1}, "gA or gB")
    add("R2", {B: -1, C: 1}, "gC")
    return m


def test_each_isozyme_alone_disables_nothing(isozyme_model):
    """The premise: both pass the single-gene eligibility check."""
    assert gpr_disabled_reactions(isozyme_model, "gA") == []
    assert gpr_disabled_reactions(isozyme_model, "gB") == []


def test_but_gating_BOTH_disables_the_reaction(isozyme_model):
    """The measured failure: single-gene safety does not compose."""
    assert disabled_by_set(isozyme_model, {"gA", "gB"}) == ["R1"]


def test_compose_drops_one_member_of_the_collision(isozyme_model):
    """The fix. One survives, the reaction lives, and exactly one gene was dropped."""
    expr = {"gA": 5.0, "gB": 1.0}
    safe, dropped = compose_safe_set(isozyme_model, {"gA", "gB"}, expr)
    assert disabled_by_set(isozyme_model, safe) == []
    assert len(safe) == 1
    assert dropped == 1


def test_compose_drops_the_HIGHEST_expressed_member(isozyme_model):
    """Pre-registered resolution rule (§2 step 3): the higher-expressed gene is the one more likely to
    be present, so it is the one restored."""
    safe, _ = compose_safe_set(isozyme_model, {"gA", "gB"}, {"gA": 5.0, "gB": 1.0})
    assert safe == {"gB"}
    safe, _ = compose_safe_set(isozyme_model, {"gA", "gB"}, {"gA": 1.0, "gB": 5.0})
    assert safe == {"gA"}


def test_compose_is_deterministic_on_an_expression_tie(isozyme_model):
    """A tie must not make the run non-reproducible — the determinism gate depends on this."""
    tie = {"gA": 1.0, "gB": 1.0}
    first = compose_safe_set(isozyme_model, {"gA", "gB"}, tie)[0]
    for _ in range(5):
        assert compose_safe_set(isozyme_model, {"gA", "gB"}, tie)[0] == first


def test_a_non_colliding_set_is_returned_untouched(isozyme_model):
    """The composition step must not shrink a set that was already safe."""
    safe, dropped = compose_safe_set(isozyme_model, {"gA"}, {"gA": 1.0})
    assert safe == {"gA"} and dropped == 0


def test_empty_candidate_set_is_safe(isozyme_model):
    assert disabled_by_set(isozyme_model, set()) == []
    assert compose_safe_set(isozyme_model, set(), {}) == (set(), 0)


def test_a_sole_route_gene_is_never_eligible_in_the_first_place(isozyme_model):
    """gC alone carries R2, so the single-gene filter excludes it before composition runs. This is why
    the selector cannot gate away capability."""
    assert gpr_disabled_reactions(isozyme_model, "gC") == ["R2"]
