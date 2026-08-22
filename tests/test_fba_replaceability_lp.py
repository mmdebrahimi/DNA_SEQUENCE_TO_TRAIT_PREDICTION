"""Tests for the shadow-reaction replaceability LP (`scripts/fba_replaceability_lp.py`).

The LP underwrites a claim quantified over ALL media ("no constraint-based method can ever call this gene
essential"), so the tests target the ways it could say YES wrongly:

  * it must find a genuine multi-step bypass (otherwise it adds nothing over the exact-duplicate screen);
  * it must return 0 when no bypass exists;
  * a NARROWER bypass must not read as full replaceability;
  * the reactions removed alongside the target must not be able to rescue it;
  * closing the exchanges must actually bite -- a route that needs the medium is not medium-independent.

The last one is tested as a CONTRAST between two models differing in exactly one thing: whether the
bypass's cofactor comes from outside the cell or inside it. Anything less would pass vacuously.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

pytest.importorskip("cobra", reason="cobrapy required")

from scripts.fba_replaceability_lp import capacity, replaceable_rate  # noqa: E402


def _model(bypass_bound: float = 1000.0, with_bypass: bool = True, cofactor: str | None = None):
    """`TARGET` converts A -> B. Optional bypasses do the same net conversion by another route.

    cofactor="medium"   bypass is `A + X_e -> B`, and X_e can ONLY come from the medium.
    cofactor="internal" bypass is `A + X_c -> B`, and X_c is made inside the cell.
    Same reaction shape either way; only the cofactor's origin differs.
    """
    import cobra

    m = cobra.Model("t")
    # An unambiguous external compartment is required: with no metabolite in "e", cobra falls back to
    # "the compartment with the most boundary reactions" and can misclassify an internal demand reaction
    # as an exchange -- silently changing what the LP closes.
    m.compartments = {"c": "cytosol", "e": "extracellular"}
    A, B, C, X_c = (cobra.Metabolite(x, compartment="c") for x in ("A", "B", "C", "X_c"))
    D_e, X_e = (cobra.Metabolite(x, compartment="e") for x in ("D_e", "X_e"))

    def add(rid, st, lb=0.0, ub=1000.0):
        r = cobra.Reaction(rid)
        r.lower_bound, r.upper_bound = lb, ub
        m.add_reactions([r])
        r.add_metabolites(st)
        return r

    add("TARGET", {A: -1, B: 1})
    add("EX_D_e", {D_e: -1}, lb=-1000.0, ub=1000.0)      # pins "e" as the external compartment
    if with_bypass:
        add("BY1", {A: -1, C: 1}, ub=bypass_bound)
        add("BY2", {C: -1, B: 1}, ub=bypass_bound)
    if cofactor == "medium":
        add("EX_X_e", {X_e: -1}, lb=-1000.0, ub=1000.0)  # the ONLY source of X_e
        add("BY_MED", {A: -1, X_e: -1, B: 1})
    elif cofactor == "internal":
        add("XSRC", {X_c: 1})                            # X made inside; a demand, not an exchange
        add("BY_INT", {A: -1, X_c: -1, B: 1})
    return m


def test_a_multi_step_bypass_is_found():
    """The whole point: functional redundancy that is NOT an exact stoichiometric duplicate."""
    m = _model()
    assert replaceable_rate(m, "TARGET", set()) >= capacity(m, "TARGET") - 1e-7


def test_no_bypass_gives_zero():
    m = _model(with_bypass=False)
    assert replaceable_rate(m, "TARGET", set()) == pytest.approx(0.0, abs=1e-7)


def test_a_narrower_bypass_is_not_full_replaceability():
    """A route that exists but cannot carry the target's capacity must not count as replaceable."""
    m = _model(bypass_bound=1.0)
    rate = replaceable_rate(m, "TARGET", set())
    assert 0.0 < rate < capacity(m, "TARGET")


def test_co_disabled_reactions_cannot_rescue_the_target():
    """A gene's own reactions must not rescue each other -- that would make every multi-reaction gene
    look replaceable by itself."""
    m = _model()
    assert replaceable_rate(m, "TARGET", {"BY1"}) == pytest.approx(0.0, abs=1e-7)


def test_a_medium_dependent_route_is_refused_but_the_same_route_internally_is_not():
    """Closing the exchanges is what makes the claim hold for ANY medium.

    Both models carry a bypass of the SAME shape (`A + X -> B`). It counts when X is made inside the cell
    and is refused when X can only come from the medium. Asserting BOTH directions is what stops this
    passing on a model that simply had no bypass at all.
    """
    med = _model(with_bypass=False, cofactor="medium")
    assert replaceable_rate(med, "TARGET", set()) == pytest.approx(0.0, abs=1e-7)

    internal = _model(with_bypass=False, cofactor="internal")
    assert replaceable_rate(internal, "TARGET", set()) >= capacity(internal, "TARGET") - 1e-7


def test_capacity_reads_the_wider_bound():
    m = _model()
    m.reactions.get_by_id("TARGET").lower_bound = -5.0
    assert capacity(m, "TARGET") == 1000.0
    m.reactions.get_by_id("TARGET").upper_bound = 2.0
    assert capacity(m, "TARGET") == 5.0
