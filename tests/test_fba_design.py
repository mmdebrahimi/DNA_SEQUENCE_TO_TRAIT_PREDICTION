"""Growth-coupled strain design — pure-logic tests (wheel-only) + one real-model gate."""
from __future__ import annotations

import math

import pytest

from dna_decode.fba.design import (
    INFEASIBLE,
    OBLIGATORY,
    POSSIBLE,
    Design,
    classify_coupling,
    improves_on_baseline,
    rank_designs,
)


# ---- the 2026-08-07 regression: a design must BEAT the wild type, not inherit its coupling ----

def test_matching_the_wildtype_guarantee_is_NOT_a_design():
    """Anaerobic E. coli already secretes a little succinate obligatorily, so the wild type is
    OBLIGATORY before any edit. Counting 'is OBLIGATORY' as 'is a design' reported 2096 of 2096
    evaluated knockouts as designs -- every one merely inheriting the baseline."""
    assert not improves_on_baseline(0.005256, 0.005256)


def test_beating_the_wildtype_guarantee_IS_a_design():
    assert improves_on_baseline(4.0, 0.005256)


def test_a_worse_guarantee_is_not_a_design():
    assert not improves_on_baseline(0.001, 0.005256)


def test_baseline_improvement_needs_more_than_solver_noise():
    # an LP-tolerance-sized "gain" is not a gain
    assert not improves_on_baseline(0.005256 + 1e-9, 0.005256)


def test_improvement_from_an_uncoupled_baseline():
    # the ordinary case: wild type guaranteed nothing, the design guarantees something
    assert improves_on_baseline(2.0, 0.0)


# ---- pure: the coupling verdict (the load-bearing distinction of the whole module) ----

def test_obligatory_requires_a_POSITIVE_minimum():
    # min > 0 => every growing flux distribution secretes the product. That IS the design.
    assert classify_coupling(2.5, 9.0) == OBLIGATORY


def test_possible_is_not_a_design():
    # The cell CAN make it and CAN avoid it -- the un-engineered case. Must never read as coupled.
    assert classify_coupling(0.0, 9.0) == POSSIBLE


def test_numerically_zero_minimum_is_not_obligatory():
    # LP tolerances put "no flux" at ~1e-9, not exactly 0; that must not be sold as a guarantee.
    assert classify_coupling(1e-9, 9.0) == POSSIBLE
    assert classify_coupling(1e-12, 9.0) == POSSIBLE


def test_no_production_capacity_is_infeasible():
    assert classify_coupling(0.0, 0.0) == INFEASIBLE
    assert classify_coupling(0.0, 1e-9) == INFEASIBLE


def test_infeasible_LP_is_not_read_as_no_production():
    # None/NaN means "cannot grow at this floor", NOT "grows but makes nothing" -- either way not a design.
    assert classify_coupling(None, None) == INFEASIBLE
    assert classify_coupling(float("nan"), 5.0) == INFEASIBLE
    assert classify_coupling(1.0, float("nan")) == INFEASIBLE


@pytest.mark.parametrize("mn,mx,expected", [
    (5.0, 5.0, OBLIGATORY),   # forced to exactly one flux
    (-1.0, 5.0, POSSIBLE),    # negative min = net UPTAKE possible -> not guaranteed secretion
])
def test_coupling_edges(mn, mx, expected):
    assert classify_coupling(mn, mx) == expected


# ---- pure: ranking ----

def _d(kos, mn, mx=99.0, growth=0.5):
    return Design(knockouts=tuple(kos), min_product=mn, max_product=mx, growth=growth,
                  coupling=classify_coupling(mn, mx))


def test_ranking_puts_guaranteed_product_first():
    ranked = rank_designs([_d(["a"], 1.0), _d(["b"], 7.0), _d(["c"], 3.0)])
    assert [d.knockouts[0] for d in ranked] == ["b", "c", "a"]


def test_ranking_breaks_ties_toward_growth_then_fewer_knockouts():
    # same guaranteed product: faster-growing wins; then the cheaper (fewer-KO) strain wins.
    ranked = rank_designs([
        _d(["a", "b"], 5.0, growth=0.4),
        _d(["c"], 5.0, growth=0.4),
        _d(["d"], 5.0, growth=0.9),
    ])
    assert ranked[0].knockouts == ("d",)
    assert ranked[1].knockouts == ("c",)
    assert ranked[2].knockouts == ("a", "b")


def test_design_dict_reports_the_guarantee_and_the_count():
    rec = _d(["x", "y"], 2.0, 8.0, 0.3).as_dict()
    assert rec["n_knockouts"] == 2
    assert rec["min_product_flux"] == 2.0 and rec["coupling"] == OBLIGATORY


# ---- real model (slow): the honest baseline -- wild type must NOT look coupled ----

@pytest.mark.slow
def test_reproduces_the_literature_anaerobic_succinate_design():
    """The end-to-end gate: the search must FIND a published growth-coupled design.

    Anaerobic succinate via removal of the competing fermentation routes (PFL / LDH_D / ALCD2x) is the
    OptKnock-lineage design. Reaching it requires all three of the fixes made 2026-08-07 -- REACTION-level
    knockouts (GPR isozymes blunt gene-level), a NEAR-OPTIMAL growth floor taken against each strain's own
    maximum, and depth 3. Drop any one and this returns zero designs.
    """
    pytest.importorskip("cobra")
    from dna_decode.fba.design import find_coupled_designs
    from dna_decode.fba.model import load_model, wildtype_growth

    m = load_model()
    r = find_coupled_designs(
        m, "succ", wt_growth=wildtype_growth(m), anaerobic=True, level="reaction",
        max_knockouts=3, growth_frac=0.9, gene_ids=["PFL", "LDH_D", "ALCD2x", "ACKr", "PTAr"],
    )
    assert r["wildtype_already_coupled"] is True          # anaerobic WT already secretes a little
    assert r["n_coupled_designs"] >= 1                    # ... and the design must still BEAT it
    best = r["designs"][0]
    assert set(best["knockouts"]) == {"PFL", "LDH_D", "ALCD2x"}
    assert best["min_product_flux"] > 5.0                 # guaranteed succinate, not a rounding artifact
    assert best["improvement_over_wildtype"] > 5.0
    assert best["growth_per_h"] > 0.0                     # the designed strain must still grow


@pytest.mark.slow
def test_wildtype_ecoli_is_not_growth_coupled_for_succinate():
    """Real iML1515: WT can secrete succinate but is never forced to. A regression that reported the
    wild type as OBLIGATORY would make every design look unnecessary."""
    pytest.importorskip("cobra")
    from dna_decode.fba.design import evaluate_knockouts, resolve_target
    from dna_decode.fba.model import load_model, wildtype_growth

    m = load_model()
    wt = wildtype_growth(m)
    d = evaluate_knockouts(m, resolve_target(m, "succ"), [], 0.1, wt)
    assert d.coupling == POSSIBLE
    assert math.isclose(d.min_product, 0.0, abs_tol=1e-6)
    assert d.max_product > 1.0
