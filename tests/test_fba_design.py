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

def test_default_pool_strategy_does_not_narrow_the_search():
    """The default must stay the exhaustive strategy, not the MEASURED-NEGATIVE `competition` heuristic.

    `competition_ranking` was built 2026-08-07 to fix the pool limitation and FAILED twice (see its
    docstring): as written it ranks by "is growth-associated" and floods the top with ion transporters
    and BIOMASS; matched-growth fixes that but drops PFL/ALCD2x to rank ~1450. Defaulting to it would
    narrow 2266 candidates to 120 on an UNVALIDATED score -- trading an honest zero for a fast one that
    could hide designs the exhaustive scan would reach.
    """
    import inspect

    from dna_decode.fba.design import find_coupled_designs

    default = inspect.signature(find_coupled_designs).parameters["pool_strategy"].default
    assert default == "single_effect", (
        f"default pool_strategy is {default!r}; 'competition' is a measured negative and must stay opt-in"
    )


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
def test_milp_finds_the_same_design_and_the_same_guarantee():
    """The MILP path must reach the SAME design as the enumeration, with the SAME verified guarantee.

    Two independent methods agreeing on 9.263835 is the cross-validation that makes either trustworthy.
    Also pins the floor semantics that cost a 30-minute run: the MILP floor is a fraction of WILD-TYPE
    growth, and at 0.9 the known design is `infeasible` by construction because it only grows at 0.0816
    (52% of the anaerobic wild type). Skips cleanly when the optional [design] extra is absent.
    """
    pytest.importorskip("cobra")
    pytest.importorskip("straindesign")
    from dna_decode.fba.design import find_coupled_designs_milp
    from dna_decode.fba.model import load_model

    m = load_model()
    r = find_coupled_designs_milp(
        m, "succ", anaerobic=True, growth_floor_frac_of_wt=0.5, max_knockouts=3,
        candidates=["PFL", "LDH_D", "ALCD2x", "ACKr", "PTAr", "FRD2", "FRD3", "PPCK", "PPC"],
        time_limit=600,
    )
    assert r["n_designs"] >= 1, f"MILP found nothing (status {r['status']})"
    best = r["designs"][0]
    assert set(best["knockouts"]) == {"PFL", "LDH_D", "ALCD2x"}
    # the guarantee is re-derived by evaluate_knockouts, so it must match the enumeration path exactly
    assert best["min_product_flux"] == pytest.approx(9.263835, abs=1e-3)
    assert best["coupling"] == OBLIGATORY


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


@pytest.mark.slow
def test_restricted_pool_finds_anaerobic_d_lactate_designs():
    """The companion gate to the succinate one: the POOL, not the budget, is what finds designs.

    Measured 2026-11-19 after four unrestricted arms found nothing overnight. Searching all 2266
    candidate reactions returns ZERO coupled D-lactate designs at every depth and budget tried --
    including a bilevel MILP given 12600 s per arm, seven times its hardcoded limit. Restricting the
    pool to the routes that COMPETE with lactate for fermentative carbon and redox returns 63
    designs exhaustively in seconds.

    Why the unrestricted search cannot find them, measured: PFL alone raises the product CEILING
    (7.256 -> 18.458) and leaves the FLOOR at 0.0, so it scores improvement_over_wildtype exactly
    0.0. Coupling first appears at DEPTH 2. Since the pair/triple pool is PRE-RANKED, a member that
    individually scores zero can never be promoted, and the triple containing it is never built --
    which is verbatim the limitation pyproject.toml records for this path.

    LDH_D is DELIBERATELY ABSENT from the pool: it is the lactate-PRODUCING enzyme, and the
    succinate gate above reaches ITS product precisely by deleting it. The two products compete.
    """
    pytest.importorskip("cobra")
    from dna_decode.fba.design import find_coupled_designs
    from dna_decode.fba.model import load_model, wildtype_growth

    pool = ["PFL", "ALCD2x", "ACALD", "FRD2", "FRD3", "ACKr", "PTAr", "PTA2",
            "POX", "MGSA", "LDH_D2", "L_LACD2", "L_LACD3"]
    assert "LDH_D" not in pool, "deleting the producing enzyme cannot yield a lactate design"

    m = load_model()
    r = find_coupled_designs(
        m, "EX_lac__D_e", wt_growth=wildtype_growth(m), anaerobic=True, level="reaction",
        max_knockouts=3, growth_frac=0.9, gene_ids=pool,
    )

    assert r["n_candidates_scanned"] == len(pool)
    assert r["wildtype_already_coupled"] is False       # unlike succinate, WT secretes none
    assert r["n_coupled_designs"] >= 10                 # 63 measured; the bar is deliberately loose
    best = r["designs"][0]
    assert len(best["knockouts"]) >= 2                  # depth 1 cannot couple -- see the docstring
    assert "PFL" in best["knockouts"]                   # the zero-scoring member every design needs
    assert best["min_product_flux"] > 5.0               # guaranteed lactate, not a rounding artifact
    assert best["growth_per_h"] > 0.0                   # the designed strain must still grow
    assert best["min_product_flux"] < 2 * 10.0          # 2 lactate per glucose at an uptake of 10


@pytest.mark.slow
def test_a_single_knockout_cannot_couple_lactate_and_that_is_why_pre_ranking_fails():
    """The mechanism behind the test above, pinned on its own.

    This is the whole reason an unrestricted search fails, so it is worth a test that does not
    depend on the search at all: PFL alone is the highest-value MEMBER of every top design and
    scores ZERO improvement by itself.
    """
    pytest.importorskip("cobra")
    from dna_decode.fba.design import evaluate_knockouts, set_anaerobic
    from dna_decode.fba.model import load_model

    # The `with m:` block is LOAD-BEARING, not hygiene, and `reuse=False` does NOT substitute for it.
    # `load_model()` hands back a SHARED CACHED object, and `set_anaerobic` mutates bounds in place.
    # The first version of this test called `set_anaerobic` on that shared object and left EX_o2_e at
    # 0.0 for every later caller in the process, which broke test_fba_zero_flux_theorem and
    # test_fba_gapfill in the FULL suite while both still passed in isolation.
    #
    # `reuse=False` looks like the fix and is NOT: model.py writes `_MODELS[mid] = model`
    # unconditionally, guarded only on `path is None` and never on `reuse`, so a reuse=False load
    # REPLACES the cache with its fresh object and mutating that pollutes exactly the same way
    # (measured). The module's own idiom is the real fix -- `evaluate_knockouts` and
    # `find_coupled_designs` both confine their condition changes to `with model:`, which is why the
    # sibling test above never polluted anything.
    m = load_model()
    with m:
        set_anaerobic(m)
        wt = m.slim_optimize()
        tgt = m.reactions.get_by_id("EX_lac__D_e")

        alone = evaluate_knockouts(m, tgt, ["PFL"], 0.9, wt, "reaction")
        pair = evaluate_knockouts(m, tgt, ["PFL", "ALCD2x"], 0.9, wt, "reaction")

    # PFL alone lifts the ceiling but NOT the floor -- invisible to any improvement-based ranking.
    assert alone.coupling == POSSIBLE
    assert math.isclose(alone.min_product, 0.0, abs_tol=1e-6)
    assert alone.max_product > 15.0

    # Add one more and the floor leaves zero. Coupling is a property of the SET, not its members.
    assert pair.coupling == OBLIGATORY
    assert pair.min_product > 5.0

    # The condition change must NOT outlive this test -- this is the regression guard.
    assert m.reactions.get_by_id("EX_o2_e").lower_bound < 0.0
