"""Kill-test for the class-A claim: a zero-flux gene is UNREACHABLE by any constraint-based method.

THE CLAIM UNDER TEST
    If some optimal flux distribution carries zero flux through every reaction of gene g, then deleting g
    leaves the optimal objective EXACTLY unchanged (ratio 1.0) -- and no tightening of reaction bounds
    (E-Flux, pFBA, regulatory restriction, ...) can change that, because such methods only reshape flux
    among reactions that carry it.

WHY IT MATTERS
    Four independent levers (gap-fill, threshold retune, pFBA, E-Flux) all failed the same way on
    conditional gene essentiality. If this claim holds, it is a MECHANISM for that shared failure rather
    than four coincidences: a measured 25.1% of true-essential gene x condition cells are zero-flux, and
    every one of them is unreachable by construction.

HOW THIS COULD FAIL (it is a real test, not a tautology)
    - Deleting a zero-flux gene could still lower growth, if the gene participates in a reaction that is
      zero in the returned solution but load-bearing in every alternative optimum -- i.e. if the LP's
      returned solution is not representative.
    - Bound-tightening could make a previously-idle reaction essential, if the tightening forces flux
      onto the idle route. If that happens the "unreachable" half is false.

These run on the MODEL ALONE (cobrapy cache is local) and need no Fitness Browser DB.
"""
from __future__ import annotations

import pytest

pytest.importorskip("cobra", reason="cobrapy required")

FLUX_EPS = 1e-9
TOL = 1e-6


@pytest.fixture(scope="module")
def model():
    from dna_decode.fba.model import load_model
    return load_model()


def _zero_flux_genes(model, fluxes, limit=40):
    """Genes whose EVERY reaction is at zero flux in this optimal solution."""
    out = []
    for g in model.genes:
        rxns = list(g.reactions)
        if not rxns:
            continue
        if all(abs(fluxes.get(r.id, 0.0)) <= FLUX_EPS for r in rxns):
            out.append(g.id)
        if len(out) >= limit:
            break
    return out


def test_zero_flux_gene_deletion_leaves_growth_exactly_unchanged(model):
    """The theorem's first half, on real genes of the real model."""
    from cobra.flux_analysis import single_gene_deletion

    with model:
        sol = model.optimize()
        assert sol.status == "optimal"
        wt = sol.objective_value
        fluxes = dict(sol.fluxes)
        genes = _zero_flux_genes(model, fluxes)
        assert genes, "no zero-flux genes found -- the premise itself would be false"
        res = single_gene_deletion(
            model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)

    offenders = []
    for _, row in res.iterrows():
        gid = next(iter(row["ids"]))
        g = row["growth"]
        if g != g:                      # NaN => infeasible => a real change
            offenders.append((gid, "infeasible"))
        elif abs(g - wt) > TOL:
            offenders.append((gid, g))
    assert not offenders, f"zero-flux deletion CHANGED growth (theorem false): {offenders[:5]}"


def test_bound_tightening_cannot_make_a_zero_flux_gene_essential(model):
    """The theorem's second half -- the part that explains the four failed levers.

    Emulates what E-Flux/pFBA-style methods do: scale every GPR-carrying reaction's bounds toward zero.
    A zero-flux gene must STILL be non-essential afterwards, because the constraint only removes capacity
    from routes, and this gene's routes were already carrying nothing.
    """
    from cobra.flux_analysis import single_gene_deletion

    with model:
        sol = model.optimize()
        wt_before = sol.objective_value
        genes = _zero_flux_genes(model, dict(sol.fluxes))
        assert genes

        exchanges = {r.id for r in model.exchanges}
        for rxn in model.reactions:
            if rxn.id in exchanges or not rxn.gene_reaction_rule.strip():
                continue
            # 0.05, NOT 0.5. MEASURED 2026-08-20: iML1515's internal bounds are entirely
            # non-binding down to a 10x tightening -- growth is bit-identical at scale 0.5, 0.2 and
            # 0.1, and only responds at 0.05 (0.877 -> 0.715). A 2x tightening is INERT, so the first
            # version of this test was vacuous and its own guard caught it.
            if rxn.upper_bound > 0:
                rxn.upper_bound *= 0.05
            if rxn.lower_bound < 0:
                rxn.lower_bound *= 0.05

        sol2 = model.optimize()
        assert sol2.status == "optimal", "tightening made the model infeasible -- test inconclusive"
        wt_after = sol2.objective_value
        # the tightening must actually BITE, or the test proves nothing
        assert wt_after < wt_before - TOL, (
            f"bound tightening did not reduce growth ({wt_before} -> {wt_after}); "
            "the perturbation was inert so this test would be vacuous")

        res = single_gene_deletion(
            model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)

    became_essential = []
    for _, row in res.iterrows():
        gid = next(iter(row["ids"]))
        g = row["growth"]
        if g != g or g < 0.99 * wt_after:
            became_essential.append((gid, g))
    assert not became_essential, (
        "bound tightening MADE a zero-flux gene matter -- the 'unreachable by constraints' half is "
        f"false: {became_essential[:5]}")
