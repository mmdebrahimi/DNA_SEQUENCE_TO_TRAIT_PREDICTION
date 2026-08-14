"""KILL-TEST for the G9 abduction: are flat genes flat because their products are not biomass precursors?

    uv run python scripts/fba_flat_gene_biomass_killtest.py

**Claim under test (candidate C-G9-1):** the conditionally-essential genes that iML1515 keeps FLAT
(deletion changes nothing) and that are NOT explained by isozyme GPR redundancy are flat because their
reaction products are **not biomass precursors** — the objective never demands what they make, so
deleting them is free by construction.

This is the best available explanation for the surprising residual: GPR isozyme structure covers only
32 of 110 flat genes (29%), leaving ~78 flat for reasons above the GPR
(`wiki/fba_constant_gene_diagnostic_2026-08-13.md`).

**POLARITY (falsification-engine convention): exit 0 = DISPROVED = the claim is KILLED.**
The claim predicts flat genes are ENRICHED for non-biomass-precursor status relative to the
conditionally-essential genes the model does NOT keep flat. Exit 0 (kill) when that enrichment is absent
or reversed; exit 1 (survive) when it holds.

Threshold, derived not asserted: the claim must explain a 71%-unexplained residual, so a difference below
10 percentage points would not do the work being claimed. `MIN_ENRICHMENT_PP = 10.0`.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.conditional_essentiality import conditionally_essential_genes  # noqa: E402
from dna_decode.fba.fitness_browser import (  # noqa: E402
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402

FRAC = 0.01
FLAT_EPS = 1e-6
MIN_ENRICHMENT_PP = 10.0


def biomass_precursors(model) -> set[str]:
    """Metabolite ids the objective (biomass) reaction CONSUMES. PURE."""
    obj = next((r for r in model.reactions if r.objective_coefficient), None)
    if obj is None:
        obj = next((r for r in model.reactions if "biomass" in r.id.lower()), None)
    if obj is None:
        raise SystemExit("no objective/biomass reaction found -- cannot evaluate the claim")
    return {m.id for m, c in obj.metabolites.items() if c < 0}


def makes_precursor(model, gene_id: str, precursors: set[str]) -> bool:
    """Does ANY reaction this gene participates in touch a biomass precursor? PURE.

    Deliberately GENEROUS (any participation, either direction): a generous definition makes the claim
    HARDER to survive, which is the right bias for a falsification attempt.
    """
    try:
        gene = model.genes.get_by_id(gene_id)
    except KeyError:
        return False
    return any(m.id in precursors for rxn in gene.reactions for m in rxn.metabolites)


def main() -> int:
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    try:
        conn = open_db()
    except FileNotFoundError as e:
        print(f"INDETERMINATE: {e}", file=sys.stderr)
        return 2                      # neither killed nor survived

    model = load_model()
    conds = carbon_conditions(conn, model)
    keys = tuple(sorted(conds))
    subset = conditionally_essential_genes(
        load_records(conn, conds, gene_filter={g.id for g in model.genes}))
    genes = [r.gene_id for r in subset]
    print(f"{model.id}: {len(keys)} conditions | {len(genes)} conditionally-essential genes", flush=True)

    ratios: dict[str, dict[str, float]] = {}
    all_ex = tuple(conds.values())
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            rat = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes])
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    gv = row["growth"]
                    rat[gid] = None if gv != gv else float(gv) / wt
            ratios[cond] = rat
        print(f"   [{n:2d}/{len(keys)}] {cond[:34]:36s}", flush=True)

    # FLAT = ratio ~1.0 in EVERY condition where the gene is truly essential (the diagnostic's definition)
    flat, not_flat = [], []
    for r in subset:
        true_ess = [c for c in keys if r.experimental.get(c, False)]
        if not true_ess:
            continue
        vals = [ratios.get(c, {}).get(r.gene_id) for c in true_ess]
        if all(v is not None and v >= 1.0 - FLAT_EPS for v in vals):
            flat.append(r.gene_id)
        else:
            not_flat.append(r.gene_id)

    prec = biomass_precursors(model)
    print(f"\nbiomass precursors: {len(prec)} metabolites")
    print(f"flat genes: {len(flat)} | non-flat conditionally-essential: {len(not_flat)}")
    if not flat or not not_flat:
        print("INDETERMINATE: one arm is empty -- no comparison possible", file=sys.stderr)
        return 2

    flat_no_prec = sum(1 for g in flat if not makes_precursor(model, g, prec)) / len(flat)
    nf_no_prec = sum(1 for g in not_flat if not makes_precursor(model, g, prec)) / len(not_flat)
    enrichment_pp = 100.0 * (flat_no_prec - nf_no_prec)

    print(f"\n   flat genes NOT touching a biomass precursor : {100 * flat_no_prec:.1f}%")
    print(f"   non-flat genes NOT touching a precursor      : {100 * nf_no_prec:.1f}%")
    print(f"   enrichment                                   : {enrichment_pp:+.1f} pp "
          f"(claim needs >= +{MIN_ENRICHMENT_PP})")

    if enrichment_pp >= MIN_ENRICHMENT_PP:
        print("\nCLAIM SURVIVES this falsification attempt (enrichment holds).")
        return 1
    print("\nCLAIM KILLED: flat genes are NOT meaningfully enriched for non-precursor status, "
          "so 'the objective never demands what they make' does not explain the flatness.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
