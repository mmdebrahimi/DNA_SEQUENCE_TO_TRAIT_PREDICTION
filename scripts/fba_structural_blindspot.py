"""Which genes can iML1515 NEVER call essential -- as a structural fact, not a per-run observation?

The flat-mechanism partition (`fba_flat_mechanism_partition.py`) classified a cell as "no flux" from ONE
returned optimal solution, which makes it basis-dependent and only a LOWER bound. This is the
basis-independent version, and it is model-only: no Fitness Browser DB required.

Two strictly nested classes, weakest first:

  BLOCKED   the gene's reactions cannot carry flux in ANY feasible solution, in any medium reachable by
            opening exchanges. A dead-end in the reconstruction. Such a gene can NEVER be called
            essential by any FBA variant, under any condition, with any constraint layer. Permanent.

  IDLE      unblocked, but zero flux is attainable at the optimum in this condition -- so deleting it
            leaves an optimal solution feasible and growth is unchanged. Condition-specific.

Why the distinction earns its keep: they imply DIFFERENT fixes. A BLOCKED gene needs the reconstruction
repaired (a missing transport/sink, a dead-end metabolite). An IDLE gene needs the OBJECTIVE to demand
its product -- a condition-specific biomass or maintenance term. Neither is fixed by tightening
constraints, which is why gap-fill, threshold retuning, pFBA and E-Flux all failed identically.

Conditions come from the COMMITTED artifact `wiki/fba_conditional_carbon_2026-08-13.json` (its
`conditions` map is name -> exchange id), so the 25-condition panel is reproducible with the Fitness
Browser DB offline. Only the experimental LABELS need that DB, and this script does not use them.

Usage:
    uv run python scripts/fba_structural_blindspot.py            # default medium only (fast)
    uv run python scripts/fba_structural_blindspot.py --panel    # all 25 carbon conditions
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.model import load_model  # noqa: E402

PANEL_ARTIFACT = Path("wiki/fba_conditional_carbon_2026-08-13.json")
FLUX_EPS = 1e-9


def panel_conditions() -> dict[str, str]:
    """{condition -> exchange} from the committed artifact. No feba.db needed."""
    d = json.loads(PANEL_ARTIFACT.read_text(encoding="utf-8"))
    conds = d.get("conditions")
    if not isinstance(conds, dict) or not conds:
        raise SystemExit(f"{PANEL_ARTIFACT} has no usable `conditions` map")
    return conds


def genes_fully_on(reaction_ids: set[str], model) -> list[str]:
    """Genes ALL of whose reactions lie in `reaction_ids`. A gene with no reactions is excluded."""
    out = []
    for g in model.genes:
        rxns = [r.id for r in g.reactions]
        if rxns and all(r in reaction_ids for r in rxns):
            out.append(g.id)
    return out


def apply_carbon(model, exchange: str, all_carbon: tuple[str, ...]) -> None:
    """Sole-carbon medium, mirroring `fitness_browser.apply_carbon_condition` without importing it
    (that module reaches for the Fitness Browser DB path at import time)."""
    medium = dict(model.medium)
    for ex in set(all_carbon) | {"EX_glc__D_e"}:
        medium.pop(ex, None)
    medium[exchange] = 10.0
    model.medium = medium


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", action="store_true", help="scan all 25 carbon conditions (slower)")
    ap.add_argument("--out", default=f"wiki/fba_structural_blindspot_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import find_blocked_reactions

    model = load_model()
    n_rxn, n_gene = len(model.reactions), len(model.genes)
    print(f"model {model.id}: {n_rxn} reactions, {n_gene} genes")

    # --- BLOCKED: open every exchange first, so "blocked" means structurally dead, not merely
    # starved by the current medium. This is the most generous medium the model admits.
    print("\nfinding structurally blocked reactions (all exchanges open) ...", flush=True)
    with model:
        for r in model.exchanges:
            r.lower_bound = min(r.lower_bound, -1000.0)
            r.upper_bound = max(r.upper_bound, 1000.0)
        blocked = set(find_blocked_reactions(model))
    blocked_genes = genes_fully_on(blocked, model)
    print(f"  blocked reactions : {len(blocked)}/{n_rxn} ({len(blocked)/n_rxn:.1%})")
    print(f"  genes entirely on blocked reactions: {len(blocked_genes)}/{n_gene} "
          f"({len(blocked_genes)/n_gene:.1%})  <- can NEVER be called essential")

    out = {
        "record": "fba-structural-blindspot-v1",
        "date": date.today().isoformat(),
        "model": model.id,
        "n_reactions": n_rxn,
        "n_genes": n_gene,
        "blocked": {
            "n_reactions": len(blocked),
            "n_genes_entirely_blocked": len(blocked_genes),
            "gene_fraction": round(len(blocked_genes) / n_gene, 4),
            "genes": sorted(blocked_genes),
            "blocked_reactions": sorted(blocked),
            "note": ("computed with EVERY exchange open, so these are dead-ends in the reconstruction "
                     "itself -- not an artefact of the medium. No FBA variant, condition or constraint "
                     "layer can make these genes essential."),
        },
        "method": "cobra.flux_analysis.find_blocked_reactions; conditions from the committed artifact",
        "needs_fitness_browser_db": False,
    }

    # --- IDLE per condition (optional, slower)
    if a.panel:
        conds = panel_conditions()
        all_ex = tuple(conds.values())
        print(f"\nper-condition idle scan over {len(conds)} conditions ...", flush=True)
        per_cond = {}
        for n, (name, ex) in enumerate(sorted(conds.items()), 1):
            with model:
                apply_carbon(model, ex, all_ex)
                sol = model.optimize()
                if sol.status != "optimal":
                    print(f"  [{n}/{len(conds)}] {name[:32]:34} NON-OPTIMAL")
                    continue
                fl = sol.fluxes
                idle_rxn = {r.id for r in model.reactions if abs(fl.get(r.id, 0.0)) <= FLUX_EPS}
            idle_genes = genes_fully_on(idle_rxn, model)
            only_idle = [g for g in idle_genes if g not in set(blocked_genes)]
            per_cond[name] = {"n_idle_genes": len(idle_genes),
                              "n_idle_not_blocked": len(only_idle)}
            print(f"  [{n}/{len(conds)}] {name[:32]:34} idle {len(idle_genes):4} "
                  f"(of which not-blocked {len(only_idle):4})", flush=True)
        out["idle_per_condition"] = per_cond
        if per_cond:
            vals = [v["n_idle_genes"] for v in per_cond.values()]
            out["idle_summary"] = {"min": min(vals), "max": max(vals),
                                   "mean": round(sum(vals) / len(vals), 1)}

    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
