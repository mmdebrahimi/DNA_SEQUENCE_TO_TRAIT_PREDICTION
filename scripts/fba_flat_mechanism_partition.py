"""WHY does the model miss a conditionally-essential gene? Partition every miss by mechanism.

Four independent levers (gap-fill, threshold retune, pFBA route restriction, E-Flux expression
constraints) all failed the SAME way, and the shared diagnosis was that 70-95% of missed essential cells
are FLAT -- the in-silico deletion changed nothing. "Flat" is a symptom, not a mechanism. This partitions
it into three mutually-exclusive, mechanistically distinct classes:

  A. the gene's reaction carries NO FLUX at the optimum in that condition
  B. it carries flux AND the GPR has an `or` (isozyme redundancy)
  C. it carries flux, no `or`

**CLASS C IS NOT A FAILURE CLASS -- do not describe it as "compensated by rerouting".** The denominator
here is ALL true-essential cells, hits included. Class A is 100% missed by the theorem below, so every
one of the model's TRUE POSITIVES lives in B or C. Class C therefore mixes cells the model got RIGHT
with cells it missed, and the split is NOT measured by this script (it records mechanism class, not
`predicted_essential`). Bounding it from the committed carbon panel: class-C MISSES lie somewhere in
[112, 624] of the 1,083 misses. Join the hit/miss column before treating C as a mechanism at all.

The classes matter because they imply DIFFERENT fixes, and one of them is a proof:

  **A class-A cell is NECESSARILY missed -- this is a theorem, not a correlation.** If optimal solution
  S carries zero flux through reaction R, then deleting R's gene (forcing R to 0) leaves S feasible and
  still optimal, so the deletion objective EQUALS the wildtype objective and the ratio is exactly 1.0.
  Within a FIXED objective and a FIXED medium, no intervention that merely tightens or removes capacity
  on reactions carrying no flux can move such a cell -- that class of intervention only reshapes flux
  among reactions that already carry it. That is a mechanism for why the four levers failed identically.

  **Do NOT broaden this to "no constraint-based method can move it."** Opening an exchange IS a
  constraint change, and `wiki/fba_structural_blindspot_2026-08-21.md` measured that 57 of 58 genes
  pinned at zero flux in glucose carry flux fine once other exchanges are opened. Changing the medium,
  the objective, or adding a demand term can all remove the zero-flux optimum.

  Because the proof needs only SOME zero-flux optimum, class A is a LOWER BOUND on the
  necessarily-missed set: a cell in class B or C may also admit an alternative optimum with zero flux.

**Evaluate each cell IN ITS OWN CONDITION.** A glucose-only reference state reported class A at 50.2%;
the correct per-condition measure is 25.1%. A gene can be idle on glucose and load-bearing on xylose, so
the reference-state shortcut overstates class A ~2x.

Usage:
    uv run python scripts/fba_flat_mechanism_partition.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.conditional_essentiality import conditionally_essential_genes  # noqa: E402
from dna_decode.fba.fitness_browser import (  # noqa: E402
    ESSENTIAL_FITNESS,
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model  # noqa: E402

FLUX_EPS = 1e-9

SCOPE_NOTE = (
    "Scope, stated because it is easy to overclaim and my own data refutes the broad version: this holds with the OBJECTIVE and the MEDIUM HELD FIXED, against interventions that only tighten or remove capacity on reactions not used by at least one optimal solution. It does NOT say the cell is beyond every constraint-based method -- opening exchanges is a constraint change, and wiki/fba_structural_blindspot_2026-08-21.md measured that 57 of 58 zero-flux-in-glucose genes carry flux fine in a richer medium."
)


def classify_cell(reactions, gpr_has_or: bool, fluxes, eps: float = FLUX_EPS) -> str:
    """Which mechanism class does this (gene, condition) cell fall into? PURE.

    Order matters: no-flux is checked FIRST, because a gene whose reactions are all idle is class A
    regardless of whether its GPR also happens to carry an `or`.
    """
    carries = any(abs(fluxes.get(r, 0.0)) > eps for r in reactions)
    if not carries:
        return "no_flux"
    return "isozyme_or" if gpr_has_or else "reroute"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=f"wiki/fba_flat_mechanism_partition_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    model = load_model()
    conn = open_db()
    conds = carbon_conditions(conn, model)
    keys = sorted(conds)
    records = load_records(conn, conds, gene_filter={g.id for g in model.genes},
                           threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    print(f"conditions {len(keys)} | conditionally-essential genes {len(subset)}")

    rxn_ids = {r.gene_id: [x.id for x in model.genes.get_by_id(r.gene_id).reactions] for r in subset}
    has_or = {r.gene_id: any(" or " in (x.gene_reaction_rule or "").lower()
                             for x in model.genes.get_by_id(r.gene_id).reactions) for r in subset}

    counts = {"no_flux": 0, "isozyme_or": 0, "reroute": 0}
    per_cond = {}
    all_ex = tuple(conds.values())
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            sol = model.optimize()
            if sol.status != "optimal":
                print(f"  [{n}/{len(keys)}] {cond[:34]:36} NON-OPTIMAL -- skipped")
                continue
            fluxes = dict(sol.fluxes)
        local = {"no_flux": 0, "isozyme_or": 0, "reroute": 0}
        for r in subset:
            if not r.experimental.get(cond):
                continue  # only TRUE-essential cells
            cls = classify_cell(rxn_ids[r.gene_id], has_or[r.gene_id], fluxes)
            counts[cls] += 1
            local[cls] += 1
        per_cond[cond] = local
        print(f"  [{n}/{len(keys)}] {cond[:34]:36} {local}", flush=True)

    total = sum(counts.values())
    if not total:
        print("no true-essential cells -- nothing to partition", file=sys.stderr)
        return 2

    out = {
        "record": "fba-flat-mechanism-partition-v1",
        "date": date.today().isoformat(),
        "model": model.id,
        "labels": f"Fitness Browser RB-TnSeq orgId=Keio, carbon source, fit<{a.threshold}",
        "n_conditions": len(keys),
        "n_conditionally_essential_genes": len(subset),
        "n_true_essential_cells": total,
        "counts": counts,
        "fractions": {k: round(v / total, 4) for k, v in counts.items()},
        "per_condition": per_cond,
        "class_A_is_necessarily_missed": True,
        "class_A_note": ("a reaction carrying no flux at the optimum cannot lower the optimal objective "
                         "when deleted -- the same flux distribution stays feasible -- so the ratio is "
                         "1.0 and the gene is predicted dispensable. " + SCOPE_NOTE + " "
                         "these cells."),
        "caveats": [
            "Each cell is evaluated IN ITS OWN CONDITION. A glucose-only reference state gives 50.2% "
            "for class A vs the correct 25.1% -- the shortcut overstates it ~2x.",
            "Denominator is TRUE-essential cells (both hits and misses), not misses alone -- so "
            "classes B and C each contain the model's true positives as well as its misses, and this "
            "script does NOT separate them.",
            "The A/C boundary is basis-dependent: `classify_cell` reads ONE optimal solution, so a cell "
            "called C may have an alternative optimum with zero flux and really be A. FVA-at-optimum "
            "(is max|v| == 0?) is the basis-independent version and is NOT done here.",
            "Class A is a LOWER BOUND on the necessarily-missed set, not an estimate of it. The "
            "proof only needs SOME optimal solution with zero flux through the reaction, and the "
            "solver returned one -- so every cell classified A is genuinely necessarily-missed. But a "
            "cell classified B or C may ALSO admit an alternative optimum with zero flux, which would "
            "make it necessarily-missed too. FVA (is zero flux attainable at optimum?) would give the "
            "true set and is NOT done here.",
        ],
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nTRUE-essential cells: {total}")
    for k, v in counts.items():
        print(f"  {k:12} {v:5} ({v/total:.1%})")
    print(f"wrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
