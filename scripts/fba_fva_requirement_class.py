"""Does the model REQUIRE a gene's flux at the optimum, or is zero merely attainable?

This is the basis-independent successor to `fba_flat_mechanism_partition.py`. That script read ONE
optimal solution, so "carries no flux" depended on which vertex the LP happened to return. FVA at 100%
of optimum removes the dependence: it reports, per reaction, the FULL range of flux consistent with
optimal growth.

Per reaction, at 100% optimality:

  REQUIRED         0 is NOT in [min, max] -- every optimal solution routes flux here. Deleting the gene
                   MUST lower growth, so the model predicts essential.
  ZERO-ATTAINABLE  0 IS in [min, max] -- some optimal solution avoids the reaction entirely, so deleting
                   it leaves growth unchanged and the model predicts dispensable.
                   Split further by capability:
                     CAPABLE   max|flux| > 0 -- the reaction CAN carry flux, it just doesn't have to.
                               An OBJECTIVE problem: biomass does not demand the product.
                     INACTIVE_IN_CONDITION  max|flux| == 0 HERE. Pinned at zero in THIS medium.
                               **NOT the same as structurally blocked**: measured 2026-08-21, 57 of 58
                               such genes carry flux fine once other exchanges are opened. Only 1 is a
                               true reconstruction dead-end. So this class is medium-induced, and the
                               fix is a condition the model was never given -- not model repair.

The split matters because it names the fix. `fba_structural_blindspot.py` already showed BLOCKED
explains only 1 of 131 experimentally-essential genes, so the prediction is that ZERO-ATTAINABLE is
dominated by CAPABLE -- i.e. the deficit is about what the objective demands, not about dead-ends.
This script is the direct test of that.

**Gene-level caveat, stated up front:** FVA is per-REACTION. For a gene whose reactions are all
individually zero-attainable, it does not follow that they can be zero SIMULTANEOUSLY. So gene-level
ZERO-ATTAINABLE is necessary, not sufficient. The ground truth for a gene is the deletion itself, which
this script also runs for cross-check on the same genes.

Model-only: the essential-gene list comes from a COMMITTED artifact, so no Fitness Browser DB.

Usage:
    uv run python scripts/fba_fva_requirement_class.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.model import load_model  # noqa: E402

RATIOS_ARTIFACT = Path("wiki/fba_eflux_bridge_2026-08-17_ratios.json")
PANEL_ARTIFACT = Path("wiki/fba_conditional_carbon_2026-08-13.json")
EPS = 1e-9


def committed_essential_genes() -> list[str]:
    """The conditionally-essential gene ids, from a committed artifact (no feba.db)."""
    d = json.loads(RATIOS_ARTIFACT.read_text(encoding="utf-8"))
    genes: set[str] = set()
    for cells in d["arms"]["baseline"].values():
        genes |= set(cells)
    return sorted(genes)


def classify_reaction(lo: float, hi: float, eps: float = EPS) -> str:
    """REQUIRED / CAPABLE / BLOCKED for one reaction's FVA interval at optimum. PURE."""
    if lo > eps or hi < -eps:          # 0 outside [lo, hi] -> flux is forced
        return "REQUIRED"
    if abs(lo) <= eps and abs(hi) <= eps:
        # Pinned at zero IN THIS MEDIUM. Deliberately not called "blocked": cross-checking against
        # find_blocked_reactions with every exchange open showed 57 of 58 such genes are fine in a
        # richer medium. Medium-induced inactivity, not a reconstruction dead-end.
        return "INACTIVE_IN_CONDITION"
    return "CAPABLE"                   # zero attainable, but flux is possible


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default="D-Glucose")
    ap.add_argument("--out", default=f"wiki/fba_fva_requirement_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion

    model = load_model()
    conds = json.loads(PANEL_ARTIFACT.read_text(encoding="utf-8"))["conditions"]
    if a.condition not in conds:
        raise SystemExit(f"unknown condition {a.condition!r}; have {sorted(conds)[:5]} ...")
    carbon = conds[a.condition]
    all_ex = tuple(conds.values())

    genes = committed_essential_genes()
    print(f"committed conditionally-essential genes: {len(genes)}")
    print(f"condition: {a.condition} ({carbon})")

    with model:
        medium = dict(model.medium)
        for ex in set(all_ex) | {"EX_glc__D_e"}:
            medium.pop(ex, None)
        medium[carbon] = 10.0
        model.medium = medium

        wt = model.slim_optimize()
        print(f"wildtype growth: {wt:.6f}")

        rxns = sorted({r for g in genes for r in model.genes.get_by_id(g).reactions},
                      key=lambda r: r.id)
        print(f"reactions carried by those genes: {len(rxns)}")
        print("running FVA at 100% of optimum ...", flush=True)
        fva = flux_variability_analysis(model, reaction_list=rxns, fraction_of_optimum=1.0)

        # gene-level ground truth on the SAME genes, for cross-check
        print("running gene deletions for cross-check ...", flush=True)
        dele = single_gene_deletion(
            model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
    del_ratio = {}
    for _, row in dele.iterrows():
        gid = next(iter(row["ids"]))
        g = row["growth"]
        del_ratio[gid] = 0.0 if g != g else float(g) / wt

    rxn_class = {r.id: classify_reaction(float(fva.loc[r.id, "minimum"]),
                                         float(fva.loc[r.id, "maximum"])) for r in rxns}

    gene_class: dict[str, str] = {}
    for gid in genes:
        cls = {rxn_class[r.id] for r in model.genes.get_by_id(gid).reactions}
        if "REQUIRED" in cls:
            gene_class[gid] = "REQUIRED"
        elif cls == {"INACTIVE_IN_CONDITION"}:
            gene_class[gid] = "INACTIVE_IN_CONDITION"
        else:
            gene_class[gid] = "CAPABLE_BUT_IDLE"

    counts: dict[str, int] = {}
    for v in gene_class.values():
        counts[v] = counts.get(v, 0) + 1

    # cross-check: REQUIRED should predict essential (ratio < 1); the others ratio == 1
    mism = [(g, gene_class[g], round(del_ratio.get(g, float("nan")), 4))
            for g in genes
            if (gene_class[g] == "REQUIRED") != (del_ratio.get(g, 1.0) < 1 - 1e-6)]

    total = len(genes)
    print(f"\ngene classes over {total} experimentally-essential genes ({a.condition}):")
    for k in ("REQUIRED", "CAPABLE_BUT_IDLE", "INACTIVE_IN_CONDITION"):
        n = counts.get(k, 0)
        print(f"  {k:18} {n:4} ({n/total:.1%})")
    print(f"\nFVA-vs-deletion disagreements: {len(mism)} (expected 0 for single-reaction genes; "
          f"multi-reaction genes can disagree -- FVA is per-reaction, joint zero is not implied)")
    for row in mism[:5]:
        print(f"   {row}")

    n_or = sum(1 for g, _c, _r in mism
               if any(" or " in (r.gene_reaction_rule or "").lower()
                      for r in model.genes.get_by_id(g).reactions))
    print(f"  of those disagreements, {n_or} carry an isozyme OR in their GPR -- the reaction IS",
          "required, but a paralog covers the deletion, so the GENE is not essential")

    out = {
        "record": "fba-fva-requirement-class-v1",
        "date": date.today().isoformat(),
        "model": model.id,
        "condition": a.condition,
        "carbon_exchange": carbon,
        "wildtype_growth": round(float(wt), 6),
        "n_genes": total,
        "counts": counts,
        "fractions": {k: round(v / total, 4) for k, v in counts.items()},
        "n_fva_vs_deletion_disagreements": len(mism),
        "n_disagreements_with_isozyme_or": n_or,
        "disagreement_examples": mism[:20],
        "gene_class": gene_class,
        "caveats": [
            "FVA is per-REACTION. Gene-level CAPABLE_BUT_IDLE means every reaction can individually "
            "reach zero at optimum; it does NOT prove they can all be zero simultaneously. The "
            "deletion cross-check is the gene-level ground truth and is reported alongside.",
            "One condition at a time. A gene idle on this carbon source may be REQUIRED on another.",
            "The essential-gene list is the committed 131-gene 11-condition E-Flux panel, not the "
            "full 217-gene set (which needs feba.db on the disconnected D:).",
            "INACTIVE_IN_CONDITION is medium-induced, NOT structural: 57 of 58 such genes are unblocked "
            "once other exchanges are opened (cross-checked against fba_structural_blindspot).",
        ],
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
