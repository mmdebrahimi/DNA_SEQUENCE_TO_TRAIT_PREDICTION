"""The remaining 100: does the model NEVER fire for them, or does it fire in the WRONG condition?

WHERE THIS SITS
The deficit accounting so far, over the 131 conditionally-essential genes:

  31   PROVABLY UNCALLABLE   no medium, objective or constraint layer can move them
                             (`wiki/fba_replaceability_lp_2026-08-22.md`)
   8   ...of which isozyme-masked, and expression-gating does NOT recover them
                             (`wiki/fba_expression_gated_gpr_result_2026-08-22.md`)
 100   CALLABLE, STILL MISSED   <- this script

"Callable" is a statement about the model's REACH, not about whether it ever actually fires. So the 100
split two ways, and the two imply completely different work:

  MIS-CONDITIONED   the model DOES predict this gene essential -- just not in the condition where the
                    experiment says it matters. The capability is there and the mapping is wrong.
                    Fix lives in condition modelling (medium composition, aeration, uptake rates).

  NEVER-FIRES       the model never predicts it essential in ANY of the 25 carbon conditions, despite
                    being structurally capable of it. Nothing in the objective ever demands its product.
                    Fix lives in the biomass objective / network realism.

The split is decisive because MIS-CONDITIONED is a *calibration* problem and NEVER-FIRES is a *modelling*
problem, and the four failed levers were aimed at neither.

GENEROSITY IS DELIBERATE. The "ever fires" test runs over all 25 carbon conditions with labels, not just
the 11 that also have expression data -- giving the model the maximum number of chances to fire. That
makes NEVER-FIRES a conservative claim: a gene in that class failed to fire anywhere we could look.

Model + `feba.db` (labels). No expression data needed.

Usage:
    uv run python scripts/fba_missed_gene_partition.py
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
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402
from dna_decode.fba.nitrogen import (  # noqa: E402
    apply_nitrogen_condition,
    load_nitrogen_records,
    nitrogen_conditions,
)

#: The axis seam. Every result in this arc was measured on CARBON only; `MIS_CONDITIONED = 0` is
#: therefore a ONE-AXIS result, which is exactly what the nitrogen arm tests. The three callables per
#: axis are drop-in mirrors of each other (`nitrogen.py` documents the parallel deliberately).
AXES = {
    "carbon": (carbon_conditions, apply_carbon_condition, load_records, "all_carbon"),
    "nitrogen": (nitrogen_conditions, apply_nitrogen_condition, load_nitrogen_records, "all_nitrogen"),
}

FRAC = 0.01
SCREEN_ARTIFACT = Path("wiki/fba_orphan_protection_2026-08-21.json")
REPL_ARTIFACT = Path("wiki/fba_replaceability_2026-08-22.json")


def uncallable_set() -> set[str]:
    """The 31 proved-unreachable genes: the exact-duplicate screen UNION the replaceability LP."""
    a = set(json.loads(SCREEN_ARTIFACT.read_text(encoding="utf-8"))
            ["structurally_uncallable"]["genes"])
    b = set(json.loads(REPL_ARTIFACT.read_text(encoding="utf-8"))["loop_audit"]["genes"])
    return a | b


def classify(pred_conds: set[str], true_conds: set[str]) -> str:
    """PURE. How does the model's firing pattern relate to the experiment's?"""
    if not pred_conds:
        return "NEVER_FIRES"
    if pred_conds & true_conds:
        return "PARTIAL_OVERLAP"          # fires in at least one right place, misses others
    return "MIS_CONDITIONED"              # fires, but never where the experiment says it matters


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--axis", choices=sorted(AXES), default="carbon",
                    help="environmental axis; MIS_CONDITIONED=0 was only ever measured on carbon")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    if a.out is None:
        suffix = "" if a.axis == "carbon" else f"_{a.axis}"
        a.out = f"wiki/fba_missed_partition{suffix}_{date.today().isoformat()}"

    from cobra.flux_analysis import single_gene_deletion

    cond_fn, apply_fn, load_fn, all_kw = AXES[a.axis]
    model = load_model()
    conn = open_db()
    conds = cond_fn(conn, model)
    all_ex = tuple(conds.values())
    keys = sorted(conds)
    model_genes = {g.id for g in model.genes}
    print(f"axis: {a.axis}")
    records = load_fn(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    genes = [r.gene_id for r in subset]
    truth = {r.gene_id: {c for c, y in r.experimental.items() if y} for r in subset}
    names = {r.gene_id: r.gene for r in subset}
    print(f"conditions with labels: {len(keys)}")
    print(f"conditionally essential genes: {len(genes)}")

    unc = uncallable_set() & set(genes)
    print(f"of which PROVABLY UNCALLABLE (prior work): {len(unc)}")

    # --- model predictions across every labelled condition
    pred: dict[str, set[str]] = {g: set() for g in genes}
    for n, cond in enumerate(keys, 1):
        with model:
            apply_fn(model, conds[cond], **{all_kw: all_ex})
            wt = wildtype_growth(model)
            if not (wt > 1e-9):
                print(f"  [{n}/{len(keys)}] {cond[:34]:36} WILDTYPE INFEASIBLE -- skipped")
                continue
            res = single_gene_deletion(
                model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
            hit = 0
            for _, row in res.iterrows():
                gid = next(iter(row["ids"]))
                g = row["growth"]
                if (g != g) or (g < FRAC * wt):
                    pred[gid].add(cond)
                    hit += 1
        print(f"  [{n}/{len(keys)}] {cond[:34]:36} wt={float(wt):.5f} predicted essential={hit}",
              flush=True)

    rows, counts = [], {}
    for g in genes:
        if g in unc:
            cls = "PROVABLY_UNCALLABLE"
        else:
            cls = classify(pred[g], truth[g])
        counts[cls] = counts.get(cls, 0) + 1
        rows.append({"gene": g, "name": names.get(g, ""), "class": cls,
                     "n_pred_conditions": len(pred[g]), "n_true_conditions": len(truth[g]),
                     "pred_conditions": sorted(pred[g]), "true_conditions": sorted(truth[g])})

    print(f"\n--- partition of the {len(genes)} two-sided genes on this {len(keys)}-condition panel ---")
    for c in sorted(counts, key=lambda k: -counts[k]):
        print(f"  {c:22} {counts[c]:4}  ({counts[c]/len(genes):5.1%})")

    never = [r for r in rows if r["class"] == "NEVER_FIRES"]
    mis = [r for r in rows if r["class"] == "MIS_CONDITIONED"]
    part = [r for r in rows if r["class"] == "PARTIAL_OVERLAP"]
    print(f"\nNEVER_FIRES     {len(never):4}  the model never calls these essential anywhere -> "
          f"objective / network problem")
    print(f"MIS_CONDITIONED {len(mis):4}  fires, but never where the experiment says -> "
          f"condition-modelling problem")
    print(f"PARTIAL_OVERLAP {len(part):4}  fires in at least one right condition")
    if part:
        tot_p = sum(len(set(r['pred_conditions']) & set(r['true_conditions'])) for r in part)
        tot_t = sum(r["n_true_conditions"] for r in part)
        print(f"                      of their {tot_t} true cells, {tot_p} are caught "
              f"({tot_p/tot_t:.1%})")

    # ---- CONTINUITY. Every prior artifact in this arc is scored against the 131 genes from the
    # 11-condition expression panel. This script's own panel is 25 conditions, which yields a LARGER
    # two-sided set -- a different denominator. Report the restricted view too so the accounting can be
    # carried forward without silently swapping denominators mid-arc.
    restricted = None
    ratios_art = Path("wiki/fba_eflux_bridge_2026-08-17_ratios.json")
    if ratios_art.exists():
        d = json.loads(ratios_art.read_text(encoding="utf-8"))
        prior: set[str] = set()
        for cells in d["arms"]["baseline"].values():
            prior |= set(cells)
        sub = [r for r in rows if r["gene"] in prior]
        rc: dict[str, int] = {}
        for r in sub:
            rc[r["class"]] = rc.get(r["class"], 0) + 1
        restricted = {"gene_set": "the 131 from the 11-condition expression panel",
                      "n_genes": len(sub), "class_counts": rc}
        print(f"\n--- restricted to the prior {len(sub)}-gene set (continuity) ---")
        for c in sorted(rc, key=lambda k: -rc[k]):
            print(f"  {c:22} {rc[c]:4}  ({rc[c]/len(sub):5.1%})")

    out = {
        "record": "fba-missed-gene-partition-v1",
        "axis": a.axis,
        "date": date.today().isoformat(),
        "model": model.id,
        "n_conditions": len(keys), "conditions": keys,
        "n_genes": len(genes),
        "class_counts": counts,
        "restricted_to_prior_gene_set": restricted,
        "denominator_warning": (
            "This panel is 25 conditions and yields a LARGER two-sided gene set than the 11-condition "
            "expression panel used by every prior artifact in this arc. Do not mix the two denominators; "
            "`restricted_to_prior_gene_set` carries the continuity view."),
        "genes": rows,
        "interpretation": {
            "PROVABLY_UNCALLABLE": "no medium/objective/constraint layer can move them (prior work)",
            "NEVER_FIRES": ("structurally capable but the model never predicts essential in ANY "
                            "labelled condition -- the objective never demands the product"),
            "MIS_CONDITIONED": ("the model DOES fire, just never in a condition the experiment calls "
                                "essential -- a condition-modelling / calibration problem"),
            "PARTIAL_OVERLAP": "fires in at least one correct condition; misses the rest",
        },
        "caveats": [
            "The 'ever fires' test spans ALL labelled carbon conditions, not just the 11 with "
            "expression -- deliberately generous, so NEVER_FIRES is a conservative claim.",
            "Condition modelling here is sole-carbon-source only. A gene whose true dependence is on "
            "aeration, pH or a nitrogen source cannot be caught by this panel and will read as "
            "NEVER_FIRES.",
            "Model reach, not biology: this partitions the MODEL's failure modes, not the cell's.",
        ],
        "needs_fitness_browser_db": True,
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
