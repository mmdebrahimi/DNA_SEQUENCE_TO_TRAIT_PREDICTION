"""WHY does iML1515 refuse to switch on 184 of 217 conditionally-essential genes?

    uv run python scripts/fba_constant_gene_diagnostic.py

This is where the conditional deficit actually lives. The model commits to a varying pattern for only
33/217 genes, and those commitments are now explained: they are sole-route catabolism detected via
ATPM-maintenance infeasibility (`wiki/fba_infeasibility_finding_2026-08-13.md`). The other **184 are
constant** — 145 predicted dispensable-everywhere, 39 predicted essential-everywhere — and nothing has
ever asked why.

The question matters because it separates two very different diagnoses, which imply opposite fixes:

  **MODEL problem** -- the knockout growth ratio is FLAT at ~1.0 in the conditions where the gene is
  truly essential. The model has an alternative route the real cell does not use. No readout change can
  recover this; it needs regulation (the pFBA result) or gene-content curation.

  **READOUT problem** -- the ratio is materially depressed but sits above the 1% cutoff. The model DOES
  see a growth defect and the binary threshold throws it away. A graded readout recovers it for free.

The 4-media run reported "64% of these genes have a perfectly flat ratio", but that was 67 genes on 4
conditions with no per-cell breakdown. This measures it properly on 25 carbon sources, split by
prediction stratum, and adds the mechanistic cross-check the earlier number lacked: **GPR structure**.
A gene whose reactions are all covered by an isozyme (`OR` in the gene-reaction rule) is redundant BY
CONSTRUCTION, which is a checkable cause rather than a guess.

Exit 0 always: this is a diagnostic, not a gate.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from datetime import date
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
from dna_decode.fba.solver_audit import audit_deletion_frame, merge_audits  # noqa: E402

FRAC = 0.01
FLAT_EPS = 1e-6          # a ratio this close to 1.0 is "the deletion did nothing at all"
NEAR_MISS_HI = 0.90      # below this, the model saw a MATERIAL growth defect


def classify_miss(ratio: float | None) -> str:
    """A gene the model called dispensable in a condition where truth says ESSENTIAL. Why? PURE.

    `flat`          -- ratio indistinguishable from 1.0: the deletion changed nothing. The model has a
                       redundant route. A MODEL problem; no readout change recovers it.
    `near_threshold`-- ratio in (FRAC, 0.10]: the model nearly called it. A readout problem, and the
                       cheapest possible one to fix.
    `material`      -- ratio in (0.10, 0.90]: a real, sizeable growth defect the binary cutoff discards.
                       A READOUT problem -- a graded metric recovers this for free.
    `slight`        -- ratio in (0.90, 1.0): a measurable but small defect. Mostly a model problem.
    `no_value`      -- no growth value (should not occur for a dispensable call).
    """
    if ratio is None or ratio != ratio:
        return "no_value"
    if ratio >= 1.0 - FLAT_EPS:
        return "flat"
    if ratio <= 0.10:
        return "near_threshold"
    if ratio <= NEAR_MISS_HI:
        return "material"
    return "slight"


def gene_is_isozyme_redundant(model, gene_id: str) -> bool | None:
    """Is EVERY reaction this gene participates in also covered by another gene? PURE-ish (reads model).

    A reaction whose gene-reaction rule contains `or` at the top level can run without this gene, so a
    deletion cannot block it. If that holds for all of the gene's reactions, the model CANNOT call it
    essential no matter the medium -- redundancy by construction, not by parameter choice.
    """
    try:
        gene = model.genes.get_by_id(gene_id)
    except KeyError:
        return None
    rxns = list(gene.reactions)
    if not rxns:
        return None
    for rxn in rxns:
        rule = (rxn.gene_reaction_rule or "").lower()
        others = {g.id for g in rxn.genes} - {gene_id}
        if not others or " or " not in rule:
            return False            # at least one reaction depends on this gene alone
    return True


def diagnose(subset, calls, ratios, keys, model=None) -> dict:
    """Split the constant-predicted genes and explain each stratum. PURE apart from optional GPR reads."""
    strata = {"predicted_all_dispensable": [], "predicted_all_essential": [], "committed": []}
    for r in subset:
        n_ess = sum(1 for c in keys if calls.get(c, {}).get(r.gene_id, False))
        if n_ess == 0:
            strata["predicted_all_dispensable"].append(r)
        elif n_ess == len(keys):
            strata["predicted_all_essential"].append(r)
        else:
            strata["committed"].append(r)

    # --- the 145: missed conditional essentiality. WHY did the model keep growing? ---
    miss_kinds: dict[str, int] = {}
    miss_ratios: list[float] = []
    per_gene_flat: list[str] = []
    for r in strata["predicted_all_dispensable"]:
        true_ess = [c for c in keys if r.experimental.get(c, False)]
        kinds = []
        for c in true_ess:
            v = ratios.get(c, {}).get(r.gene_id)
            k = classify_miss(v)
            miss_kinds[k] = miss_kinds.get(k, 0) + 1
            kinds.append(k)
            if v is not None and v == v:
                miss_ratios.append(v)
        if kinds and all(k == "flat" for k in kinds):
            per_gene_flat.append(r.gene_id)

    total_miss = sum(miss_kinds.values())
    recoverable = miss_kinds.get("near_threshold", 0) + miss_kinds.get("material", 0)

    out = {
        "n_scored": len(subset),
        "n_predicted_all_dispensable": len(strata["predicted_all_dispensable"]),
        "n_predicted_all_essential": len(strata["predicted_all_essential"]),
        "n_committed": len(strata["committed"]),
        "missed_essential_cells": total_miss,
        "missed_cells_by_cause": miss_kinds,
        "n_genes_flat_in_EVERY_true_essential_condition": len(per_gene_flat),
        "readout_recoverable_cells": recoverable,
        "readout_recoverable_fraction": (round(recoverable / total_miss, 4) if total_miss else None),
        "missed_ratio_median": (round(statistics.median(miss_ratios), 4) if miss_ratios else None),
        "example_flat_genes": sorted(per_gene_flat)[:15],
    }

    if model is not None:
        flat_set = set(per_gene_flat)
        red = [g for g in flat_set if gene_is_isozyme_redundant(model, g)]
        out["n_flat_genes_isozyme_redundant"] = len(red)
        out["isozyme_share_of_flat_genes"] = (round(len(red) / len(flat_set), 4) if flat_set else None)
        out["example_isozyme_redundant"] = sorted(red)[:15]

    # --- the 39: over-called essential everywhere. Where truth says dispensable, what did FBA see? ---
    over = {"infeasible": 0, "sub_threshold": 0}
    for r in strata["predicted_all_essential"]:
        for c in keys:
            if r.experimental.get(c, False):
                continue                                    # truth agrees here
            v = ratios.get(c, {}).get(r.gene_id)
            over["infeasible" if (v is None or v != v) else "sub_threshold"] += 1
    out["overcalled_dispensable_cells_by_cause"] = over
    return out


def verdict_for(d: dict) -> str:
    """Which diagnosis dominates the deficit? PURE."""
    frac = d.get("readout_recoverable_fraction")
    if frac is None:
        return "NO_MISSES_TO_EXPLAIN"
    if frac >= 0.5:
        return "DEFICIT_IS_MOSTLY_A_READOUT_PROBLEM"
    if frac >= 0.2:
        return "DEFICIT_IS_MIXED_MODEL_AND_READOUT"
    return "DEFICIT_IS_A_MODEL_PROBLEM_NOT_A_READOUT_ONE"


def main(argv: list[str] | None = None) -> int:
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None)
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    try:
        conn = open_db(a.db)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 2

    model = load_model()
    conds = carbon_conditions(conn, model)
    keys = tuple(sorted(conds))
    subset = conditionally_essential_genes(
        load_records(conn, conds, gene_filter={g.id for g in model.genes}))
    genes = [r.gene_id for r in subset]
    print(f"{model.id}: {len(keys)} carbon sources | {len(genes)} conditionally-essential genes")

    calls: dict[str, dict[str, bool]] = {}
    ratios: dict[str, dict[str, float]] = {}
    audits = {}
    all_ex = tuple(conds.values())
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            d, rat = {}, {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes])
                audits[cond] = audit_deletion_frame(res, cond)
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    gv = row["growth"]
                    rat[gid] = None if gv != gv else float(gv) / wt
                    d[gid] = (gv != gv) or (gv < FRAC * wt)
            calls[cond] = d
            ratios[cond] = rat
        print(f"   [{n:2d}/{len(keys)}] {cond[:38]:40s} wt {wt:.4f}", flush=True)

    diag = diagnose(subset, calls, ratios, keys, model=model)
    v = verdict_for(diag)

    print(f"\n=== WHY THE MODEL DOES NOT SWITCH ({diag['n_scored']} genes) ===")
    print(f"   predicted all-dispensable : {diag['n_predicted_all_dispensable']}")
    print(f"   predicted all-essential   : {diag['n_predicted_all_essential']}")
    print(f"   commits to a pattern      : {diag['n_committed']}")
    print(f"\n   MISSED essential cells    : {diag['missed_essential_cells']}")
    for k, n in sorted(diag["missed_cells_by_cause"].items(), key=lambda x: -x[1]):
        print(f"      {k:16s} {n:5d}")
    print(f"   median missed ratio       : {diag['missed_ratio_median']}")
    print(f"   READOUT-recoverable       : {diag['readout_recoverable_cells']} "
          f"({diag['readout_recoverable_fraction']}) -- a graded metric would catch these")
    print(f"   genes FLAT in every true-essential condition: "
          f"{diag['n_genes_flat_in_EVERY_true_essential_condition']}")
    if "n_flat_genes_isozyme_redundant" in diag:
        print(f"      of which isozyme-redundant by GPR: {diag['n_flat_genes_isozyme_redundant']} "
              f"({diag['isozyme_share_of_flat_genes']})")
        print(f"      examples: {', '.join(diag['example_isozyme_redundant'][:8])}")
    print(f"\n   over-called (essential where truth says dispensable): "
          f"{diag['overcalled_dispensable_cells_by_cause']}")
    print(f"\nVERDICT: {v}")

    result = {
        "record": "fba-constant-gene-diagnostic-v1", "date": a.date, "model": model.id,
        "n_conditions": len(keys), "essentiality_frac": FRAC,
        "diagnosis": diag, "verdict": v,
        "solver_audit": merge_audits(audits) if audits else None,
        "caveats": [
            "`flat` uses an exact-1.0 tolerance of 1e-6; degenerate LP optima can shift a mid-range "
            "ratio between runs, but a ratio AT 1.0 is stable (the deletion changed nothing).",
            "Isozyme redundancy is read from the GPR `or` structure, which is a MODEL property -- it "
            "says the model cannot call the gene essential, not that the real cell is redundant.",
            "A `material` ratio being READOUT-recoverable is an upper bound: recovering it needs a "
            "threshold that does not simultaneously destroy precision elsewhere, which is untested here.",
            "All 25 conditions are aerobic carbon sources; no oxygen axis.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_constant_gene_diagnostic_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
