"""Does the ~5% conditional-essentiality finding hold at 25 carbon sources instead of 4 media?

    uv run python scripts/fba_conditional_carbon_validate.py

`wiki/fba_conditional_essentiality_2026-08-12.md` measured, on the Orth 2011 4-media screen, that iML1515
reproduces the medium-dependent essentiality switch for **3/67 genes (4.5%)** and predicts a CONSTANT
pattern for **94%** of the genes whose essentiality actually depends on the medium — barely above a
constant-predictor null.

That was 268 gene x condition cells. This re-asks the same question against the Fitness Browser RB-TnSeq
compendium: **25 carbon sources, ~1,300 genes, ~33,000 cells** — roughly 125x the resolution. Either the
finding holds at scale or it does not, and both answers are worth having.

Same metric, same nulls, same honesty rails as the 4-media cell (`switch_accuracy` / `constant_baselines`
/ `pattern_distribution` are reused verbatim, generalised only to accept a condition set).

**Two differences from the 4-media substrate, both recorded in the artifact:**
  1. **No reproduction gate.** Orth ships the paper's own iJO1366 FBA columns; the Fitness Browser does
     not. There is no published prediction to check the pipeline against before trusting a new number.
  2. **Aerobicity is not varied.** These are all aerobic carbon-source assays, so the anaerobic axis that
     the 4-media set carries is absent. This tests carbon-source specificity, not oxygen response.

Exit 0 always: this is an experiment, not a gate.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    conditionally_essential_genes,
    confusion_from_calls,
    constant_baselines,
    mcc,
    pattern_distribution,
    switch_accuracy,
)
from dna_decode.fba.fitness_browser import (  # noqa: E402
    ESSENTIAL_FITNESS,
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402
from dna_decode.fba.solver_audit import audit_deletion_frame, merge_audits  # noqa: E402

FRAC = 0.01


def commit_strata(subset, calls: dict[str, dict[str, bool]], keys: tuple[str, ...],
                  suspect: set[tuple[str, str]] | None = None) -> dict:
    """Partition the scored genes by WHAT THE MODEL PREDICTED, not by whether it was right.

    "FBA commits rarely but is accurate when it commits" is arithmetically true and still misleading on
    its own, because 25 of the 33 commitments predict essentiality in exactly ONE of 25 conditions --
    the concentrated end, where the true labels are sparse the same way. Splitting
    constant / 1-of-N / 2+-of-N is what makes that visible instead of buried in a caveat.

    Three strata, exactly partitioning the gene set:
      predicted_constant  -- all-dispensable or all-essential; CANNOT exact-match a two-sided gene
      predicted_1_of_n    -- commits to essentiality in exactly one condition
      predicted_2plus     -- commits to a genuinely multi-condition pattern
    """
    suspect = suspect or set()
    strata = {"predicted_constant": [], "predicted_1_of_n": [], "predicted_2plus": []}
    for r in subset:
        n_ess = sum(1 for c in keys if calls.get(c, {}).get(r.gene_id, False))
        if n_ess == 0 or n_ess == len(keys):
            strata["predicted_constant"].append(r)
        elif n_ess == 1:
            strata["predicted_1_of_n"].append(r)
        else:
            strata["predicted_2plus"].append(r)

    out = {}
    for name, recs in strata.items():
        exact, right, total, touched = 0, 0, 0, 0
        for r in recs:
            pred_set = {c for c in keys if calls.get(c, {}).get(r.gene_id, False)}
            true_set = {c for c in keys if r.experimental.get(c, False)}
            if pred_set == true_set:
                exact += 1
            if any((r.gene_id, c) in suspect for c in keys):
                touched += 1
            for c in keys:
                total += 1
                if calls.get(c, {}).get(r.gene_id, False) == r.experimental.get(c, False):
                    right += 1
        out[name] = {
            "n_genes": len(recs),
            "n_exact_set_match": exact,
            "per_cell_agreement": round(right / total, 4) if total else None,
            "n_genes_touching_a_nonoptimal_cell": touched,
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None, help="path to feba.db (default D:/dna_decode_cache/...)")
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
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
    print(f"{model.id}: {len(conds)} mappable carbon sources")

    model_genes = {g.id for g in model.genes}
    records = load_records(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    print(f"genes with a complete row across all {len(keys)} conditions: {len(records)}")
    print(f"CONDITIONALLY ESSENTIAL (fit<{a.threshold} in >=1, not all): {len(subset)}"
          f"   [the 4-media Orth set had 68]")
    if not subset:
        print("no two-sided genes -- nothing to score", file=sys.stderr)
        return 2

    genes = [r.gene_id for r in subset]
    all_ex = tuple(conds.values())
    print(f"\nrunning FBA deletions: {len(genes)} genes x {len(keys)} conditions "
          f"= {len(genes) * len(keys):,} knockouts ...", flush=True)

    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415
    calls: dict[str, dict[str, bool]] = {}
    wt_by_cond: dict[str, float] = {}
    audits = {}
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            wt_by_cond[cond] = round(wt, 4)
            d: dict[str, bool] = {}
            if wt > 1e-9:
                res = single_gene_deletion(model, gene_list=[model.genes.get_by_id(g) for g in genes])
                # Per-CELL, not per-condition counts: the earlier count-only audit found 39 non-optimal
                # solves but could never answer whether they concentrate in the genes where the model
                # actually commits -- one such cell can create or destroy a 1-of-25 exact-set match.
                audits[cond] = audit_deletion_frame(res, cond)
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    d[gid] = (g != g) or (g < FRAC * wt)
            calls[cond] = d
        print(f"   [{n:2d}/{len(keys)}] {cond[:38]:40s} wt {wt:.4f} | {sum(d.values()):4d} called essential",
              flush=True)

    sw = switch_accuracy(subset, calls, conditions=keys)
    nulls = constant_baselines(subset, conditions=keys)
    best_null = max(g["per_condition_agreement"] for g in nulls.values())
    pat = pattern_distribution(subset, calls, conditions=keys)
    true_pat = pattern_distribution(subset, None, conditions=keys)
    mccs = []
    for c in keys:
        cm = confusion_from_calls({r.gene_id: r.experimental[c] for r in subset}, calls[c])
        mccs.append(mcc(cm))

    print(f"\n=== CONDITIONAL SWITCH over {len(keys)} carbon sources ===")
    print(f"   exact-set match      {sw['exact_set_match']}/{sw['n_conditionally_essential']} "
          f"({sw['exact_set_match_rate']})")
    print(f"   per-cell agreement   {sw['per_condition_agreement']}")
    print(f"   best constant null   {best_null}  -> lift "
          f"{sw['per_condition_agreement'] - best_null:+.4f}")
    print(f"   mean per-condition MCC {sum(mccs) / len(mccs):.4f}")
    print(f"   model predicts a CONSTANT pattern for {pat['n_constant_pattern']}/{pat['n_genes']} "
          f"({pat['constant_pattern_fraction']}) across {pat['n_distinct_patterns']} shapes")
    committed = pat["n_genes"] - pat["n_constant_pattern"]
    print(f"   -> commits to a VARYING pattern for {committed}/{pat['n_genes']} genes; of those, "
          f"{sw['exact_set_match']} are exactly right "
          f"({100 * sw['exact_set_match'] / committed if committed else 0:.0f}% where it commits)")
    print(f"   true patterns: {true_pat['n_distinct_patterns']} shapes, "
          f"{true_pat['n_constant_pattern']} constant (0 by definition)")
    audit = merge_audits(audits) if audits else None
    suspect = {tuple(x) for x in ((audit or {}).get("suspect_cells") or [])}
    if audit and audit["n_suspect_total"]:
        print(f"   SOLVER AUDIT: {audit['n_suspect_total']} suspect cells of {audit['n_rows_total']} "
              f"across {audit['n_conditions_with_nonoptimal']} conditions (f={audit['suspect_fraction']})")

    strata = commit_strata(subset, calls, keys, suspect)
    print("\n   COMMIT STRATA (what the model predicted, not whether it was right):")
    for name, s in strata.items():
        print(f"      {name:20s} n={s['n_genes']:4d} | exact {s['n_exact_set_match']:3d} | "
              f"per-cell {s['per_cell_agreement']} | touching a suspect cell "
              f"{s['n_genes_touching_a_nonoptimal_cell']}")

    # M4 cross-check: a global 0.7% suspect rate does not clear a CONCENTRATED subset.
    committed_genes = [r for r in subset
                       if 0 < sum(1 for c in keys if calls.get(c, {}).get(r.gene_id, False)) < len(keys)]
    exact_genes = [r for r in subset
                   if {c for c in keys if calls.get(c, {}).get(r.gene_id, False)}
                   == {c for c in keys if r.experimental.get(c, False)}]
    def pred_ess(r):
        return {c for c in keys if calls.get(c, {}).get(r.gene_id, False)}

    # "Touches a suspect cell" is the weak form -- the suspect cell might sit in a condition the model
    # called dispensable, in which case it did not manufacture the prediction. The DECISIVE question is
    # whether the ESSENTIAL calls themselves are the failed solves.
    enrichment = {
        "n_committed_genes": len(committed_genes),
        "n_committed_touching_a_suspect_cell": sum(
            1 for r in committed_genes if any((r.gene_id, c) in suspect for c in keys)),
        "n_committed_whose_essential_calls_are_ALL_suspect": sum(
            1 for r in committed_genes if pred_ess(r) and all((r.gene_id, c) in suspect
                                                              for c in pred_ess(r))),
        "n_exact_set_matches": len(exact_genes),
        "n_exact_matches_touching_a_suspect_cell": sum(
            1 for r in exact_genes if any((r.gene_id, c) in suspect for c in keys)),
        "n_exact_matches_whose_essential_calls_are_ALL_suspect": sum(
            1 for r in exact_genes if pred_ess(r) and all((r.gene_id, c) in suspect
                                                          for c in pred_ess(r))),
        "n_exact_matches_that_are_constant_predictions": sum(
            1 for r in exact_genes if not pred_ess(r) or pred_ess(r) == set(keys)),
    }
    print(f"   non-optimal enrichment: {enrichment['n_committed_touching_a_suspect_cell']}"
          f"/{enrichment['n_committed_genes']} committed genes and "
          f"{enrichment['n_exact_matches_touching_a_suspect_cell']}"
          f"/{enrichment['n_exact_set_matches']} exact matches TOUCH a suspect cell")
    print(f"   DECISIVE: {enrichment['n_committed_whose_essential_calls_are_ALL_suspect']}"
          f"/{enrichment['n_committed_genes']} committed genes and "
          f"{enrichment['n_exact_matches_whose_essential_calls_are_ALL_suspect']}"
          f"/{enrichment['n_exact_set_matches']} exact matches have essential calls that are ALL "
          f"failed solves")
    print("   (4-media baseline: exact-set 3/67 = 0.0448, per-cell 0.5709, lift +0.0121, constant 94.0%)")

    result = {
        "record": "fba-conditional-carbon-v1",
        "date": a.date, "model": model.id,
        "substrate": "Fitness Browser RB-TnSeq (figshare 25236931, CC BY 4.0), orgId=Keio",
        "essentiality_threshold_fitness": a.threshold,
        "n_conditions": len(keys), "conditions": {k: conds[k] for k in keys},
        "n_genes_complete_rows": len(records),
        "n_conditionally_essential": len(subset),
        "wildtype_growth_per_condition": wt_by_cond,
        "switch": sw, "null_controls": nulls,
        "lift_over_best_constant_null": round(sw["per_condition_agreement"] - best_null, 4),
        "mean_per_condition_mcc": round(sum(mccs) / len(mccs), 4),
        "pattern_distribution_predicted": pat,
        "pattern_distribution_experimental": true_pat,
        "solver_audit": audit,
        "commit_strata": strata,
        "nonoptimal_enrichment": enrichment,
        "four_media_baseline": {"exact_set_match": "3/67", "per_condition_agreement": 0.5709,
                                "lift": 0.0121, "constant_pattern_fraction": 0.9403},
        "caveats": [
            "NO reproduction gate: unlike the Orth substrate there is no published FBA column here, so the "
            "pipeline cannot be checked against a prior published prediction before trusting this number.",
            "All conditions are AEROBIC carbon-source assays -- the anaerobic axis the 4-media set carries "
            "is absent. This measures carbon-source specificity, not oxygen response.",
            "RB-TnSeq fitness is a POOLED competition readout, not a pure growth/no-growth call; fit<-2 is "
            "inherited from the shipped Keio validation (Bernstein 2023) for comparability.",
            "Replicate experiments per carbon source are averaged (62 experiments over 25 sources).",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_conditional_carbon_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
