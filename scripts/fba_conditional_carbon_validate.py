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

FRAC = 0.01


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
    nonoptimal: dict[str, int] = {}
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            wt_by_cond[cond] = round(wt, 4)
            d: dict[str, bool] = {}
            if wt > 1e-9:
                res = single_gene_deletion(model, gene_list=[model.genes.get_by_id(g) for g in genes])
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    # cobrapy returns a `status` column; a non-optimal solve is otherwise
                    # indistinguishable from a real growth value, so count them as an audit field.
                    st = row.get("status") if hasattr(row, "get") else None
                    if st is not None and st != "optimal":
                        nonoptimal[cond] = nonoptimal.get(cond, 0) + 1
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
    if nonoptimal:
        print(f"   NON-OPTIMAL solver statuses: {nonoptimal}")
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
        "n_solver_nonoptimal_by_condition": nonoptimal,
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
