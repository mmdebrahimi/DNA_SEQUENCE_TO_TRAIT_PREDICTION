"""Validate FBA CONDITIONAL gene essentiality: same gene, four media, does the model get the switch right?

    uv run python scripts/fba_conditional_essentiality_validate.py

Every prior essentiality number in this repo is single-condition, which cannot measure the property a
strain designer depends on -- that a gene is dispensable on one carbon source and required on another.
This scores E. coli iML1515 against the Orth 2011 iJO1366 screen: 1,075 K-12 genes x 4 minimal media,
**68 of them conditionally essential**.

REPRODUCTION GATE FIRST. The supplement ships the paper's OWN iJO1366 FBA calls, so before any new number
is trusted the pipeline scores those against the same experimental labels and checks the result is sane.
The gate has already caught two unit/wiring errors elsewhere in this repo; it is not ceremony.

HONEST scope: iML1515 is the SUCCESSOR of the iJO1366 the paper scored, so a difference is a model
difference, not a reproduction failure -- the gate checks the paper's OWN columns, never ours against
theirs. Labels are the experimental columns only; the FBA columns are never treated as truth.

Exit 0 = scored; 2 = labels missing.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    CONDITIONS,
    apply_condition,
    conditionally_essential_genes,
    confusion_from_calls,
    constant_baselines,
    load_labels,
    mcc,
    pattern_distribution,
    switch_accuracy,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402


def _metrics(cm: dict) -> dict:
    tp, fp, fn, tn = cm["tp"], cm["fp"], cm["fn"], cm["tn"]
    return {
        "n": cm["n"], "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "accuracy": round((tp + tn) / cm["n"], 4) if cm["n"] else None,
        "precision": round(tp / (tp + fp), 4) if (tp + fp) else None,
        "recall": round(tp / (tp + fn), 4) if (tp + fn) else None,
        "mcc": round(mcc(cm), 4),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", default=None, help="gold-standard TSV (defaults to the committed one)")
    ap.add_argument("--organism", default="ecoli")
    ap.add_argument("--frac", type=float, default=0.01, help="essentiality threshold, fraction of wild type")
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    try:
        records = load_labels(a.labels)
    except FileNotFoundError:
        print("conditional-essentiality labels not found; see the module docstring for provenance",
              file=sys.stderr)
        return 2
    cond_ess = conditionally_essential_genes(records)
    print(f"labels: {len(records)} genes x {len(CONDITIONS)} media | "
          f"conditionally essential: {len(cond_ess)}")
    flagged = sum(1 for r in records if r.conditionally_essential)
    print(f"   recomputed {len(cond_ess)} vs supplement's own flag {flagged}"
          f"{'  (AGREE)' if flagged == len(cond_ess) else '  <-- DISAGREE, recomputed value is used'}")

    # ---- reproduction gate: the PAPER's own iJO1366 calls vs the experimental labels ----
    paper_per_cond = {}
    for c in CONDITIONS:
        cm = confusion_from_calls({r.gene_id: r.experimental[c] for r in records},
                                  {r.gene_id: r.paper_fba[c] for r in records})
        paper_per_cond[c] = _metrics(cm)
    paper_switch = switch_accuracy(records, {c: {r.gene_id: r.paper_fba[c] for r in records}
                                             for c in CONDITIONS})
    print("\nREPRODUCTION GATE -- the paper's OWN iJO1366 FBA calls vs the same experimental labels:")
    for c, m in paper_per_cond.items():
        print(f"   {c:20s} acc {m['accuracy']:.4f} | MCC {m['mcc']:.4f} | "
              f"TP {m['tp']:3d} FP {m['fp']:3d} FN {m['fn']:3d} TN {m['tn']:4d}")
    print(f"   conditional switch: exact-set {paper_switch['exact_set_match']}/"
          f"{paper_switch['n_conditionally_essential']} "
          f"({paper_switch['exact_set_match_rate']}) | per-cell {paper_switch['per_condition_agreement']}")

    # NULL CONTROLS -- the switch metric is meaningless without them. Most conditionally-essential genes
    # are essential in only 1-2 of 4 media, so "dispensable everywhere" already scores ~0.56 per-cell.
    nulls = constant_baselines(records)
    print("\nNULL CONTROLS (constant predictors on the same conditionally-essential subset):")
    for name, g in nulls.items():
        print(f"   {name:20s} exact-set {g['exact_set_match']}/{g['n_conditionally_essential']} | "
              f"per-cell {g['per_condition_agreement']}")

    true_pat = pattern_distribution(records)
    paper_pat = pattern_distribution(records, {c: {r.gene_id: r.paper_fba[c] for r in records}
                                               for c in CONDITIONS})
    print(f"\nPATTERN SHAPE on the {true_pat['n_genes']} conditionally-essential genes "
          f"(order {true_pat['conditions_order']}):")
    print(f"   TRUE  : {true_pat['n_distinct_patterns']} distinct shapes; "
          f"constant-pattern genes {true_pat['n_constant_pattern']} (0 by definition)")
    print(f"   iJO1366: constant-pattern genes {paper_pat['n_constant_pattern']}/{paper_pat['n_genes']} "
          f"({paper_pat['constant_pattern_fraction']}) <- it is not switching at all")

    # ---- our model, same four conditions ----
    model = load_model(organism=a.organism)
    model_gene_ids = {g.id for g in model.genes}
    scored = [r for r in records if r.gene_id in model_gene_ids]
    print(f"\n{model.id}: {len(model.genes)} genes | gold-standard genes present in model: "
          f"{len(scored)}/{len(records)}")

    ours: dict[str, dict[str, bool]] = {}
    growth_by_condition = {}
    for c in CONDITIONS:
        with model:
            apply_condition(model, c)
            wt = wildtype_growth(model)
            growth_by_condition[c] = round(wt, 4)
            calls = {}
            if wt > 1e-6:
                from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415
                res = single_gene_deletion(model, gene_list=[model.genes.get_by_id(r.gene_id)
                                                             for r in scored])
                # cobrapy returns an `ids` COLUMN holding a frozenset per row (a RangeIndex, not the ids);
                # older versions put the frozenset in the index. Handle both rather than assume one.
                for idx, row in res.iterrows():
                    ids = row["ids"] if "ids" in res.columns else idx
                    gid = next(iter(ids)) if not isinstance(ids, str) else ids
                    g = row["growth"]
                    calls[gid] = (g != g) or (g < a.frac * wt)     # NaN (infeasible) counts as essential
            ours[c] = calls
        print(f"   {c:20s} wild-type growth {growth_by_condition[c]:.4f} | "
              f"{sum(calls.values())} genes called essential")

    our_per_cond = {}
    for c in CONDITIONS:
        cm = confusion_from_calls({r.gene_id: r.experimental[c] for r in scored}, ours[c])
        our_per_cond[c] = _metrics(cm)
    our_switch = switch_accuracy(scored, ours)

    print(f"\n{model.id} vs the experimental labels:")
    for c, m in our_per_cond.items():
        print(f"   {c:20s} acc {m['accuracy']:.4f} | MCC {m['mcc']:.4f} | "
              f"TP {m['tp']:3d} FP {m['fp']:3d} FN {m['fn']:3d} TN {m['tn']:4d}")
    print(f"   CONDITIONAL SWITCH: exact-set {our_switch['exact_set_match']}/"
          f"{our_switch['n_conditionally_essential']} "
          f"({our_switch['exact_set_match_rate']}) | per-cell {our_switch['per_condition_agreement']}")
    our_pat = pattern_distribution(scored, ours)
    print(f"   pattern shape: constant-pattern genes {our_pat['n_constant_pattern']}/{our_pat['n_genes']} "
          f"({our_pat['constant_pattern_fraction']}) across {our_pat['n_distinct_patterns']} shapes")
    best_null = max(g["per_condition_agreement"] for g in nulls.values())
    print(f"   vs best constant-predictor null {best_null} -> "
          f"lift {our_switch['per_condition_agreement'] - best_null:+.4f}")

    result = {
        "record": "fba-conditional-essentiality-v1",
        "date": a.date,
        "organism": a.organism,
        "model": model.id,
        "labels": {
            "source": "Orth 2011 Mol Syst Biol 7:535 (iJO1366) Supplementary Table 1",
            "n_genes": len(records),
            "conditions": sorted(CONDITIONS),
            "n_conditionally_essential_recomputed": len(cond_ess),
            "n_conditionally_essential_supplement_flag": flagged,
        },
        "reproduction_gate_paper_own_fba": {"per_condition": paper_per_cond, "switch": paper_switch},
        "null_controls": nulls,
        "pattern_distribution": {"experimental": true_pat, "paper_fba": paper_pat},
        "model_scored": {
            "n_genes_scored": len(scored),
            "wildtype_growth_per_condition": growth_by_condition,
            "per_condition": our_per_cond,
            "switch": our_switch,
            "pattern_distribution": our_pat,
        },
        "caveats": [
            "The paper scored iJO1366; this scores its SUCCESSOR iML1515, so a difference vs the "
            "reproduction gate is a MODEL difference, not a reproduction failure.",
            "Labels are the EXPERIMENTAL columns only; the paper's FBA columns are a gate, never truth.",
            "Media are the reconstruction's own M9 mineral background with the carbon source swapped and "
            "the oxygen bound set; they are not independently calibrated to the assay's exact medium.",
            "Lactate is scored as L-lactate (EX_lac__L_e); the assay medium may have been D,L-lactate.",
            "In-distribution vs a published knowledge baseline; not an independent-lab claim.",
            "per_condition_agreement MUST be read against `null_controls`: predicting dispensable "
            "everywhere already scores ~0.56 on this subset, so the models' ~0.57 is ~1 point of lift.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_conditional_essentiality_{a.organism}_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
