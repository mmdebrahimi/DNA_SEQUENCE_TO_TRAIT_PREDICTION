"""Nitrogen conditional essentiality — the first non-carbon substrate axis.

Bar frozen in `wiki/fba_nitrogen_prereg_2026-08-17.md` BEFORE any solve, with three pre-registered
carbon->nitrogen predictions (P1 bimodality, P2 flat-deficit, P3 beats-constant-null) and a MANDATORY
determinism requirement: the full panel runs TWICE at processes=1 and the two runs must agree exactly,
or no verdict is reported.

Usage:
    uv run python scripts/fba_conditional_nitrogen.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    conditionally_essential_genes,
    confusion_from_calls,
    mcc,
)
from dna_decode.fba.fitness_browser import ESSENTIAL_FITNESS, open_db  # noqa: E402
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402
from dna_decode.fba.nitrogen import (  # noqa: E402
    NITROGEN_UNMAPPABLE,
    apply_nitrogen_condition,
    determinism_verdict,
    load_nitrogen_records,
    nitrogen_conditions,
    redact_unverified,
)

FRAC = 0.01


def score_panel(ratios, subset, keys, genes):
    """(per-cell agreement, per-condition confusion) for ONE pass.

    Split out so the determinism gate can recompute the headline metric on the SECOND pass and require
    the two to be identical -- the check that would have caught today's retracted result, where the
    conclusion moved between runs while every individual solve looked fine.
    """
    right = total = 0
    per_cond = {}
    for c in keys:
        e = {r.gene_id: r.experimental[c] for r in subset}
        p = {g: (ratios[c].get(g, 0.0) <= FRAC) for g in genes}
        cm = confusion_from_calls(e, p)
        right += cm["tp"] + cm["tn"]
        total += cm["n"]
        per_cond[c] = {"agreement": round((cm["tp"] + cm["tn"]) / cm["n"], 4) if cm["n"] else None,
                       "mcc": round(mcc(cm), 4), **cm}
    return (round(right / total, 4) if total else None), per_cond


def run_panel(model, conds, keys, genes, single_gene_deletion):
    """One full deterministic pass. Returns {cond: {gene: ratio}} + wildtype growth per condition."""
    all_n = tuple(conds.values())
    ratios, wt_by = {}, {}
    for n, cond in enumerate(keys, 1):
        with model:
            apply_nitrogen_condition(model, conds[cond], all_nitrogen=all_n)
            wt = wildtype_growth(model)
            wt_by[cond] = round(float(wt), 6)
            rt = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    rt[gid] = 0.0 if g != g else max(0.0, float(g) / wt)
            else:
                rt = {g: 0.0 for g in genes}
            ratios[cond] = rt
        print(f"  [{n}/{len(keys)}] {cond:30} wt={wt_by[cond]:.4f}", flush=True)
    return ratios, wt_by


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=f"wiki/fba_conditional_nitrogen_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    model = load_model()
    conn = open_db()
    conds = nitrogen_conditions(conn, model)
    keys = sorted(conds)
    print(f"nitrogen conditions mapped: {len(keys)} (excluded: {len(NITROGEN_UNMAPPABLE)})")
    for k in keys:
        print(f"   {k:30} {conds[k]}")

    model_genes = {g.id for g in model.genes}
    records = load_nitrogen_records(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    print(f"\ngenes with complete rows: {len(records)} | conditionally essential (two-sided): {len(subset)}")
    if not subset:
        print("no two-sided genes on this panel -- nothing to score", file=sys.stderr)
        return 2
    genes = [r.gene_id for r in subset]

    # --- DETERMINISM GATE: two full passes must agree exactly (pre-registered, checked FIRST) ---
    print(f"\n=== PASS 1 ({len(genes)} genes x {len(keys)} conditions) ===", flush=True)
    r1, wt1 = run_panel(model, conds, keys, genes, single_gene_deletion)
    print(f"\n=== PASS 2 (determinism check) ===", flush=True)
    r2, wt2 = run_panel(model, conds, keys, genes, single_gene_deletion)
    ratios = r1
    calls = {c: {g: (ratios[c].get(g, 0.0) <= FRAC) for g in genes} for c in keys}

    # --- metrics (computed for BOTH passes; the gate requires them to be identical) ---
    per_cell, per_cond = score_panel(r1, subset, keys, genes)
    per_cell_2, _ = score_panel(r2, subset, keys, genes)

    gate = determinism_verdict(r1, r2, keys, genes, frac=FRAC,
                               metric_a=per_cell, metric_b=per_cell_2)
    deterministic = gate["deterministic_at_claim_level"] and wt1 == wt2
    print(f"\nDETERMINISM GATE: {'PASS' if deterministic else 'FAIL'}")
    print(f"  call flips across the FRAC={FRAC} line : {gate['n_call_flips']}")
    print(f"  largest numerical drift                : {gate['max_abs_delta']:.3e}")
    print(f"  nearest cell to the threshold          : {gate['min_margin_to_threshold']:.3e}")
    print(f"  safety factor (margin / drift)         : {gate['safety_factor']:.3g} "
          f"(bar {gate['min_safety_factor']:g})")
    print(f"  headline metric both passes            : {gate['headline_metric']}")

    # best-constant null (P3): always-essential / always-dispensable, whichever scores higher
    n_ess = sum(1 for c in keys for r in subset if r.experimental[c])
    n_cells = len(keys) * len(subset)
    null_all_e = n_ess / n_cells
    null_all_d = 1 - null_all_e
    best_const = round(max(null_all_e, null_all_d), 4)

    # P1 bimodality
    flat_vals = [ratios[c].get(g, 0.0) for c in keys for g in genes]
    band = sum(1 for v in flat_vals if 0.001 <= v < 0.05)
    p1_frac = band / len(flat_vals)

    # P2 flat-deficit: of MISSED essential cells, how many had a deletion that changed nothing
    missed = [(c, g) for c in keys for r in subset if (g := r.gene_id)
              and r.experimental[c] and not calls[c][g]]
    flat_missed = sum(1 for c, g in missed if ratios[c].get(g, 0.0) >= 0.999)
    p2_frac = (flat_missed / len(missed)) if missed else None

    verdicts = {
        "P1_bimodality": {"frac_in_band_0.001_0.05": round(p1_frac, 5), "bar": "< 0.01",
                          "result": "REPLICATES" if p1_frac < 0.01 else "FALSIFIED"},
        "P2_flat_deficit": {"frac_missed_essential_that_are_flat": None if p2_frac is None else round(p2_frac, 4),
                            "n_missed": len(missed), "bar": ">= 0.50",
                            "result": None if p2_frac is None else ("REPLICATES" if p2_frac >= 0.50 else "FALSIFIED")},
        "P3_beats_constant_null": {"per_cell": per_cell, "best_constant_null": best_const,
                                   "lift": None if per_cell is None else round(per_cell - best_const, 4),
                                   "bar": "> 0",
                                   "result": None if per_cell is None else
                                   ("REPLICATES" if per_cell > best_const else "FALSIFIED")},
    }

    out = {
        "record": "fba-conditional-nitrogen-v1",
        "date": date.today().isoformat(),
        "prereg": "wiki/fba_nitrogen_prereg_2026-08-17.md",
        "model": "iML1515",
        "labels": f"Fitness Browser RB-TnSeq orgId=Keio, expGroup='nitrogen source', fit<{a.threshold}",
        "n_conditions": len(keys),
        "conditions": conds,
        "excluded": NITROGEN_UNMAPPABLE,
        "n_genes_complete_rows": len(records),
        "n_conditionally_essential": len(subset),
        "deterministic": deterministic,
        "determinism_gate": {**gate, "wildtype_identical": wt1 == wt2},
        "wildtype_growth": wt1,
        "per_cell_agreement": per_cell,
        "best_constant_null": best_const,
        "per_condition": per_cond,
        "predictions": verdicts,
        "verdict": ("NON_DETERMINISTIC_NO_VERDICT" if not deterministic else
                    "ALL_THREE_REPLICATE" if all(
                        v.get("result") == "REPLICATES" for v in verdicts.values())
                    else "MIXED"),
        "caveats": [
            "Several N sources (alanine/serine/aspartate/glutamine) also supply CARBON -- true of the "
            "real assay too (glucose minimal medium + test compound as sole N source), so these are NOT "
            "nitrogen-only perturbations.",
            "13 of 16 assay sources; 2 dipeptides + casamino acids are unmappable to iML1515 exchanges.",
            "2 replicates per source, averaged -- thinner than carbon's ~2.5.",
        ],
    }
    # The numbers are REMOVED, not flagged, when the panel is not reproducible. The first run of this
    # script wrote `deterministic: false` beside a full set of quotable verdicts; a flag next to the
    # numbers it invalidates is not a control.
    out = redact_unverified(out, deterministic)

    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    if not deterministic:
        print("\nDETERMINISM GATE FAILED -- every solve-derived number has been WITHHELD from the "
              "artifact per the pre-registration. Nothing here may be quoted.")
    else:
        print(f"\nper-cell {per_cell} | best-constant null {best_const}")
        for k, v in verdicts.items():
            print(f"  {k:24} {v.get('result')}")
    print(f"VERDICT: {out['verdict']}\nwrote {a.out}.json")
    return 0 if deterministic else 1


if __name__ == "__main__":
    raise SystemExit(main())
