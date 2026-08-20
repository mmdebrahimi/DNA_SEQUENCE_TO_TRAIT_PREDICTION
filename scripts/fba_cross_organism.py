"""Cross-organism conditional essentiality — the first non-E.-coli test.

Bar frozen in `wiki/fba_cross_organism_prereg_2026-08-20.md` BEFORE any solve. Same three predictions as
carbon and nitrogen (P1 bimodality, P2 flat-deficit, P3 beats-constant-null), now on a DIFFERENT organism
and a DIFFERENT genome-scale model.

Reuses, unchanged: `apply_carbon_condition` (already organism-agnostic), and the nitrogen axis's
claim-level `determinism_verdict` + `redact_unverified`.

Usage:
    uv run python scripts/fba_cross_organism.py
    uv run python scripts/fba_cross_organism.py --org-id Keio --organism escherichia_coli   # control
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.carbon_growth import build_exchange_name_index, match_carbon_exchange  # noqa: E402
from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    GeneRecord,
    conditionally_essential_genes,
    confusion_from_calls,
    mcc,
)
from dna_decode.fba.fitness_browser import (  # noqa: E402
    ESSENTIAL_FITNESS,
    apply_carbon_condition,
    open_db,
)
from dna_decode.fba.model import load_model, organism_for, wildtype_growth  # noqa: E402
from dna_decode.fba.nitrogen import determinism_verdict, redact_unverified  # noqa: E402

FRAC = 0.01


def org_conditions(conn: sqlite3.Connection, model, org_id: str, exp_group: str) -> dict[str, str]:
    """{condition -> exchange} for one organism, using the shared carbon name index."""
    idx = build_exchange_name_index(model)
    out = {}
    for (cond,) in conn.execute(
            "SELECT DISTINCT condition_1 FROM Experiment WHERE orgId=? AND expGroup=?",
            (org_id, exp_group)):
        ex = match_carbon_exchange(cond, idx)
        if ex:
            out[cond] = ex
    return out


def load_org_records(conn: sqlite3.Connection, org_id: str, exp_group: str,
                     conditions: dict[str, str], gene_filter: set[str] | None = None,
                     threshold: float = ESSENTIAL_FITNESS) -> list[GeneRecord]:
    """Generalised `load_records`: the SAME contract, for any organism.

    `fitness_browser.load_records` hardcodes `orgId='Keio'` and `expGroup='carbon source'`. Rather than
    edit that load-bearing, heavily-pinned function, this mirrors it for an arbitrary organism.
    Replicates are averaged; a gene is kept only if it has a value in EVERY condition, so a missing
    measurement can never look like a dispensability call.
    """
    keys = tuple(sorted(conditions))
    exp_to_cond = {}
    for name, cond in conn.execute(
            "SELECT expName, condition_1 FROM Experiment WHERE orgId=? AND expGroup=?",
            (org_id, exp_group)):
        if cond in conditions:
            exp_to_cond[name] = cond

    agg: dict[tuple[str, str], list[float]] = {}
    for sysname, expname, fit in conn.execute(
            "SELECT g.sysName, f.expName, f.fit FROM GeneFitness f "
            "JOIN Gene g ON g.orgId=f.orgId AND g.locusId=f.locusId WHERE f.orgId=?", (org_id,)):
        cond = exp_to_cond.get(expname)
        if cond is None or not sysname:
            continue
        if gene_filter is not None and sysname not in gene_filter:
            continue
        agg.setdefault((sysname, cond), []).append(float(fit))

    by_gene: dict[str, dict[str, float]] = {}
    for (gene, cond), vals in agg.items():
        by_gene.setdefault(gene, {})[cond] = sum(vals) / len(vals)

    out: list[GeneRecord] = []
    for gene, per_cond in by_gene.items():
        if len(per_cond) != len(keys):
            continue
        exp = {c: (per_cond[c] < threshold) for c in keys}
        out.append(GeneRecord(gene_id=gene, gene=gene, experimental=exp,
                              paper_fba={c: False for c in keys},
                              conditionally_essential=any(exp.values()) and not all(exp.values())))
    return out


def run_panel(model, conds, keys, genes, single_gene_deletion):
    all_c = tuple(conds.values())
    ratios, wt_by = {}, {}
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_c)
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
        print(f"  [{n}/{len(keys)}] {cond[:38]:40} wt={wt_by[cond]:.4f}", flush=True)
    return ratios, wt_by


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--org-id", default="Putida", help="Fitness Browser orgId")
    ap.add_argument("--organism", default="pseudomonas_putida", help="BiGG registry organism alias")
    ap.add_argument("--exp-group", default="carbon source")
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    out_base = a.out or f"wiki/fba_cross_organism_{a.org_id.lower()}_{date.today().isoformat()}"

    from cobra.flux_analysis import single_gene_deletion

    model = load_model(organism=a.organism)
    conn = open_db()
    print(f"model {model.id} -> {organism_for(model.id)}")
    print(f"   reactions {len(model.reactions)} | genes {len(model.genes)}")

    conds = org_conditions(conn, model, a.org_id, a.exp_group)
    keys = sorted(conds)
    print(f"\n{a.exp_group} conditions mapped: {len(keys)}")
    for k in keys:
        print(f"   {k[:40]:42} {conds[k]}")

    model_genes = {g.id for g in model.genes}
    records = load_org_records(conn, a.org_id, a.exp_group, conds,
                               gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    print(f"\ncomplete-row genes: {len(records)} | conditionally essential: {len(subset)}")
    if not subset:
        print("no two-sided genes -- nothing to score", file=sys.stderr)
        return 2
    genes = [r.gene_id for r in subset]

    print(f"\n=== PASS 1 ({len(genes)} genes x {len(keys)} conditions) ===", flush=True)
    r1, wt1 = run_panel(model, conds, keys, genes, single_gene_deletion)
    print(f"\n=== PASS 2 (determinism check) ===", flush=True)
    r2, wt2 = run_panel(model, conds, keys, genes, single_gene_deletion)

    calls = {c: {g: (r1[c].get(g, 0.0) <= FRAC) for g in genes} for c in keys}
    right = total = 0
    per_cond = {}
    for c in keys:
        e = {r.gene_id: r.experimental[c] for r in subset}
        cm = confusion_from_calls(e, {g: calls[c][g] for g in genes})
        right += cm["tp"] + cm["tn"]
        total += cm["n"]
        per_cond[c] = {"agreement": round((cm["tp"] + cm["tn"]) / cm["n"], 4) if cm["n"] else None,
                       "mcc": round(mcc(cm), 4), **cm}
    per_cell = round(right / total, 4) if total else None

    # second pass metric, for the determinism gate's metric-equality condition
    r2_right = r2_total = 0
    for c in keys:
        e = {r.gene_id: r.experimental[c] for r in subset}
        cm = confusion_from_calls(e, {g: (r2[c].get(g, 0.0) <= FRAC) for g in genes})
        r2_right += cm["tp"] + cm["tn"]
        r2_total += cm["n"]
    per_cell_2 = round(r2_right / r2_total, 4) if r2_total else None

    gate = determinism_verdict(r1, r2, keys, genes, frac=FRAC,
                               metric_a=per_cell, metric_b=per_cell_2)
    deterministic = gate["deterministic_at_claim_level"] and wt1 == wt2
    print(f"\nDETERMINISM GATE: {'PASS' if deterministic else 'FAIL'}")
    print(f"  call flips {gate['n_call_flips']} | max drift {gate['max_abs_delta']:.3g} | "
          f"safety {gate['safety_factor']:.3g} (bar {gate['min_safety_factor']:g}) | "
          f"metric {gate['headline_metric']}")

    n_ess = sum(1 for c in keys for r in subset if r.experimental[c])
    n_cells = len(keys) * len(subset)
    best_const = round(max(n_ess / n_cells, 1 - n_ess / n_cells), 4)

    vals = [r1[c].get(g, 0.0) for c in keys for g in genes]
    p1 = sum(1 for v in vals if 0.001 <= v < 0.05) / len(vals)
    missed = [(c, r.gene_id) for c in keys for r in subset
              if r.experimental[c] and not calls[c][r.gene_id]]
    flat = sum(1 for c, g in missed if r1[c].get(g, 0.0) >= 0.999)
    p2 = (flat / len(missed)) if missed else None

    preds = {
        "P1_bimodality": {"frac_in_band": round(p1, 5), "bar": "< 0.01",
                          "result": "REPLICATES" if p1 < 0.01 else "FALSIFIED"},
        "P2_flat_deficit": {"frac_flat": None if p2 is None else round(p2, 4), "n_missed": len(missed),
                            "bar": ">= 0.50",
                            "result": None if p2 is None else ("REPLICATES" if p2 >= 0.50 else "FALSIFIED")},
        "P3_beats_constant_null": {"per_cell": per_cell, "best_constant_null": best_const,
                                   "lift": None if per_cell is None else round(per_cell - best_const, 4),
                                   "bar": "> 0",
                                   "result": None if per_cell is None else
                                   ("REPLICATES" if per_cell > best_const else "FALSIFIED")},
    }

    payload = {
        "record": "fba-cross-organism-v1",
        "date": date.today().isoformat(),
        "prereg": "wiki/fba_cross_organism_prereg_2026-08-20.md",
        "org_id": a.org_id, "organism": organism_for(model.id), "model": model.id,
        "exp_group": a.exp_group,
        "labels": f"Fitness Browser RB-TnSeq orgId={a.org_id}, fit<{a.threshold}",
        "n_conditions": len(keys), "conditions": conds,
        "n_genes_complete_rows": len(records), "n_conditionally_essential": len(subset),
        "deterministic": deterministic, "determinism": gate,
        "wildtype_growth": wt1,
        "per_cell_agreement": per_cell, "best_constant_null": best_const,
        "per_condition": per_cond, "predictions": preds,
        "verdict": ("NON_DETERMINISTIC_NO_VERDICT" if not deterministic else
                    "ALL_THREE_REPLICATE" if all(v.get("result") == "REPLICATES" for v in preds.values())
                    else "MIXED"),
        "caveats": [
            "n=1 new organism -- a replication, not a survey.",
            "13 of 47 assay carbon sources; the panel is biased toward central metabolism.",
            "iJN1463 is less validated than iML1515; a weak result may reflect model quality.",
        ],
    }
    payload = redact_unverified(payload, deterministic)
    Path(out_base + ".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"\nper-cell {payload.get('per_cell_agreement')} | best-constant null {best_const}")
    for k, v in (payload.get("predictions") or {}).items():
        print(f"  {k:24} {v.get('result')}")
    print(f"VERDICT: {payload['verdict']}\nwrote {out_base}.json")
    return 0 if deterministic else 1


if __name__ == "__main__":
    raise SystemExit(main())
