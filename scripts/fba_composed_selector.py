"""Composed expression selector: per-gene threshold, single-gene eligibility, JOINT verification.

RUNS `wiki/fba_composed_selector_prereg_2026-08-22.md` EXACTLY (frozen at commit d08a210, BEFORE this
scored anything). Constants trace to numbered sections:

  §2 step 1  per-gene percentile against the gene's OWN distribution over all 1,035 PRECISE-1K samples
             (NOT a genome-wide percentile -- that selector is already falsified)
  §2 step 2  keep only genes whose knockout disables no reaction
  §2 step 3  JOINT verification -- the step both prior attempts lacked. Single-gene safety does NOT
             compose: gate one member of an isozyme pair and the reaction survives; gate BOTH and it dies.
             Drop the highest-expressed member of a collision (tie-break: lexicographic id) until the SET
             disables nothing.
  §2         sensitivity {5, 10, 20}; PRIMARY is 10; best-of reporting is forbidden
  §4 Gate-0  wildtype unchanged (<=1e-6) in EVERY scored condition, asserted at run time. A constructed
             invariant that is not checked is an assumption. Failing it invalidates the whole run.
  §5         processes=1, pipeline twice, identical calls required
  §6         seven pre-committed verdicts

Usage:
    uv run python scripts/fba_composed_selector.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

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
from scripts.fba_eflux_bridge import P1K_DIR, build_condition_expression  # noqa: E402
from scripts.fba_orphan_protection_screen import gpr_disabled_reactions  # noqa: E402

FRAC = 0.01
SCREEN_ARTIFACT = Path("wiki/fba_orphan_protection_2026-08-21.json")
PRIMARY_PCTL = 10.0
SENSITIVITY_PCTLS = (5.0, 10.0, 20.0)
PRIMARY_BAR = 4
GUARDRAIL_REL_INCREASE = 0.20
WT_TOL = 1e-6


def per_gene_reference(model_gene_ids: set[str], pctl: float) -> dict[str, float]:
    """§2 step 1 -- each gene's own pctl over the FULL compendium, in linear TPM."""
    tpm = pd.read_csv(P1K_DIR / "log_tpm_qc.csv", index_col=0)
    lin = 2.0 ** tpm
    return {g: float(np.percentile(lin.loc[g].values, pctl))
            for g in lin.index if g in model_gene_ids}


def disabled_by_set(model, genes: set[str]) -> list[str]:
    """Reactions that become non-functional when the WHOLE set is knocked out."""
    if not genes:
        return []
    with model:
        for g in genes:
            model.genes.get_by_id(g).knock_out()
        seen, out = set(), []
        for g in genes:
            for r in model.genes.get_by_id(g).reactions:
                if r.id not in seen:
                    seen.add(r.id)
                    if not r.functional:
                        out.append(r.id)
    return out


def compose_safe_set(model, candidates: set[str], expr: dict[str, float]) -> tuple[set[str], int]:
    """§2 step 3 -- shrink the candidate set until it disables NOTHING jointly.

    Returns (safe_set, n_dropped). Deterministic: on each collision drop the highest-expressed member
    among the genes of a disabled reaction, tie-broken lexicographically.
    """
    sel = set(candidates)
    dropped = 0
    while True:
        bad = disabled_by_set(model, sel)
        if not bad:
            return sel, dropped
        # genes in `sel` implicated in the first disabled reaction
        culprits = sorted({g.id for g in model.reactions.get_by_id(bad[0]).genes} & sel)
        if not culprits:                      # nothing removable -> drop the whole reaction's members
            sel -= {g.id for g in model.reactions.get_by_id(bad[0]).genes}
            dropped += 1
            continue
        victim = max(culprits, key=lambda g: (expr.get(g, 0.0), g))
        sel.discard(victim)
        dropped += 1


def run_arm(model, conds, genes, pctl, ref, expr_by_cond, single_gene_deletion, label):
    all_ex = tuple(conds.values())
    calls, stats = {}, {}
    gate0_ok = True
    for n, cond in enumerate(sorted(conds), 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt0 = float(wildtype_growth(model))
            sel: set[str] = set()
            ndrop = 0
            if pctl is not None:
                e = expr_by_cond[cond]
                cand = {g for g, v in e.items() if g in ref and v < ref[g]}
                cand = {g for g in cand if not gpr_disabled_reactions(model, g)}
                sel, ndrop = compose_safe_set(model, cand, e)
                for g in sel:
                    model.genes.get_by_id(g).knock_out()
            wt = float(wildtype_growth(model))
            ok = (pctl is None) or (wt > 1e-9 and abs(wt - wt0) <= WT_TOL)
            gate0_ok = gate0_ok and ok
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
                d = {}
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    d[gid] = (g != g) or (g < FRAC * wt)
            else:
                d = {g: True for g in genes}
            calls[cond] = d
            stats[cond] = {"wt_base": round(wt0, 6), "wt_gated": round(wt, 6),
                           "n_selected": len(sel), "n_dropped_by_joint_check": ndrop,
                           "gate0_ok": ok, "selected": sorted(sel)}
        print(f"  [{n:2}/{len(conds)}] {label:9} {cond[:30]:32} wt {wt0:.5f}->{wt:.5f} "
              f"sel={len(sel):4} dropped={ndrop:3} gate0={'OK' if ok else 'FAIL'}", flush=True)
    return calls, stats, gate0_ok


def confusion(calls, truth, genes, conds):
    tp = fp = fn = tn = 0
    for c in conds:
        for g in genes:
            p = bool(calls.get(c, {}).get(g, False))
            y = bool(truth.get((g, c), False))
            if p and y:
                tp += 1
            elif p:
                fp += 1
            elif y:
                fn += 1
            else:
                tn += 1
    return tp, fp, fn, tn


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=f"wiki/fba_composed_selector_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    targets = json.loads(SCREEN_ARTIFACT.read_text(encoding="utf-8"))[
        "impact_on_experimental_deficit"]["genes_isozyme_masked"]
    if len(targets) != 8:
        raise SystemExit(f"frozen target set changed size ({len(targets)}) -- refusing to run")

    model = load_model()
    conn = open_db()
    conds_all = carbon_conditions(conn, model)
    expr_by_cond, _ = build_condition_expression(conds_all)
    keys = sorted(expr_by_cond)
    conds = {k: conds_all[k] for k in keys}
    model_genes = {g.id for g in model.genes}
    # DENOMINATOR: build the two-sided gene set from the SCORED conditions (the 11 with expression),
    # which yields the 131 used by every prior artifact in this arc. Passing `conds_all` here instead
    # yields 217 -- a different denominator, and mixing the two silently breaks comparability.
    records = load_records(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    genes = [r.gene_id for r in subset]
    truth = {(r.gene_id, c): bool(y) for r in subset for c, y in r.experimental.items()}
    print(f"conditions with expression: {len(keys)}/{len(conds_all)} | genes {len(genes)} | "
          f"targets {len(targets)}")

    refs = {p: per_gene_reference(model_genes, p) for p in SENSITIVITY_PCTLS}
    print(f"per-gene reference built for {len(refs[PRIMARY_PCTL])} model genes\n")

    arms, gate0 = {}, {}
    for run_i in (1, 2):
        print(f"===== RUN {run_i} =====")
        arms.setdefault("baseline", {})[run_i] = run_arm(
            model, conds, genes, None, {}, expr_by_cond, single_gene_deletion, "baseline")
        for p in SENSITIVITY_PCTLS:
            k = f"sel_p{int(p)}"
            arms.setdefault(k, {})[run_i] = run_arm(
                model, conds, genes, p, refs[p], expr_by_cond, single_gene_deletion, f"sel p{int(p)}")
            gate0[k] = arms[k][run_i][2]

    det = {}
    for arm, runs in arms.items():
        c1, c2 = runs[1][0], runs[2][0]
        diff = sum(1 for c in conds for g in genes
                   if bool(c1.get(c, {}).get(g)) != bool(c2.get(c, {}).get(g)))
        det[arm] = {"n_cells_differing": diff, "identical": diff == 0}
    determinism_ok = all(v["identical"] for v in det.values())
    print(f"\n[determinism] {determinism_ok}")

    base = arms["baseline"][1][0]
    b_tp, b_fp, b_fn, b_tn = confusion(base, truth, genes, conds)
    print(f"[baseline] TP={b_tp} FP={b_fp} FN={b_fn} TN={b_tn}")

    results = {}
    for p in SENSITIVITY_PCTLS:
        k = f"sel_p{int(p)}"
        gc, gs, ok = arms[k][1]
        tp, fp, fn, tn = confusion(gc, truth, genes, conds)
        rec, mech_ok, mech_tot = [], 0, 0
        for t in targets:
            hits = [c for c in conds
                    if truth.get((t, c)) and gc.get(c, {}).get(t) and not base.get(c, {}).get(t)]
            if hits:
                rec.append(t)
            for c in hits:
                mech_tot += 1
                partners = {x.id for r in model.genes.get_by_id(t).reactions for x in r.genes} - {t}
                if partners & set(gs[c]["selected"]):
                    mech_ok += 1
        fp_rel = (fp - b_fp) / b_fp if b_fp else (float("inf") if fp else 0.0)
        results[k] = {
            "percentile": p, "gate0_ok": ok,
            "confusion": {"TP": tp, "FP": fp, "FN": fn, "TN": tn},
            "recall": round(tp / (tp + fn), 4) if (tp + fn) else None,
            "n_recovered_of_8": len(rec), "recovered": rec,
            "mechanism_partner_gated": f"{mech_ok}/{mech_tot}" if mech_tot else "n/a",
            "fp_relative_increase": round(fp_rel, 4),
            "guardrail_breached": fp_rel > GUARDRAIL_REL_INCREASE,
            "mean_selected": round(sum(v["n_selected"] for v in gs.values()) / max(1, len(gs)), 1),
            "total_dropped_by_joint_check": sum(v["n_dropped_by_joint_check"] for v in gs.values()),
        }
        r = results[k]
        print(f"[{k}] gate0={'OK' if ok else 'FAIL'} recovered {r['n_recovered_of_8']}/8 {rec} "
              f"mech={r['mechanism_partner_gated']} FP {fp_rel:+.1%} "
              f"breached={r['guardrail_breached']} recall {r['recall']}")

    prim = results[f"sel_p{int(PRIMARY_PCTL)}"]
    if not prim["gate0_ok"]:
        verdict = "INVALID_WILDTYPE_PERTURBED"
    elif not determinism_ok:
        verdict = "INDETERMINATE"
    elif prim["guardrail_breached"]:
        verdict = "FAILURE_GUARDRAIL_BREACHED"
    elif prim["n_recovered_of_8"] >= PRIMARY_BAR:
        mo, mt = (prim["mechanism_partner_gated"].split("/") + ["0"])[:2]
        verdict = ("H1_SUPPORTED" if mt != "0" and mo == mt
                   else "H1_SUPPORTED_MECHANISM_DISCONFIRMED")
    elif prim["n_recovered_of_8"] >= 1:
        verdict = "H1_WEAKLY_SUPPORTED"
    else:
        verdict = "H1_FALSIFIED"

    out = {
        "record": "fba-composed-selector-v1", "date": date.today().isoformat(),
        "prereg": "wiki/fba_composed_selector_prereg_2026-08-22.md", "prereg_commit": "d08a210",
        "model": model.id, "n_conditions_scored": len(keys), "conditions": keys,
        "n_conditions_without_expression": len(conds_all) - len(keys),
        "frozen_targets": targets, "n_genes": len(genes),
        "determinism": {"runs": 2, "per_arm": det, "passed": determinism_ok},
        "baseline_confusion": {"TP": b_tp, "FP": b_fp, "FN": b_fn, "TN": b_tn},
        "arms": results, "primary_percentile": PRIMARY_PCTL, "verdict": verdict,
        "per_condition": {k: arms[k][1][1] for k in results},
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {verdict}\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
