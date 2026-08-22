"""Expression-gated GPR: does treating an unexpressed isozyme as ABSENT recover the masked genes?

RUNS THE PRE-REGISTRATION AT `wiki/fba_expression_gated_gpr_prereg_2026-08-22.md` EXACTLY.
That document was frozen and committed (2e95b7d) while `D:` was physically disconnected, so nothing
below could have been informed by a result. Every constant here traces to a numbered section of it:

  §2  target set        the 8 genes read from the committed artifact, NOT re-derived here
  §3  intervention      mark a gene ABSENT below its condition's Nth expression percentile
  §3  primary gate 20   sensitivity range {10, 20, 30}, ALL reported, best-of is FORBIDDEN
  §4  primary endpoint  >= 4 of 8 recovered
  §4  guardrail         false positives must not rise > 20% relative -> else FAILURE
  §5  determinism       processes=1, whole pipeline twice, identical calls required

WHY THIS IS NOT E-FLUX AGAIN. E-Flux scales a reaction's BOUNDS. An `or` of two isozymes keeps the
reaction functional however the bounds move, so a single-gene deletion is untouched by it. Gating the
GPR BOOLEAN changes the operator, not the parameter: cobrapy's `Reaction.functional` evaluates
`self._gpr.eval({g.id for g in self.genes if not g.functional})`, so marking the partner non-functional
is what actually lets the target's deletion bite.

Expression normalisation, join key and condition set are IMPORTED from `fba_eflux_bridge` (§3: inherited
unchanged, including its join-key fix). No re-derivation.

Usage:
    uv run python scripts/fba_expression_gated_gpr.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

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
from scripts.fba_eflux_bridge import build_condition_expression  # noqa: E402

FRAC = 0.01
SCREEN_ARTIFACT = Path("wiki/fba_orphan_protection_2026-08-21.json")

#: §3 primary gate + pre-declared sensitivity range. Reporting the best of these as the headline is
#: FORBIDDEN by the pre-registration; the 20th is the primary and the other two are context.
PRIMARY_PCTL = 20.0
SENSITIVITY_PCTLS = (10.0, 20.0, 30.0)

#: §4 endpoints
PRIMARY_BAR = 4                 # of 8
GUARDRAIL_REL_INCREASE = 0.20   # false positives may not rise more than this, relative


def frozen_target_set() -> list[str]:
    """§2 -- read the 8 genes from the committed artifact. Never hardcoded, never re-derived."""
    d = json.loads(SCREEN_ARTIFACT.read_text(encoding="utf-8"))
    genes = d["impact_on_experimental_deficit"]["genes_isozyme_masked"]
    if len(genes) != 8:
        raise SystemExit(f"frozen target set changed size ({len(genes)} != 8) -- refusing to run")
    return sorted(genes)


def gated_genes(expr: dict[str, float], model_gene_ids: set[str], pctl: float) -> set[str]:
    """§3 -- genes to mark ABSENT: measured, in the model, below this condition's own percentile.

    An UNMEASURED gene is never gated. That mirrors the bridge's `eval_gpr` returning None rather than
    0.0 for an unmeasured gene: absence of evidence must not silently knock a gene out.
    """
    measured = {g: v for g, v in expr.items() if g in model_gene_ids}
    if not measured:
        return set()
    cut = float(np.percentile(list(measured.values()), pctl))
    return {g for g, v in measured.items() if v < cut}


def run_arm(model, conds, all_ex, genes, expr_by_cond, pctl, single_gene_deletion, label):
    """One arm. pctl=None -> baseline (no gating). Returns {cond: {gene: essential_bool}} + stats."""
    calls: dict[str, dict[str, bool]] = {}
    stats: dict[str, dict] = {}
    model_gene_ids = {g.id for g in model.genes}
    for n, cond in enumerate(sorted(conds), 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            off: set[str] = set()
            if pctl is not None:
                off = gated_genes(expr_by_cond[cond], model_gene_ids, pctl)
                for gid in off:
                    model.genes.get_by_id(gid).knock_out()
            wt = wildtype_growth(model)
            if wt > 1e-9:
                # processes=1 is PINNED (§5, and the bridge's own hard-won note: the parallel path
                # returned different numbers run-to-run with nothing else changed).
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
                d = {}
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    # NaN == genuine essentiality (ATPM floor unmet), not a solver failure. Same
                    # coding as every shipped FBA script in this repo.
                    d[gid] = (g != g) or (g < FRAC * wt)
            else:
                d = {g: True for g in genes}
            calls[cond] = d
            stats[cond] = {"wt": round(float(wt), 6), "n_gated_off": len(off),
                           "gated_off": sorted(off) if pctl is not None else []}
        print(f"  [{n}/{len(conds)}] {label:14} {cond[:34]:36} wt={stats[cond]['wt']:.5f} "
              f"gated_off={len(off):4} essential={sum(d.values()):4}", flush=True)
    return calls, stats


def confusion(calls, truth, genes, conds):
    """(TP, FP, FN, TN) over the full gene x condition grid."""
    tp = fp = fn = tn = 0
    for c in conds:
        for g in genes:
            p = bool(calls.get(c, {}).get(g, False))
            y = bool(truth.get((g, c), False))
            if p and y:
                tp += 1
            elif p and not y:
                fp += 1
            elif y:
                fn += 1
            else:
                tn += 1
    return tp, fp, fn, tn


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=f"wiki/fba_expression_gated_gpr_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    targets = frozen_target_set()
    model = load_model()
    conn = open_db()
    conds_all = carbon_conditions(conn, model)
    expr_by_cond, prov = build_condition_expression(conds_all)
    keys = sorted(expr_by_cond)
    print(f"conditions with PRECISE-1K expression: {len(keys)}/{len(conds_all)}")
    if len(keys) < 2:
        print("too few matched conditions", file=sys.stderr)
        return 2

    conds = {k: conds_all[k] for k in keys}
    all_ex = tuple(conds_all.values())
    model_genes = {g.id for g in model.genes}
    records = load_records(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    genes = [r.gene_id for r in subset]
    print(f"conditionally essential on this panel: {len(genes)}")

    # `GeneRecord.experimental` is {condition -> essential?}, already thresholded by `load_records`
    # at `fit < -2`. Read it directly rather than re-deriving -- re-thresholding here would silently
    # fork the label definition away from every other FBA artifact in this repo.
    truth = {(r.gene_id, c): bool(y)
             for r in subset for c, y in r.experimental.items()}

    missing = [t for t in targets if t not in set(genes)]
    print(f"frozen targets present on this panel: {8 - len(missing)}/8" +
          (f"  MISSING {missing}" if missing else ""))

    # ---- arms. §5: the WHOLE pipeline runs twice and the calls must be identical.
    arms: dict[str, dict] = {}
    for run_i in (1, 2):
        print(f"\n===== RUN {run_i} =====")
        base_calls, base_stats = run_arm(model, conds, all_ex, genes, expr_by_cond, None,
                                         single_gene_deletion, "baseline")
        arms.setdefault("baseline", {})[run_i] = (base_calls, base_stats)
        for p in SENSITIVITY_PCTLS:
            gc, gs = run_arm(model, conds, all_ex, genes, expr_by_cond, p,
                             single_gene_deletion, f"gated p{int(p)}")
            arms.setdefault(f"gated_p{int(p)}", {})[run_i] = (gc, gs)

    det = {}
    for arm, runs in arms.items():
        c1, c2 = runs[1][0], runs[2][0]
        diff = sum(1 for c in conds for g in genes
                   if bool(c1.get(c, {}).get(g)) != bool(c2.get(c, {}).get(g)))
        det[arm] = {"n_cells_differing": diff, "identical": diff == 0}
    determinism_ok = all(v["identical"] for v in det.values())
    print(f"\n[§5 determinism] identical calls across both runs: {determinism_ok}  {det}")

    base = arms["baseline"][1][0]
    b_tp, b_fp, b_fn, b_tn = confusion(base, truth, genes, conds)
    print(f"[baseline] TP={b_tp} FP={b_fp} FN={b_fn} TN={b_tn}")

    results = {}
    for p in SENSITIVITY_PCTLS:
        arm = f"gated_p{int(p)}"
        gc, gs = arms[arm][1]
        tp, fp, fn, tn = confusion(gc, truth, genes, conds)
        # A gate that kills the WILDTYPE sends every gene down the `wt <= 0 -> all essential` path,
        # which manufactures recoveries that mean nothing. Track those conditions explicitly.
        collapsed = {c for c in conds if gs[c]["wt"] <= 1e-9}
        recovered, evidence, recovered_feasible = [], {}, []
        for t in targets:
            hits = [c for c in conds
                    if truth.get((t, c)) and gc.get(c, {}).get(t) and not base.get(c, {}).get(t)]
            if hits:
                recovered.append(t)
                if any(c not in collapsed for c in hits):
                    recovered_feasible.append(t)
                # §2 mechanistic check: a recovery is only ON-mechanism if the PARTNER was gated off.
                evidence[t] = [{"condition": c,
                                "target_gated_off": t in set(gs[c]["gated_off"]),
                                "n_gated_off": gs[c]["n_gated_off"],
                                "wildtype_collapsed": c in collapsed} for c in hits]
        fp_rel = (fp - b_fp) / b_fp if b_fp else (float("inf") if fp else 0.0)
        results[arm] = {
            "percentile": p,
            "confusion": {"TP": tp, "FP": fp, "FN": fn, "TN": tn},
            "recall": round(tp / (tp + fn), 4) if (tp + fn) else None,
            "baseline_recall": round(b_tp / (b_tp + b_fn), 4) if (b_tp + b_fn) else None,
            "n_recovered_of_8": len(recovered),
            "recovered": recovered,
            # POST-HOC DIAGNOSTIC, not pre-registered. It restricts the primary to conditions where the
            # model can still grow, and it makes the negative STRONGER, never weaker -- which is the only
            # safe direction for an unplanned analysis. Reported alongside the pre-registered number,
            # never in place of it.
            "n_recovered_growth_feasible_only": len(recovered_feasible),
            "recovered_growth_feasible_only": recovered_feasible,
            "n_conditions_wildtype_collapsed": len(collapsed),
            "conditions_wildtype_collapsed": sorted(collapsed),
            "recovery_evidence": evidence,
            "fp_relative_increase": round(fp_rel, 4),
            "guardrail_breached": fp_rel > GUARDRAIL_REL_INCREASE,
            "mean_genes_gated_off": round(
                sum(v["n_gated_off"] for v in gs.values()) / max(1, len(gs)), 1),
        }
        r = results[arm]
        print(f"[{arm}] recovered {r['n_recovered_of_8']}/8 {recovered}  "
              f"FP {b_fp}->{fp} ({fp_rel:+.1%}) guardrail_breached={r['guardrail_breached']}  "
              f"recall {r['baseline_recall']}->{r['recall']}")
        print(f"        wildtype COLLAPSED in {len(collapsed)}/{len(conds)} conditions; "
              f"recovered in a GROWING condition: {len(recovered_feasible)}/8 {recovered_feasible}")

    # ---- §6 pre-committed verdict, on the PRIMARY gate only.
    prim = results[f"gated_p{int(PRIMARY_PCTL)}"]
    n_rec = prim["n_recovered_of_8"]
    if not determinism_ok:
        verdict = "INDETERMINATE"
    elif prim["guardrail_breached"]:
        verdict = "FAILURE_GUARDRAIL_BREACHED"
    elif n_rec >= PRIMARY_BAR:
        verdict = "H1_SUPPORTED"
    elif n_rec >= 1:
        verdict = "H1_WEAKLY_SUPPORTED"
    else:
        verdict = "H1_FALSIFIED"

    out = {
        "record": "fba-expression-gated-gpr-v1",
        "date": date.today().isoformat(),
        "prereg": "wiki/fba_expression_gated_gpr_prereg_2026-08-22.md",
        "prereg_commit": "2e95b7d",
        "model": model.id,
        "n_conditions": len(keys), "conditions": keys,
        "expression_provenance": prov,
        "frozen_targets": targets, "targets_missing_from_panel": missing,
        "n_conditionally_essential_genes": len(genes),
        "determinism": {"runs": 2, "per_arm": det, "passed": determinism_ok},
        "baseline_confusion": {"TP": b_tp, "FP": b_fp, "FN": b_fn, "TN": b_tn},
        "arms": results,
        "primary_gate_percentile": PRIMARY_PCTL,
        "primary_endpoint_bar": PRIMARY_BAR,
        "verdict": verdict,
        "caveats": [
            "PRECISE-1K is K-12 MG1655; Fitness Browser labels are Keio (BW25113 parent).",
            "A gate at the Nth percentile marks ~N% of measured genes absent BY CONSTRUCTION -- which "
            "is exactly why the false-positive guardrail is the binding endpoint, not the recovery count.",
            "Low mRNA does not prove no protein. The claim is 'expression-gated GPR recovers n of 8', "
            "never 'the isozyme is absent in the cell'.",
            "Only conditions with matched PRECISE-1K expression are scored; a target cannot be "
            "recovered in a condition that has no expression data.",
        ],
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {verdict}")
    print(f"wrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
