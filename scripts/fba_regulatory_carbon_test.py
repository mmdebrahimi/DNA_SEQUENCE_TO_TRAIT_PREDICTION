"""Does the pFBA regulatory restriction survive the WIDE panel and the STRONG null?

    uv run python scripts/fba_regulatory_carbon_test.py

The pFBA restriction is the only surviving candidate direction for the conditional deficit:

  * gap-filling  -- RULED OUT (154 flips of 5,425, exact-set -1)
  * threshold retuning -- RULED OUT quantitatively (<=11% of the deficit;
    `wiki/fba_constant_gene_diagnostic_2026-08-13.md`)
  * constraining which routes are available -- the direction the diagnostic points at: 76.9% of missed
    essential cells are FLAT, i.e. the deletion changed nothing, which is redundancy

But the pFBA *method* got weaker the same day it was defended. Against the rate-matched null (which fixes
only the grand TOTAL of essential calls) it scored p = 0.0. Against the margin-preserving null (which
fixes every gene's AND every condition's call count) it scored **p = 0.06** -- suggestive, not
significant, on 268 cells across 4 media.

This re-runs the identical intervention on the **25-carbon-source panel: 217 genes x 25 conditions =
5,425 cells**, ~20x the data, scored against BOTH nulls. With that much more power the question resolves
either way:

  * clears p < 0.05 comfortably -> the 4-media p = 0.06 was an underpowered read of a real effect
  * fails again -> the regulatory *method* is not supported, and the diagnostic's DIRECTION must stand on
    its own evidence (which it does -- 76.9% flatness owes nothing to pFBA)

**PRE-REGISTERED EXPECTATION (written before the run):** it will NOT clear p < 0.05. On 4 media the
margin-preserving null's MAX exactly equalled the observed 0.6157, which reads like a small-sample
ceiling rather than a real margin. Recording this so a post-hoc rationalisation of either outcome is
visible as one.

Exit 0 always: this is an experiment, not a gate.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    conditionally_essential_genes,
    confusion_from_calls,
    constant_baselines,
    mcc,
    switch_accuracy,
)
from dna_decode.fba.fitness_browser import (  # noqa: E402
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402
from dna_decode.fba.nulls import margin_preserving_null  # noqa: E402
from dna_decode.fba.solver_audit import audit_deletion_frame, merge_audits  # noqa: E402

ESSENTIAL_FRAC = 0.01


def rate_matched_null(records, keys, n_called: int, n_draws: int = 200, seed0: int = 0) -> dict:
    """The WEAK null: fixes only the grand total of essential calls. Kept for the comparison."""
    subset = conditionally_essential_genes(records)
    cells = [(r.gene_id, c) for r in subset for c in keys]
    if not cells or n_called > len(cells):
        return {"n_draws": 0, "mean": None, "max": None}
    scores = []
    for s in range(n_draws):
        pick = set(random.Random(seed0 + s).sample(cells, n_called))
        pred = {c: {r.gene_id: ((r.gene_id, c) in pick) for r in subset} for c in keys}
        v = switch_accuracy(records, pred, conditions=keys)["per_condition_agreement"]
        if v is not None:
            scores.append(v)
    return {"n_draws": len(scores), "n_called_essential": n_called,
            "mean": round(statistics.mean(scores), 4),
            "sd": round(statistics.pstdev(scores), 4),
            "max": round(max(scores), 4), "scores": scores}


def score_arm(model, records, gene_ids, conds, keys, restrict: bool) -> dict:
    """One arm on the carbon panel. Mirrors the 4-media `score_model` exactly in coding and gating."""
    from cobra.flux_analysis import pfba, single_gene_deletion  # noqa: PLC0415

    calls, mccs, wts, n_off, pfba_status, audits = {}, [], {}, {}, {}, {}
    all_ex = tuple(conds.values())
    for n, c in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[c], all_carbon=all_ex)
            off = 0
            if restrict:
                sol = pfba(model)
                pfba_status[c] = str(getattr(sol, "status", "unknown"))
                for rxn in model.reactions:
                    # gene-associated only: forcing off an exchange would change the MEDIUM, not the
                    # regulatory state, and would silently confound the whole experiment
                    if rxn.gene_reaction_rule and abs(sol.fluxes[rxn.id]) < 1e-9:
                        rxn.bounds = (0.0, 0.0)
                        off += 1
            wt = wildtype_growth(model)
            wts[c] = round(wt, 4)
            n_off[c] = off
            d = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in gene_ids])
                audits[c] = audit_deletion_frame(res, c)
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    # NaN == a non-optimal solve == the ATPM floor cannot be met == ESSENTIAL. Verified
                    # correct 2026-08-13; do NOT "fix" this into an abstention.
                    d[gid] = (g != g) or (g < ESSENTIAL_FRAC * wt)
            calls[c] = d
        cm = confusion_from_calls({r.gene_id: r.experimental[c] for r in records}, calls[c])
        mccs.append(mcc(cm))
        print(f"   [{n:2d}/{len(keys)}] {c[:34]:36s} wt {wt:.4f} | off {off:5d} | "
              f"{sum(calls[c].values()):4d} essential", flush=True)

    sw = switch_accuracy(records, calls, conditions=keys)
    n_called = sum(1 for c in calls for v in calls[c].values() if v)
    truth = {r.gene_id: r.experimental for r in records}
    tp = sum(1 for c in calls for g, v in calls[c].items() if v and truth[g][c])
    return {
        "exact_set_match": sw["exact_set_match"],
        "n_conditionally_essential": sw["n_conditionally_essential"],
        "per_condition_agreement": sw["per_condition_agreement"],
        "mean_per_condition_mcc": round(sum(mccs) / len(mccs), 4),
        "n_cells_called_essential": n_called,
        "tp": tp, "fp": n_called - tp,
        "precision": round(tp / n_called, 4) if n_called else None,
        "reactions_forced_off": n_off,
        "pfba_status": pfba_status,
        "wildtype_growth": wts,
        "solver_audit": merge_audits(audits) if audits else None,
        "_calls": calls,
    }


def verdict_for(observed: float | None, weak_p: float | None, strong_p: float | None,
                baseline: float | None = None) -> str:
    """PURE. Two conditions, and the FIRST one is the one this function originally forgot.

    An INTERVENTION must (1) beat the arm it replaces, and (2) beat a null. The first version checked
    only (2) and fired `CONFIRMED` on a run where the restriction made per-cell agreement WORSE
    (0.7368 -> 0.6839): beating a null built from the restricted arm's OWN margins says the calls are
    well-placed *given how many were made*, which is fully compatible with making more calls being a bad
    idea. On the 4-media substrate the restriction happened to improve the baseline too, so the missing
    check was invisible.

    `baseline=None` preserves the old null-only semantics for callers that genuinely have no comparator.
    """
    if observed is None or strong_p is None:
        return "INDETERMINATE"
    if baseline is not None and observed <= baseline:
        return "REGULATORY_RESTRICTION_MAKES_IT_WORSE"
    if strong_p < 0.01:
        return "REGULATORY_LIFT_CONFIRMED_ON_WIDE_PANEL"
    if strong_p < 0.05:
        return "REGULATORY_LIFT_SURVIVES_MARGINALLY"
    if weak_p is not None and weak_p < 0.05:
        return "REGULATORY_LIFT_IS_A_WEAK_NULL_ARTIFACT"
    return "REGULATORY_LIFT_NOT_SUPPORTED"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None)
    ap.add_argument("--null-draws", type=int, default=200)
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
    records = conditionally_essential_genes(
        load_records(conn, conds, gene_filter={g.id for g in model.genes}))
    genes = [r.gene_id for r in records]
    print(f"{model.id}: {len(keys)} carbon sources | {len(genes)} conditionally-essential genes "
          f"| {len(genes) * len(keys):,} cells\n")

    print("BASELINE (all routes available)")
    base = score_arm(model, records, genes, conds, keys, restrict=False)
    print("\npFBA-RESTRICTED (only routes used in that medium)")
    reg = score_arm(model, records, genes, conds, keys, restrict=True)

    for tag, s in (("baseline  ", base), ("restricted", reg)):
        print(f"\n   {tag}: exact-set {s['exact_set_match']}/{s['n_conditionally_essential']} | "
              f"per-cell {s['per_condition_agreement']} | mean MCC {s['mean_per_condition_mcc']} | "
              f"TP {s['tp']} FP {s['fp']} (prec {s['precision']})")

    obs = reg["per_condition_agreement"]
    nulls = constant_baselines(records, conditions=keys)
    best_const = max(g["per_condition_agreement"] for g in nulls.values())

    weak = rate_matched_null(records, keys, reg["n_cells_called_essential"], n_draws=a.null_draws)
    w_scores = weak.pop("scores", [])
    weak["p_empirical_vs_observed"] = (round(sum(1 for x in w_scores if x >= (obs or 0))
                                             / len(w_scores), 4) if w_scores else None)

    reg_calls = reg.get("_calls") or {}
    gene_ids = sorted({g for c in reg_calls for g in reg_calls[c]})

    def _score(cs):
        return switch_accuracy(records, cs, conditions=keys)["per_condition_agreement"]

    print("\nmargin-preserving null (this is the binding test) ...", flush=True)
    strong = margin_preserving_null(gene_ids, keys, reg_calls, _score, n_draws=a.null_draws)
    s_scores = strong.pop("scores", [])
    s_above = sum(1 for x in s_scores if x >= (obs or 0))
    strong["p_empirical_vs_observed"] = round(s_above / len(s_scores), 4) if s_scores else None

    print("\nCONTROLS:")
    print(f"   best constant predictor  per-cell {best_const}")
    print(f"   rate-matched (WEAK)      mean {weak['mean']} max {weak['max']} -> "
          f"p {weak['p_empirical_vs_observed']}")
    print(f"   margin-preserving (STRONG) mean {strong['mean']} max {strong.get('max')} "
          f"p95 {strong.get('p95')} -> p {strong['p_empirical_vs_observed']} "
          f"({s_above}/{len(s_scores)} draws reach it)")
    print(f"   observed pFBA-restricted per-cell {obs}")

    v = verdict_for(obs, weak["p_empirical_vs_observed"], strong["p_empirical_vs_observed"],
                    baseline=base["per_condition_agreement"])
    delta = round((obs or 0) - (base["per_condition_agreement"] or 0), 4)
    print(f"\n   baseline per-cell {base['per_condition_agreement']} vs restricted {obs} "
          f"-> delta {delta:+}")
    print(f"\nVERDICT: {v}")
    print("   (4-media prior: weak-null p 0.0, strong-null p 0.06 on 268 cells)")
    print("   PRE-REGISTERED EXPECTATION was: will NOT clear p<0.05")

    result = {
        "record": "fba-regulatory-carbon-test-v1", "date": a.date, "model": model.id,
        "n_conditions": len(keys), "n_cells": len(genes) * len(keys),
        "intervention": "pFBA restriction -- force off every gene-associated reaction carrying no flux "
                        "in the parsimonious solution for that carbon source",
        "preregistered_expectation": "will NOT clear p<0.05 (recorded before the run)",
        "baseline": {k: val for k, val in base.items() if k != "_calls"},
        "regulatory": {k: val for k, val in reg.items() if k != "_calls"},
        "controls": {"best_constant_predictor_per_cell": best_const,
                     "rate_matched_random_WEAK": weak,
                     "margin_preserving_STRONG": strong},
        "verdict": v,
        "delta_per_cell_vs_baseline": delta,
        "four_media_prior": {"weak_null_p": 0.0, "strong_null_p": 0.06, "n_cells": 268,
                             "observed_per_cell": 0.6157},
        "caveats": [
            "The margin-preserving null is the BINDING test; the rate-matched one fixes only the grand "
            "total and is reported for the contrast that motivated building the stronger null.",
            "pFBA picks ONE optimal-flux solution; alternate optima of equal cost would force off a "
            "different route and could change which genes look essential. Unaddressed here.",
            "Forcing off ~69% of gene-associated reactions is a CRUDE proxy for regulation -- evidence "
            "about WHERE the deficit lives, never a deployable method.",
            "All 25 conditions are aerobic carbon sources; no oxygen axis.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_regulatory_carbon_test_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
