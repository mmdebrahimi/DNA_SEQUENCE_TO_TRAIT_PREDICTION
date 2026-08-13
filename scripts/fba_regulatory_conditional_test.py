"""Does a REGULATORY constraint move the conditional metric, where gap-filling could not?

    uv run python scripts/fba_regulatory_conditional_test.py

`wiki/fba_gapfill_conditional_answer_2026-08-12.md` ruled out adding reactions -- every arm, up to adding
all 1,125 donor reactions, flipped exactly zero calls. Its closing argument was that the deficit is not
missing biochemistry but missing REGULATION: the model keeps every route available in every medium, while
a real cell represses and induces. This tests that argument instead of asserting it.

The intervention is a **parsimonious (pFBA) restriction**: in each medium, solve for the minimal-flux way
to grow, then force off every gene-associated reaction carrying no flux in that solution. The surviving
route differs BY MEDIUM, which is exactly the condition-specificity the flat knockout ratios lack.

LABEL-BLIND: pFBA never sees the essentiality labels. It only answers "what is the cheapest way to grow
here?".

**THE CONTROL IS MANDATORY AND SHIPS WITH THE NUMBER.** Forcing a unique route makes many more genes look
essential, so a naive per-cell gain could be nothing but a better-matched base rate. The run therefore
scores a RATE-MATCHED RANDOM predictor -- one calling the same NUMBER of cells essential, at random, many
times -- and reports the observed value's percentile against it. A result that a rate-matched shuffle can
reach is not a result.

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
    CONDITIONS,
    apply_condition,
    conditionally_essential_genes,
    confusion_from_calls,
    constant_baselines,
    continuous_readout,
    load_labels,
    mcc,
    switch_accuracy,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402
from dna_decode.fba.nulls import margin_preserving_null  # noqa: E402
from dna_decode.fba.solver_audit import (  # noqa: E402
    audit_deletion_frame,
    merge_audits,
    suspect_cell_set,
)

ESSENTIAL_FRAC = 0.01


def _load_probe(path: str | None, root: Path, run_date: str) -> dict | None:
    """The infeasibility probe's finding, which decides how to read a non-optimal solve.

    Absent -> None, and the verdict says INDETERMINATE rather than guessing. Assuming either reading
    is exactly the mistake the pre-committed rule made.
    """
    p = Path(path) if path else root / f"wiki/fba_infeasibility_probe_{run_date}.json"
    if not p.exists():
        for cand in sorted(root.glob("wiki/fba_infeasibility_probe_*.json"), reverse=True):
            p = cand
            break
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def rate_matched_null(records, n_called_essential: int, n_draws: int = 200, seed0: int = 0,
                      conditions: tuple[str, ...] | None = None,
                      exclude_cells: set[tuple[str, str]] | None = None) -> dict:
    """Score a predictor calling `n_called_essential` cells essential AT RANDOM, `n_draws` times.

    This is the control that decides whether a per-cell gain is signal or just a better base rate. It is a
    pure function of the labels -- no model, no solver -- so it is cheap and deterministic.

    `conditions` was hardcoded to the 4 media TWICE over: explicitly here, and again by letting
    `switch_accuracy` fall through to its own 4-media default. Correct for this script, invisible to its
    tests, and the exact shape of the bug that forced the 84.8% retraction -- so both are now closed.

    `exclude_cells` abstains the same cells the measurement abstained. The random draw is taken from the
    SURVIVING cells only, so the null's base rate matches the arm it controls.
    """
    keys = sorted(conditions) if conditions is not None else sorted(CONDITIONS)
    excl = exclude_cells or set()
    subset = conditionally_essential_genes(records)
    cells = [(r.gene_id, c) for r in subset for c in keys if (r.gene_id, c) not in excl]
    if not cells or n_called_essential > len(cells):
        return {"n_draws": 0, "mean": None, "max": None}
    scores = []
    for s in range(n_draws):
        pick = set(random.Random(seed0 + s).sample(cells, n_called_essential))
        pred = {c: {r.gene_id: ((r.gene_id, c) in pick) for r in subset} for c in keys}
        v = switch_accuracy(records, pred, conditions=tuple(keys),
                            exclude_cells=excl)["per_condition_agreement"]
        if v is not None:
            scores.append(v)
    return {
        "n_draws": len(scores), "n_called_essential": n_called_essential,
        "mean": round(statistics.mean(scores), 4),
        "sd": round(statistics.pstdev(scores), 4),
        "max": round(max(scores), 4),
        "scores": scores,
    }


def score_model(model, records, gene_ids: list[str], restrict: bool,
                abstain_nonoptimal: bool = False) -> dict:
    """Score one arm. `abstain_nonoptimal` was built as an honesty switch; it turned out to be a TRAP.

    DEFAULT (False) keeps the original coding: a NaN growth becomes ratio 0.0, below ESSENTIAL_FRAC,
    therefore ESSENTIAL. cobrapy returns NaN when a solve is non-optimal, and this was SUSPECTED of
    silently coding "the solver failed" as "the gene is required".

    **That suspicion was tested and refuted** (`scripts/fba_infeasibility_probe.py`). A non-optimal
    solve here is the LP correctly reporting that the ATPM maintenance floor (lb=6.86) cannot be met
    without that gene -- deterministic on re-solve, and 38 of 39 such cells on the carbon panel are
    experimentally essential in exactly that condition. The default coding is CORRECT.

    ABSTENTION (True) therefore removes TRUE POSITIVES, not noise. Its output is a biased LOWER BOUND
    and must never be quoted as the cleaner number. It is kept because the comparison is what proves
    the point: TP collapses 56 -> 12 when the genuine essentiality calls are discarded.
    """
    from cobra.flux_analysis import pfba, single_gene_deletion  # noqa: PLC0415

    keys = tuple(sorted(CONDITIONS))
    calls, ratios, mccs, wts, n_off_by_cond = {}, {}, [], {}, {}
    audits, pfba_status = {}, {}
    for c in keys:
        with model:
            apply_condition(model, c)
            n_off = 0
            if restrict:
                sol = pfba(model)
                # Record it: an inspected-and-clean pFBA solve should be provable, not assumed.
                pfba_status[c] = str(getattr(sol, "status", "unknown"))
                for rxn in model.reactions:
                    # only gene-associated reactions: forcing off an exchange or a maintenance
                    # pseudo-reaction would change the medium rather than the regulatory state
                    if rxn.gene_reaction_rule and abs(sol.fluxes[rxn.id]) < 1e-9:
                        rxn.bounds = (0.0, 0.0)
                        n_off += 1
            wt = wildtype_growth(model)
            wts[c] = round(wt, 4)
            n_off_by_cond[c] = n_off
            d = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in gene_ids])
                audits[c] = audit_deletion_frame(res, c)
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    d[gid] = 0.0 if g != g else g / wt
            ratios[c] = d
            calls[c] = {g: v <= ESSENTIAL_FRAC for g, v in d.items()}

    excl = suspect_cell_set(audits) if abstain_nonoptimal else set()
    for c in keys:
        truth = {r.gene_id: r.experimental[c] for r in records if (r.gene_id, c) not in excl}
        pred = {g: v for g, v in calls[c].items() if (g, c) not in excl}
        cm = confusion_from_calls(truth, pred)
        mccs.append(mcc(cm))

    sw = switch_accuracy(records, calls, conditions=keys, exclude_cells=excl)
    cont = continuous_readout(records, ratios, conditions=keys)
    n_called = sum(1 for c in calls for g, v in calls[c].items() if v and (g, c) not in excl)
    tp = sum(1 for c in calls for g, v in calls[c].items()
             if v and (g, c) not in excl
             and next(r for r in records if r.gene_id == g).experimental[c])
    return {
        "exact_set_match": sw["exact_set_match"],
        "n_conditionally_essential": sw["n_conditionally_essential"],
        "n_scored_exact_set": sw["n_scored_exact_set"],
        "per_condition_agreement": sw["per_condition_agreement"],
        "mean_per_condition_mcc": round(sum(mccs) / len(mccs), 4),
        "auroc_threshold_free": cont["auroc"],
        "n_cells_called_essential": n_called,
        "tp": tp, "fp": n_called - tp,
        "precision": round(tp / n_called, 4) if n_called else None,
        "wildtype_growth": wts,
        "reactions_forced_off": n_off_by_cond,
        "abstain_nonoptimal": abstain_nonoptimal,
        "n_cells_scored": sw["n_cells_scored"],
        "n_cells_abstained": sw["n_cells_abstained"],
        "pfba_status": pfba_status,
        "solver_audit": merge_audits(audits) if audits else None,
        "_calls": calls,          # internal: the margin-preserving null shuffles these; stripped below
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", default="ecoli")
    ap.add_argument("--null-draws", type=int, default=200)
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--infeasibility-probe", default=None,
                    help="JSON from scripts/fba_infeasibility_probe.py; decides how a "
                         "non-optimal solve is read")
    a = ap.parse_args(argv)

    model = load_model(organism=a.organism)
    present = {g.id for g in model.genes}
    records = conditionally_essential_genes([r for r in load_labels() if r.gene_id in present])
    genes = [r.gene_id for r in records]
    print(f"{model.id}: scoring {len(genes)} conditionally-essential genes")

    keys = tuple(sorted(CONDITIONS))
    base = score_model(model, records, genes, restrict=False)
    reg = score_model(model, records, genes, restrict=True)
    reg_abs = score_model(model, records, genes, restrict=True, abstain_nonoptimal=True)
    for tag, s in (("BASELINE (all routes available)", base),
                   ("pFBA-RESTRICTED (only routes used in that medium)", reg),
                   ("pFBA-RESTRICTED + ABSTAIN on non-optimal solves", reg_abs)):
        print(f"   {tag}")
        print(f"      exact-set {s['exact_set_match']}/{s['n_scored_exact_set']} | per-cell "
              f"{s['per_condition_agreement']} | mean MCC {s['mean_per_condition_mcc']} | "
              f"AUROC {s['auroc_threshold_free']} | TP {s['tp']} FP {s['fp']} "
              f"(precision {s['precision']}) | abstained {s['n_cells_abstained']}")
    print(f"   reactions forced off per condition: {reg['reactions_forced_off']}")
    print(f"   pFBA solution status: {reg['pfba_status']}")

    audit = reg["solver_audit"] or {}
    f_suspect = audit.get("suspect_fraction", 0.0)
    print(f"\nSOLVER AUDIT (restricted arm): {audit.get('n_suspect_total', 0)} suspect cells of "
          f"{audit.get('n_rows_total', 0)} -> f = {f_suspect}")
    print(f"   non-optimal {audit.get('n_nonoptimal_total', 0)} | "
          f"NaN growth {audit.get('n_nan_growth_total', 0)} | "
          f"conditions affected {audit.get('n_conditions_with_nonoptimal', 0)}")

    # The null MUST be recomputed on the abstained denominator: a null over the full cell set compared
    # against a metric over a reduced one is not a control.
    excl = {tuple(x) for x in (audit.get("suspect_cells") or [])}
    nulls = constant_baselines(records, conditions=keys)
    best_const = max(g["per_condition_agreement"] for g in nulls.values())
    rm = rate_matched_null(records, reg["n_cells_called_essential"], n_draws=a.null_draws,
                           conditions=keys)
    scores = rm.pop("scores", [])
    above = sum(1 for s in scores if s >= (reg["per_condition_agreement"] or 0))
    rm["p_empirical_vs_observed"] = round(above / len(scores), 4) if scores else None

    nulls_abs = constant_baselines(records, conditions=keys, exclude_cells=excl)
    best_const_abs = max(g["per_condition_agreement"] for g in nulls_abs.values())
    rm_abs = rate_matched_null(records, reg_abs["n_cells_called_essential"], n_draws=a.null_draws,
                               conditions=keys, exclude_cells=excl)
    scores_abs = rm_abs.pop("scores", [])
    above_abs = sum(1 for s in scores_abs if s >= (reg_abs["per_condition_agreement"] or 0))
    rm_abs["p_empirical_vs_observed"] = (round(above_abs / len(scores_abs), 4)
                                         if scores_abs else None)

    # THE STRONGER NULL. The rate-matched null preserves only the grand TOTAL, so a predictor can be
    # credited for merely matching the base rate. This one preserves EVERY gene's number of essential
    # conditions AND every condition's number of essential genes, so the only way to beat it is to place
    # the calls on the RIGHT cells. This was the pFBA result's named open weakness.
    reg_calls = reg.get("_calls") or {}
    gene_ids = sorted({g for c in reg_calls for g in reg_calls[c]})

    def _score(cs):
        return switch_accuracy(records, cs, conditions=keys)["per_condition_agreement"]

    mp = (margin_preserving_null(gene_ids, keys, reg_calls, _score, n_draws=a.null_draws)
          if gene_ids else {"n_draws": 0, "mean": None})
    mp_scores = mp.pop("scores", [])
    mp_above = sum(1 for x in mp_scores if x >= (reg["per_condition_agreement"] or 0))
    mp["p_empirical_vs_observed"] = round(mp_above / len(mp_scores), 4) if mp_scores else None

    print("\nCONTROLS:")
    print(f"   best constant predictor      per-cell {best_const}")
    print(f"   rate-matched random (k={rm['n_called_essential']}) per-cell mean {rm['mean']} "
          f"sd {rm['sd']} max {rm['max']}")
    print(f"   observed pFBA-restricted     per-cell {reg['per_condition_agreement']} -> "
          f"empirical p {rm['p_empirical_vs_observed']} ({above}/{len(scores)} draws reach it)")
    print(f"   MARGIN-PRESERVING null       mean {mp['mean']} max {mp.get('max')} "
          f"p95 {mp.get('p95')}")
    print(f"   observed vs MARGIN null      {reg['per_condition_agreement']} -> "
          f"empirical p {mp['p_empirical_vs_observed']} ({mp_above}/{len(mp_scores)} draws reach it)")
    print("   -- on the ABSTAINED denominator --")
    print(f"   best constant predictor      per-cell {best_const_abs}")
    print(f"   rate-matched random          mean {rm_abs['mean']} max {rm_abs['max']}")
    print(f"   observed abstained arm       per-cell {reg_abs['per_condition_agreement']} -> "
          f"empirical p {rm_abs['p_empirical_vs_observed']}")

    # ---- VERDICT ----
    # The PRE-COMMITTED rule (authored before the run) was:
    #     A <= null_max  ->  REGULATORY_LIFT_IS_A_SOLVER_ARTIFACT
    # It fired. Its PREMISE was then FALSIFIED by scripts/fba_infeasibility_probe.py: a non-optimal
    # solve here is not a solver failure, it is the LP correctly reporting that the ATPM maintenance
    # floor (lb=6.86) cannot be met without that gene. So abstaining those cells removes the TRUE
    # POSITIVES, and the abstained arm is a biased LOWER BOUND, not a cleaner measurement.
    # The rule is therefore superseded rather than silently deleted -- both are reported.
    A = reg_abs["per_condition_agreement"] or 0.0
    n_max = rm_abs.get("max")
    published = 0.6157
    if n_max is None:
        precommitted = "REGULATORY_LIFT_INDETERMINATE_NO_NULL"
    elif A <= n_max:
        precommitted = "REGULATORY_LIFT_IS_A_SOLVER_ARTIFACT"
    elif f_suspect < 0.05 and abs(A - published) <= 0.02:
        precommitted = "REGULATORY_LIFT_CONFIRMED"
    else:
        precommitted = "REGULATORY_LIFT_PARTIALLY_SURVIVES"

    probe = _load_probe(a.infeasibility_probe, Path(__file__).resolve().parent.parent, a.date)
    probe_verdict = (probe or {}).get("verdict")
    if probe_verdict == "INFEASIBLE_IS_DETERMINISTIC_GENUINE_ESSENTIALITY":
        beats = rm["p_empirical_vs_observed"] is not None and rm["p_empirical_vs_observed"] < 0.05
        verdict = ("REGULATORY_LIFT_STANDS_ABSTENTION_IS_A_BIASED_LOWER_BOUND" if beats
                   else "NO_MOVEMENT_BEYOND_CONTROLS")
    elif probe_verdict is None:
        verdict = "INDETERMINATE_RUN_THE_INFEASIBILITY_PROBE"
    else:
        verdict = precommitted

    print(f"\nPRE-COMMITTED RULE SAID: {precommitted}")
    print(f"   (A={A} vs abstained null max={n_max}; f={f_suspect}; published={published})")
    print(f"INFEASIBILITY PROBE:     {probe_verdict}")
    print(f"VERDICT: {verdict}")
    if probe_verdict == "INFEASIBLE_IS_DETERMINISTIC_GENUINE_ESSENTIALITY":
        print("   The pre-committed rule's PREMISE is falsified: non-optimal here means the ATPM floor "
              "cannot be met, i.e. genuine essentiality. Abstention drops true positives.")

    reproduces = (reg["per_condition_agreement"] == published and reg["tp"] == 56
                  and reg["fp"] == 40)
    print(f"   default-coding arm reproduces the committed baseline: {reproduces}")
    if not reproduces:
        print("   !! WIRING CHANGED BEHAVIOUR -- the abstention arm is NOT interpretable. Fix first.")

    result = {
        "record": "fba-regulatory-conditional-recheck-v1",
        "date": a.date, "model": model.id,
        "supersedes": "wiki/fba_regulatory_conditional_test_2026-08-12.json",
        "intervention": "pFBA restriction -- force off every gene-associated reaction carrying no flux "
                        "in the parsimonious solution for that medium",
        "baseline": {k: v for k, v in base.items() if k != "_calls"},
        "regulatory": {k: v for k, v in reg.items() if k != "_calls"},
        "regulatory_abstained": {k: v for k, v in reg_abs.items() if k != "_calls"},
        "suspect_fraction_restricted_arm": f_suspect,
        "default_coding_reproduces_committed_baseline": reproduces,
        "controls": {"best_constant_predictor_per_cell": best_const, "rate_matched_random": rm,
                     "margin_preserving": mp},
        "controls_abstained": {"best_constant_predictor_per_cell": best_const_abs,
                               "rate_matched_random": rm_abs},
        "verdict": verdict,
        "verdict_precommitted_rule_said": precommitted,
        "infeasibility_probe_verdict": probe_verdict,
        "precommitted_rule_premise_falsified": (
            probe_verdict == "INFEASIBLE_IS_DETERMINISTIC_GENUINE_ESSENTIALITY"),
        "verdict_rule_precommitted": {
            "REGULATORY_LIFT_CONFIRMED": "f < 0.05 AND A > null_max AND |A - 0.6157| <= 0.02",
            "REGULATORY_LIFT_IS_A_SOLVER_ARTIFACT": "A <= null_max (lift does not survive abstention)",
            "REGULATORY_LIFT_PARTIALLY_SURVIVES": "anything else",
            "note": "Authored BEFORE the run, per the project's verdict-vs-budget discipline. A = "
                    "abstained per-cell agreement; null_max = max rate-matched draw recomputed on the "
                    "ABSTAINED denominator; f = suspect-cell fraction of the restricted arm.",
        },
        "caveats": [
            "pFBA restriction is a CRUDE proxy for regulation -- it forces off ~69% of gene-associated "
            "reactions. It is evidence about WHERE the deficit lives, not a deployable method.",
            "The gain is RECALL, not precision: precision does not improve, and the threshold-free AUROC "
            "gets WORSE, so the continuous ranking degrades even as the binary calls improve.",
            "Still only a handful of exact-set matches -- the switch is mostly still not reproduced.",
            "pFBA picks ONE optimal-flux solution; alternate optima of equal cost would force off a "
            "different route and could change which genes look essential.",
            "The DEFAULT arm codes a NaN growth (= a non-optimal solve) as ESSENTIAL. It is kept only "
            "to reproduce the committed numbers; the abstained arm is the honest one.",
            "The margin-preserving null (2026-08-13) CLOSES the previously-named weakness: it preserves "
            "every gene's and every condition's essential-call count, so matching the marginal "
            "shape earns no credit. Both nulls ship; the margin-preserving one is the binding test.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_regulatory_conditional_recheck_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
