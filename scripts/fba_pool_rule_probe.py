"""Are the acetate and formate zeros ONE defect? Measure the producer-detection method itself.

    uv run python scripts/fba_pool_rule_probe.py --self-check
    uv run python scripts/fba_pool_rule_probe.py

**WHY THIS EXISTS.** `wiki/fba_pool_restriction_generality_2026-11-19.json` (commit a79f9f6f) found
that restricting the candidate pool turns an unsearchable 2266-reaction space into an exhaustive
search that finds real growth-coupled designs -- 3 of 5 products. Acetate and formate returned ZERO,
and that artifact attributed formate's zero to the rule while recording acetate's as UNDIAGNOSED.

This probe tests the DETECTION METHOD rather than the products, because the rule has two independent
ways to be wrong and only one was recorded:

  1. OVER-EXCLUSION. The rule drops every reaction that "produces" the target. For a product whose
     production runs through ONE dominant fermentative route, that deletes the intervention instead
     of protecting it (formate: it cut PFL).
  2. MIS-IDENTIFICATION. Detection was done by stoichiometric SIGN (`r.metabolites[met] > 0`). A
     REVERSIBLE reaction written in the consuming direction carries a NEGATIVE coefficient for the
     product it can in fact make, so sign detection misses it entirely.

If (2) holds, the set the rule actually excluded and the set sign-detection returns are DIFFERENT
sets, and the two zeros stop being two separate puzzles.

**READ-ONLY.** Every model touch is inside `with model:`; nothing here mutates the shared cached
model. That is not hygiene -- `load_model` hands back a shared object and `reuse=False` does NOT
isolate it (model.py writes the cache unconditionally, guarded only on `path is None`), which is a
trap this repo measured the same day this probe was written.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

#: The pool and the per-target exclusions are copied VERBATIM from the committed generality run so
#: the comparison is against what actually executed, not against a re-derivation of it.
ROUTES = ["PFL", "ALCD2x", "ACALD", "FRD2", "FRD3", "ACKr", "PTAr", "PTA2",
          "POX", "MGSA", "LDH_D", "LDH_D2", "L_LACD2", "L_LACD3"]

RULE_EXCLUSIONS = {
    "EX_lac__D_e": ["LDH_D"],
    "EX_etoh_e": ["ALCD2x", "ACALD"],
    "EX_ac_e": ["ACKr", "PTAr", "PTA2", "POX"],
    "EX_succ_e": ["FRD2", "FRD3"],
    "EX_for_e": ["PFL"],
}

#: The cytosolic metabolite each exchange drains, for the sign-detection reproduction.
TARGET_METABOLITE = {
    "EX_lac__D_e": "lac__D_c",
    "EX_etoh_e": "etoh_c",
    "EX_ac_e": "ac_c",
    "EX_succ_e": "succ_c",
    "EX_for_e": "for_c",
}


# ----------------------------------------------------------------------------- pure


def compare_sets(sign_detected, rule_excluded) -> dict:
    """PURE. How do the two sets relate? The symmetric difference IS the finding.

    No model, no FBA. `sign_detected` is what `r.metabolites[met] > 0` returns; `rule_excluded` is
    what the committed run actually removed from the pool.
    """
    s, r = set(sign_detected), set(rule_excluded)
    return {
        "n_sign_detected": len(s),
        "n_rule_excluded": len(r),
        "sign_detected": sorted(s),
        "rule_excluded": sorted(r),
        "excluded_but_NOT_sign_detected": sorted(r - s),
        "sign_detected_but_NOT_excluded": sorted(s - r),
        "symmetric_difference": sorted(s ^ r),
        "sets_agree": s == r,
        "reading": (
            "the two sets AGREE, so sign detection is not the problem for this target"
            if s == r else
            "the two sets DISAGREE, so the rule excluded reactions sign detection does not call "
            "producers and/or vice versa -- the detection method is implicated, not just the breadth"),
    }


# --------------------------------------------------------------- model-dependent


def sign_detected_producers(model, met_id: str, restrict_to=None) -> list[str]:
    """Reproduce the INVALID method so the artifact can show exactly what it returns."""
    met = model.metabolites.get_by_id(met_id)
    out = [r.id for r in met.reactions if r.metabolites[met] > 0]
    if restrict_to is not None:
        keep = set(restrict_to)
        out = [r for r in out if r in keep]
    return sorted(out)


def classify_flux_role(intact_max, deleted_max, tol) -> str:
    """PURE. FOUR states, and the NaN one is the whole point.

    A NaN deleted_max means the LP was INFEASIBLE at the growth floor once this reaction was removed.
    This repo has a standing finding that an infeasible FBA deletion solve is GENUINE ESSENTIALITY
    rather than a solver failure, so it is the STRONGEST form of load-bearing -- yet `nan <= tol` is
    False in Python, so a two-state test silently files it as irrelevant. The first version of this
    probe did exactly that and reported an empty load-bearing set for all five targets.

    `reduced` exists because `collapsed` alone is too strict on a genome-scale model: deleting LDH_D
    drops lactate max from 7.256 to 0.806, a 9x reduction that a `max <= tol` rule scores as nothing.
    """
    if deleted_max is None or deleted_max != deleted_max:          # None or NaN
        return "infeasible"
    if deleted_max <= tol:
        return "collapsed"
    if intact_max is not None and intact_max == intact_max and intact_max > tol:
        if deleted_max < intact_max * 0.5:
            return "reduced"
    return "unaffected"


#: Which states mean the pool must KEEP the reaction open rather than offering it as a knockout.
LOAD_BEARING_STATES = ("infeasible", "collapsed")


def flux_role_load_bearing(model, target_rxn, candidates, growth_frac, wt_growth, tol) -> dict:
    """Which candidates are LOAD-BEARING: what does single deletion do to the achievable product?

    This is the replacement for sign detection. It asks what the reaction DOES for the target rather
    than which side of the arrow it was written on. Each probe is confined to `with model:`.
    """
    from dna_decode.fba.design import apply_knockout, product_range_at_growth

    # The intact ceiling, so a REDUCED state is measurable rather than asserted.
    _mn0, intact_max, _gr0 = product_range_at_growth(model, target_rxn, growth_frac, wt_growth)
    intact_max = None if intact_max is None else float(intact_max)

    rows = {}
    for rid in candidates:
        if rid not in {r.id for r in model.reactions}:
            rows[rid] = {"present_in_model": False}
            continue
        try:
            with model:
                apply_knockout(model, rid, "reaction")
                # THREE values: (min, max, growth). The first version of this unpacked TWO, so every
                # single reaction raised ValueError, the except below swallowed it, and `load_bearing`
                # was never set -- producing an EMPTY load-bearing set for all five targets that read
                # exactly like a clean negative and would have killed this plan's Step 2 on a bug.
                mn, mx, gr = product_range_at_growth(model, target_rxn, growth_frac, wt_growth)
        except Exception as e:                                   # noqa: BLE001
            rows[rid] = {"present_in_model": True, "error": "%s: %s" % (type(e).__name__, str(e)[:120])}
            continue
        hi_f = None if mx is None else float(mx)
        state = classify_flux_role(intact_max, hi_f, tol)
        rows[rid] = {
            "present_in_model": True,
            "intact_max_product": intact_max,
            "min_product_after_deletion": None if mn is None else float(mn),
            "max_product_after_deletion": hi_f,
            "growth_after_deletion": None if gr is None else float(gr),
            "state": state,
            "load_bearing": state in LOAD_BEARING_STATES,
        }
    return rows


#: Step 3's verdicts, FIXED HERE before any validation result was read. The plan pre-registered
#: these precisely because this module has a recorded instance of a pre-committed rule mislabelling
#: a zero, so the branch is chosen mechanically rather than argued for after the fact.
VALIDATE_VERDICTS = (
    "fix-recovers-designs",                 # formate and/or acetate now yield >= 1 coupled design
    "fix-changes-pools-but-not-outcomes",   # pools provably differ, both still zero -> GENUINE negative
    "no-pool-change",                       # constructor reproduces the old pools -> retract the diagnosis
)


def choose_validate_verdict(old_results, new_results) -> dict:
    """PURE. Decide which pre-registered branch the measurement landed in.

    `old_results` / `new_results` map target -> {"pool": int, "designs": int}. No model, no FBA.
    """
    pools_changed = {
        t: {"old": old_results[t]["pool"], "new": new_results[t]["pool"]}
        for t in new_results if t in old_results
        and old_results[t]["pool"] != new_results[t]["pool"]
    }
    was_zero = [t for t in ("EX_ac_e", "EX_for_e") if old_results.get(t, {}).get("designs") == 0]
    recovered = [t for t in was_zero if (new_results.get(t, {}).get("designs") or 0) > 0]

    if not pools_changed:
        verdict = "no-pool-change"
        reading = ("The constructor reproduced the previous pools exactly, so the diagnosis that the "
                   "old rule mis-identified producers is NOT supported by this arm and must be "
                   "retracted rather than reinterpreted.")
    elif recovered:
        verdict = "fix-recovers-designs"
        reading = ("%s now yield coupled designs where the old rule returned zero, so the rule was "
                   "the wall for those products." % ", ".join(recovered))
    else:
        verdict = "fix-changes-pools-but-not-outcomes"
        reading = ("The pools provably differ yet the previously-zero products are still zero. The "
                   "fermentative route was preserved this time, so those zeros are now a GENUINE "
                   "negative for those products rather than an artifact of the rule -- a stronger "
                   "result than the original zero, not a weaker one.")

    assert verdict in VALIDATE_VERDICTS, "verdict escaped the closed enum"
    return {
        "verdict": verdict,
        "verdict_enum_closed": list(VALIDATE_VERDICTS),
        "pools_changed": pools_changed,
        "previously_zero": was_zero,
        "recovered": recovered,
        "reading": reading,
    }


def _validate(growth_frac: float, out_path) -> int:
    """Step 3: re-run all five products on CONSTRUCTED pools and pick a pre-registered branch."""
    from dna_decode.fba.design import build_candidate_pool, find_coupled_designs, resolve_target, set_anaerobic
    from dna_decode.fba.model import load_model, wildtype_growth

    # Params read VERBATIM from the committed run, never re-derived from prose.
    committed = json.loads(
        Path("wiki/fba_pool_restriction_generality_2026-11-19.json").read_text(encoding="utf-8"))
    old_results = {t: {"pool": v["pool"], "designs": v["designs"], "best_min": v["best_min"],
                       "best_kos": v["best_kos"]}
                   for t, v in committed["results"].items()}

    model = load_model()
    assert len(model.reactions) == 2712 and len(model.genes) == 1516, "not iML1515"
    wt_aer = wildtype_growth(model)

    new_results, pool_audit = {}, {}
    with model:
        set_anaerobic(model)
        for tgt in old_results:
            trx = resolve_target(model, tgt)
            built = build_candidate_pool(model, trx, ROUTES, growth_frac=growth_frac)
            pool_audit[tgt] = {k: built[k] for k in
                               ("pool", "n_routes_in", "n_pool_out", "excluded_product_essential",
                                "retained_though_lethal", "excluded_absent_from_model", "floor_basis")}
            r = find_coupled_designs(
                model, tgt, wt_growth=wt_aer, anaerobic=True, level="reaction",
                max_knockouts=3, growth_frac=growth_frac, gene_ids=built["pool"])
            ds = r["designs"] or []
            best = ds[0] if ds else {}
            new_results[tgt] = {
                "pool": r["n_candidates_scanned"],
                "designs": r["n_coupled_designs"],
                "best_min": best.get("min_product_flux"),
                "best_growth": best.get("growth_per_h"),
                "best_kos": best.get("knockouts", []),
                "n_singles": r["n_singles_evaluated"], "n_pairs": r["n_pairs_evaluated"],
                "n_triples": r["n_triples_evaluated"],
            }
        o2_in = model.reactions.get_by_id("EX_o2_e").lower_bound
    o2_out = model.reactions.get_by_id("EX_o2_e").lower_bound

    decision = choose_validate_verdict(old_results, new_results)

    # NO-REGRESSION on the three that already worked. Checked, not assumed.
    regressions = {}
    for t in ("EX_lac__D_e", "EX_etoh_e", "EX_succ_e"):
        if old_results.get(t, {}).get("designs", 0) > 0 and (new_results.get(t, {}).get("designs") or 0) == 0:
            regressions[t] = {"old": old_results[t]["designs"], "new": new_results[t]["designs"]}

    art = {
        "schema": "fba-pool-rule-fix-v1",
        "date": str(date.today()),
        "question": ("Does the pool constructor built on the MUTANT-relative floor change the pools, "
                     "and does it recover designs for the two products that returned zero?"),
        "method_correction": ("Step 1 refuted the original flux-role method as a pure growth confound "
                             "(identical answers across five unrelated products, equal to the "
                             "growth-infeasible set). The corrected method judges the product range at "
                             "each mutant's OWN maximum growth, mirroring evaluate_knockouts."),
        "params_read_verbatim_from": "wiki/fba_pool_restriction_generality_2026-11-19.json",
        "params": {"anaerobic": True, "level": "reaction", "max_knockouts": 3,
                   "growth_frac": growth_frac},
        "old_results": old_results,
        "new_results": new_results,
        "pool_audit": pool_audit,
        "decision": decision,
        "no_regression_check": {
            "checked": ["EX_lac__D_e", "EX_etoh_e", "EX_succ_e"],
            "regressions": regressions,
            "clean": not regressions,
        },
        "non_mutation_check": {"EX_o2_e_inside": o2_in, "EX_o2_e_after": o2_out,
                               "condition_reverted": o2_out < 0.0},
        "honest_limits": [
            "Exhaustive over each CONSTRUCTED pool, never over the model's 2266 reactions.",
            "ONE model, ONE condition, ONE growth floor, depth 3.",
            "STOICHIOMETRIC predictions. A coupled design is a hypothesis for the bench, not a "
            "validated strain; FBA models no regulation, kinetics, toxicity or burden.",
        ],
    }
    out = Path(out_path) if out_path else Path("wiki/fba_pool_rule_fix_%s.json" % date.today())
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")

    r = json.loads(out.read_text(encoding="utf-8"))
    assert r["decision"]["verdict"] in VALIDATE_VERDICTS
    assert r["non_mutation_check"]["condition_reverted"], "the validation leaked its condition"
    assert sum(v["n_singles"] for v in r["new_results"].values()) > 0, "nothing was actually searched"

    print("params (verbatim): %s" % r["params"])
    for t, v in r["new_results"].items():
        o = r["old_results"][t]
        print("  %-13s pool %2d->%2d   designs %3d->%-3d  best_min %s -> %s" % (
            t, o["pool"], v["pool"], o["designs"], v["designs"], o["best_min"], v["best_min"]))
    print("\n  pools changed: %s" % list(r["decision"]["pools_changed"]))
    print("  VERDICT: %s" % r["decision"]["verdict"])
    print("  %s" % r["decision"]["reading"])
    print("  no-regression clean: %s %s" % (r["no_regression_check"]["clean"],
                                            r["no_regression_check"]["regressions"] or ""))
    print("  condition reverted: %s" % r["non_mutation_check"]["condition_reverted"])
    print("wrote %s" % out)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-check", action="store_true", help="exercise the pure comparison and exit")
    ap.add_argument("--validate", action="store_true",
                    help="Step 3: re-run all five products on CONSTRUCTED pools")
    ap.add_argument("--growth-frac", type=float, default=0.9)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.self_check:
        return _self_check()
    if args.validate:
        return _validate(args.growth_frac, args.out)

    from dna_decode.fba.design import FLUX_TOL, resolve_target, set_anaerobic
    from dna_decode.fba.model import load_model

    model = load_model()
    assert len(model.reactions) == 2712 and len(model.genes) == 1516, (
        "expected iML1515 (2712/1516), got %d/%d -- refusing to compare against the committed run "
        "using a different model" % (len(model.reactions), len(model.genes)))

    per_target = {}
    # ONE context for the whole probe: the committed run was anaerobic, so the flux-role test must be
    # too, and the condition must not outlive this function.
    with model:
        set_anaerobic(model)
        wt = model.slim_optimize()
        for tgt, excluded in RULE_EXCLUSIONS.items():
            trx = resolve_target(model, tgt)
            met_id = TARGET_METABOLITE[tgt]

            sign_all = sign_detected_producers(model, met_id)
            sign_in_pool = sign_detected_producers(model, met_id, restrict_to=ROUTES)
            cmp_pool = compare_sets(sign_in_pool, excluded)

            # Flux-role over the union: every route the rule touched plus every route sign detection
            # flagged. A reaction load-bearing by this test is one the pool must KEEP OPEN.
            union = sorted(set(ROUTES) | set(excluded))
            roles = flux_role_load_bearing(model, trx, union, args.growth_frac, wt, FLUX_TOL)
            load_bearing = sorted(k for k, v in roles.items() if v.get("load_bearing"))

            per_target[tgt] = {
                "metabolite": met_id,
                "n_sign_detected_in_whole_model": len(sign_all),
                "sign_detected_in_whole_model": sign_all,
                "comparison_within_the_pool": cmp_pool,
                "flux_role": roles,
                "load_bearing_by_flux_role": load_bearing,
                "rule_excluded_a_load_bearing_route": sorted(set(excluded) & set(load_bearing)),
                "rule_excluded_a_NON_load_bearing_route": sorted(set(excluded) - set(load_bearing)),
            }
        # THE CONTROL THAT INVALIDATES THE ARM ABOVE. If the per-target load-bearing sets come back
        # IDENTICAL across five unrelated products, the statistic is measuring the setup rather than
        # the target. This checks the obvious candidate: does the deletion simply fail to reach the
        # growth floor, independent of any product objective?
        from dna_decode.fba.design import apply_knockout as _ko
        from dna_decode.fba.design import biomass_reaction as _bio

        _ = _bio(model)
        floor_abs = args.growth_frac * wt
        growth_only = {}
        for rid in ROUTES:
            with model:
                _ko(model, rid, "reaction")
                g = model.slim_optimize()
            g = None if g is None else float(g)
            growth_only[rid] = {
                "max_growth_after_deletion": g,
                "reaches_floor": bool(g is not None and g == g and g >= floor_abs),
            }
        o2_after = model.reactions.get_by_id("EX_o2_e").lower_bound
    o2_outside = model.reactions.get_by_id("EX_o2_e").lower_bound

    # Is every per-target load-bearing set the SAME set? And does it equal the growth-infeasible set?
    lb_sets = {t: tuple(d["load_bearing_by_flux_role"]) for t, d in per_target.items()}
    all_identical = len(set(lb_sets.values())) == 1
    growth_infeasible = tuple(sorted(k for k, v in growth_only.items() if not v["reaches_floor"]))
    one_set = next(iter(set(lb_sets.values()))) if all_identical else None
    confound = {
        "per_target_load_bearing_sets": {t: list(v) for t, v in lb_sets.items()},
        "all_five_targets_identical": all_identical,
        "growth_infeasible_at_floor": list(growth_infeasible),
        "growth_only_probe": growth_only,
        "floor_absolute": floor_abs,
        "load_bearing_set_EQUALS_growth_infeasible_set": bool(
            all_identical and one_set is not None and tuple(sorted(one_set)) == growth_infeasible),
        "verdict": None,   # filled below
    }
    if confound["load_bearing_set_EQUALS_growth_infeasible_set"]:
        confound["verdict"] = "FLUX_ROLE_ARM_IS_A_PURE_GROWTH_CONFOUND"
        confound["reading"] = (
            "The per-target load-bearing set is IDENTICAL across five unrelated products AND equals "
            "exactly the set of deletions that cannot reach the growth floor. So this arm measures "
            "GROWTH ESSENTIALITY at the chosen floor, with zero product-specificity. It cannot "
            "distinguish 'load-bearing for this product' from 'required to grow at this floor'.")
        confound["root_cause"] = (
            "The floor was passed as growth_frac x WILDTYPE growth. design.py's own evaluate_knockouts "
            "docstring warns against exactly this: the floor must be growth_frac x the MUTANT'S OWN "
            "maximum, because 'a floor taken relative to the wild type is unreachable for any knockout "
            "that slows the cell' (verified in-repo 2026-08-07). Every one of these five deletions "
            "slows the cell, so each is infeasible by construction rather than by product dependence.")
        confound["consequence_for_the_plan"] = (
            "This plan's Step 2 would build build_candidate_pool on this arm. Its premise is REFUTED "
            "before any code was written: the constructor would exclude the growth-essential set for "
            "every product alike. Step 2 must be REDESIGNED, not executed -- the mutant-relative floor "
            "is the documented fix and belongs in a re-planned step.")
    elif all_identical:
        confound["verdict"] = "IDENTICAL_SETS_BUT_NOT_THE_GROWTH_SET"
        confound["reading"] = ("All five targets share one load-bearing set but it is NOT the "
                               "growth-infeasible set, so some other shared cause is in play.")
    else:
        confound["verdict"] = "TARGET_SPECIFIC"
        confound["reading"] = "The load-bearing sets differ across targets, so the arm is target-specific."

    art = {
        "schema": "fba-pool-rule-diagnosis-v1",
        "date": str(date.today()),
        "question": ("The committed generality run left acetate UNDIAGNOSED and blamed only formate on "
                     "the rule. Are both zeros ONE defect, and is the producer-DETECTION method itself "
                     "wrong rather than merely too broad?"),
        "method": ("Read-only. For each target: reproduce sign detection (r.metabolites[met] > 0), "
                   "compare it against what the committed run actually excluded, and run a FLUX-ROLE "
                   "test (does single deletion collapse max product to <= FLUX_TOL) over the union. "
                   "Anaerobic, growth_frac from the committed run, every touch inside `with model:`."),
        "growth_frac": args.growth_frac,
        "flux_tol": FLUX_TOL,
        "pool_and_exclusions_copied_verbatim_from": "wiki/fba_pool_restriction_generality_2026-11-19.json (a79f9f6f)",
        "routes": ROUTES,
        "rule_exclusions": RULE_EXCLUSIONS,
        "per_target": per_target,
        "growth_confound_control": confound,
        "non_mutation_check": {
            "EX_o2_e_inside_context": o2_after,
            "EX_o2_e_after_context": o2_outside,
            "condition_reverted": o2_outside < 0.0,
        },
    }

    out = Path(args.out) if args.out else Path("wiki/fba_pool_rule_diagnosis_%s.json" % date.today())
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")

    # Re-parse, and assert the NON-VACUITY the plan specified: if the acetate sets agreed there would
    # be nothing to fix, so an empty symmetric difference there means this probe found nothing.
    r = json.loads(out.read_text(encoding="utf-8"))
    ac = r["per_target"]["EX_ac_e"]["comparison_within_the_pool"]
    assert ac["symmetric_difference"], (
        "acetate symmetric difference is EMPTY -- sign detection and the rule agreed, so the "
        "mis-identification half of the diagnosis is NOT supported and must not be claimed")
    assert r["non_mutation_check"]["condition_reverted"], "the probe leaked its anaerobic condition"

    # NON-VACUITY ON THE FLUX-ROLE ARM, and this guard exists because its absence already cost a
    # false conclusion. The first run unpacked (min,max,growth) as two values, so EVERY reaction
    # raised ValueError, the per-reaction except swallowed it, and `load_bearing` was never set --
    # yielding an empty load-bearing set for all five targets that reads exactly like a measured
    # negative. An arm whose every cell errored must REFUSE, never report emptiness as a finding.
    errored = {t: [k for k, v in d["flux_role"].items() if "error" in v]
               for t, d in r["per_target"].items()}
    n_err = sum(len(v) for v in errored.values())
    n_cells = sum(len(d["flux_role"]) for d in r["per_target"].values())
    assert n_err == 0, (
        "%d of %d flux-role cells ERRORED (%s) -- an empty load-bearing set here would be a plumbing "
        "result, not a measurement. Refusing to write a verdict." % (n_err, n_cells, errored))
    scored = sum(1 for d in r["per_target"].values()
                 for v in d["flux_role"].values() if v.get("max_product_after_deletion") is not None)
    assert scored > 0, "no flux-role cell produced a number; nothing was actually measured"

    print("growth_frac=%s  flux_tol=%s" % (args.growth_frac, FLUX_TOL))
    for tgt, d in r["per_target"].items():
        c = d["comparison_within_the_pool"]
        print("\n%s  (%s)" % (tgt, d["metabolite"]))
        print("  sign-detected in whole model : %d" % d["n_sign_detected_in_whole_model"])
        print("  sign-detected within pool    : %s" % c["sign_detected"])
        print("  rule excluded                : %s" % c["rule_excluded"])
        print("  symmetric difference         : %s" % c["symmetric_difference"])
        print("  load-bearing by flux role    : %s" % d["load_bearing_by_flux_role"])
        print("  rule cut a LOAD-BEARING route: %s" % d["rule_excluded_a_load_bearing_route"])
    c = r["growth_confound_control"]
    print("\n=== growth-confound control ===")
    print("  all five load-bearing sets identical : %s" % c["all_five_targets_identical"])
    print("  growth-infeasible at the floor       : %s" % c["growth_infeasible_at_floor"])
    print("  load-bearing set == growth set       : %s" % c["load_bearing_set_EQUALS_growth_infeasible_set"])
    print("  VERDICT: %s" % c["verdict"])
    print("  %s" % c.get("reading", ""))
    if c.get("root_cause"):
        print("  ROOT CAUSE: %s" % c["root_cause"])

    print("\ncondition reverted: %s" % r["non_mutation_check"]["condition_reverted"])
    print("wrote %s" % out)
    return 0


def _self_check() -> int:
    """Exercise the pure comparison on synthetic sets. No model, no FBA, no substrate."""
    ok = True

    agree = compare_sets(["A", "B"], ["B", "A"])
    print("  identical sets           -> agree=%s diff=%s" % (agree["sets_agree"], agree["symmetric_difference"]))
    ok &= agree["sets_agree"] and agree["symmetric_difference"] == []

    # The acetate shape: the rule excluded things sign detection never called producers.
    ac = compare_sets(["POX"], ["ACKr", "PTAr", "PTA2", "POX"])
    print("  acetate-shaped           -> excluded_not_detected=%s" % ac["excluded_but_NOT_sign_detected"])
    ok &= ac["excluded_but_NOT_sign_detected"] == ["ACKr", "PTA2", "PTAr"]
    ok &= not ac["sets_agree"]
    ok &= "DISAGREE" in ac["reading"]

    # The formate shape: they agree, so breadth is the only charge there.
    fo = compare_sets(["PFL"], ["PFL"])
    print("  formate-shaped           -> agree=%s" % fo["sets_agree"])
    ok &= fo["sets_agree"]

    empty = compare_sets([], [])
    print("  both empty               -> agree=%s" % empty["sets_agree"])
    ok &= empty["sets_agree"] and empty["symmetric_difference"] == []

    one_sided = compare_sets(["X"], [])
    print("  detected but not excluded-> %s" % one_sided["sign_detected_but_NOT_excluded"])
    ok &= one_sided["sign_detected_but_NOT_excluded"] == ["X"]

    # classify_flux_role -- the NaN rung FIRST, because treating it as not-load-bearing is the bug
    # that made this probe report an empty set for all five targets.
    nan = float("nan")
    cases = [
        (10.0, nan, "infeasible"),     # deletion broke the LP -> strongest load-bearing
        (10.0, None, "infeasible"),    # no solution object -> same treatment
        (10.0, 0.0, "collapsed"),      # deletion zeroed the ceiling
        (10.0, 1e-9, "collapsed"),     # under tolerance counts as zero
        (10.0, 0.8, "reduced"),        # the LDH_D shape: 9x down but not zero
        (10.0, 9.9, "unaffected"),     # barely moved
        (10.0, 25.0, "unaffected"),    # deletion RAISED the ceiling (the PFL shape)
    ]
    for intact, deleted, want in cases:
        got = classify_flux_role(intact, deleted, 1e-6)
        shown = "nan" if (isinstance(deleted, float) and deleted != deleted) else deleted
        print("  classify(%s -> %-5s) -> %-11s (want %s)" % (intact, shown, got, want))
        ok &= got == want
    ok &= classify_flux_role(nan, nan, 1e-6) == "infeasible"
    ok &= "infeasible" in LOAD_BEARING_STATES and "collapsed" in LOAD_BEARING_STATES
    ok &= "reduced" not in LOAD_BEARING_STATES     # graded, reported but not pool-excluding
    print("  LOAD_BEARING_STATES     -> %s" % (LOAD_BEARING_STATES,))

    # choose_validate_verdict -- Step 3's branch picker, exercised before the expensive run so a
    # logic error cannot masquerade as a measurement.
    OLD = {"EX_lac__D_e": {"pool": 13, "designs": 63}, "EX_ac_e": {"pool": 10, "designs": 0},
           "EX_for_e": {"pool": 13, "designs": 0}}
    same = choose_validate_verdict(OLD, {k: dict(v) for k, v in OLD.items()})
    print("  verdict(identical pools) -> %s" % same["verdict"])
    ok &= same["verdict"] == "no-pool-change"

    changed_zero = choose_validate_verdict(OLD, {
        "EX_lac__D_e": {"pool": 14, "designs": 63}, "EX_ac_e": {"pool": 14, "designs": 0},
        "EX_for_e": {"pool": 14, "designs": 0}})
    print("  verdict(pools moved, still 0) -> %s" % changed_zero["verdict"])
    ok &= changed_zero["verdict"] == "fix-changes-pools-but-not-outcomes"
    ok &= "GENUINE" in changed_zero["reading"]

    recovered = choose_validate_verdict(OLD, {
        "EX_lac__D_e": {"pool": 14, "designs": 63}, "EX_ac_e": {"pool": 14, "designs": 7},
        "EX_for_e": {"pool": 14, "designs": 0}})
    print("  verdict(acetate recovered)    -> %s recovered=%s" % (
        recovered["verdict"], recovered["recovered"]))
    ok &= recovered["verdict"] == "fix-recovers-designs"
    ok &= recovered["recovered"] == ["EX_ac_e"]

    # A pool change on an ALREADY-WORKING product must not count as recovery.
    not_recovery = choose_validate_verdict(OLD, {
        "EX_lac__D_e": {"pool": 14, "designs": 99}, "EX_ac_e": {"pool": 10, "designs": 0},
        "EX_for_e": {"pool": 13, "designs": 0}})
    print("  verdict(only lactate moved)   -> %s" % not_recovery["verdict"])
    ok &= not_recovery["verdict"] == "fix-changes-pools-but-not-outcomes"
    ok &= not_recovery["recovered"] == []

    print("SELF-CHECK %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
