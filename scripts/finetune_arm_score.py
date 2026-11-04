"""Score the fine-tuned-encoder arm against its FROZEN bar.

WHY A SCORER EXISTS AT ALL. `finetune_arm_run.py` emits per-seed, per-epoch validation curves and a
frozen-budget value per seed. The frozen bar classifies on a BOOTSTRAP INTERVAL
(`classify(ci_lo, ci_hi)`), and that interval is nowhere in the run artifact. So something must
compute the paired bootstrap of per-seed frozen-budget values against the committed atlas per-seed
and run the bar's gates in its declared order. The plan that preceded this called the step a
"collection"; it is a scorer, and naming it wrongly is what left it unbuilt.

THE COMPARATOR COMES FROM THE RECONCILE ARTIFACT, NOT FROM THE HEADS BAR. This was measured rather
than assumed, and it dissolved an open question instead of answering it. The fine-tune bar carries
`comparator.median_spearman` (0.5076) but NO `atlas_per_seed` and NO `n_test_per_seed` — a paired
bootstrap needs per-seed pairs, so both are required. The heads bar has them; but its own
`comparator.source` field literally names `data/promoteratlas/local_reconcile_2026-11-01.json`, and
its copies are byte-identical to that file's. So the heads bar is a DERIVED copy and the reconcile
artifact is the PRIMARY. Reading the primary avoids coupling two frozen bars to each other, and the
fine-tune bar's own `verify()` already re-reads that artifact and refuses on drift.

THE INTERVAL WIDTH IS DERIVED FROM THE BAR, NOT CHOSEN HERE. The fine-tune bar states
`multiplicity.correction = "none"` with one declared comparison but never writes the percentiles
out — a gap in the bar I froze. One uncorrected two-sided comparison at alpha 0.05 is
`[2.5, 97.5]`, which is determined by the frozen fields; `ci_percentiles` derives it and REFUSES if
the correction is anything but `none` or the comparison count anything but 1, because either would
change the right answer. Deriving from frozen fields is defensible; inventing a number is not, and
the derivation is recorded in the artifact so a reader can check it.

THE FROZEN RULE IS READ, NEVER RECOMPUTED. Each seed's `frozen_budget.value` is taken from the run
artifact, which computed it with `finetune_arm_run.frozen_budget_value`. Re-deriving the
early-stopping rule here would put two implementations of one rule in the repo, and the failure mode
is silent disagreement on exactly the seeds that stopped early.

A BOUNDED RUN IS NOT AUTOMATICALLY VOID, and the first draft of this got that wrong. The plan said
`run_completed` should refuse a `BOUNDED_BY_WALL_CLOCK` status outright. But the horizon is 36
epochs while the budget is 7, so a run bounded during the INELIGIBLE DIAGNOSTIC TAIL has every seed
past its budget and a perfectly readable verdict. Voiding it would discard a result over the loss of
epochs the bar already declared can never be the arm's result. So `run_completed` asks whether there
is anything readable, and `budget_reached` asks whether the budget was actually reached — which is
what catches a bound that cost seeds or epochs.

Usage:
    python scripts/finetune_arm_score.py [--run-artifact PATH] [--out PATH]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.freeze_finetune_arm_bar import (
    ARM,
    CHECK_ORDER,
    OUTCOMES,
    classify,
    reachable_outcomes,
)
from scripts.freeze_expression_heads_bar import NULL_SIGNAL_RATIO
from scripts.positional_arm_score import ScoreError, median_of, paired_bootstrap

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/finetune_arm_acceptance_bar.json"
RUN_PATH = REPO / "wiki/finetune_arm_run_2026-11-13.json"
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
OUT = REPO / "wiki/finetune_arm_result_2026-11-13.json"

N_BOOT = 10_000

#: Gate -> outcome, DERIVED from the frozen bar's own check order (the string-outcome gates only).
#: Derived rather than restated so the scorer cannot run a different set than the bar froze.
GATE_SEQUENCE = tuple((g, o) for g, o in CHECK_ORDER if isinstance(o, str))


def ci_percentiles(bar: dict) -> list:
    """The bootstrap percentiles, DERIVED from the bar's frozen multiplicity block.

    The bar states `correction: none` and one declared comparison but omits the percentiles. One
    uncorrected two-sided comparison at family alpha 0.05 is [2.5, 97.5]. REFUSES on any other
    correction or comparison count, because either changes the right answer and this function would
    then be inventing one.
    """
    mult = bar.get("multiplicity") or {}
    corr = mult.get("correction")
    n = mult.get("n_declared_comparisons")
    if corr != "none" or n != 1:
        raise ScoreError(
            "the bar declares correction=%r over %r comparisons; this scorer only derives "
            "percentiles for a single UNCORRECTED comparison ([2.5, 97.5] at alpha 0.05). Any "
            "other shape needs the percentiles frozen in the bar rather than derived here."
            % (corr, n))
    return [2.5, 97.5]


# --------------------------------------------------------------------------- gates


def gate_run_completed(ctx: dict) -> dict:
    """Is there anything readable? NOT 'did the run finish its horizon'.

    A bound that landed in the ineligible diagnostic tail leaves the verdict intact, so this gate
    deliberately does not fail on `BOUNDED_BY_WALL_CLOCK` alone — `budget_reached` owns that.
    """
    run = ctx["run"]
    job_state = ctx.get("job_state")
    problems = []
    if run.get("schema") != "finetune-arm-run-v1":
        problems.append("run schema is %r" % run.get("schema"))
    if not run.get("per_seed"):
        problems.append("the run records no per-seed results")
    if job_state == "DEAD_NO_SENTINEL":
        problems.append("the launcher reports DEAD_NO_SENTINEL, so the run did not terminate and "
                        "whatever is on disk is whatever survived")
    if job_state == "REFUSED_LOCK_HELD":
        problems.append("the launcher reports exit 9: the payload never ran")
    return {
        "passed": not problems,
        "run_status": run.get("status"),
        "job_state": job_state,
        "bounded_is_not_fatal_here": (
            "a run bounded during the INELIGIBLE diagnostic tail still has every seed past the "
            "frozen budget, so its verdict is readable. budget_reached owns the case where a bound "
            "actually cost seeds or epochs."),
        "detail": ("the run terminated with a readable per-seed record (status %r)"
                   % run.get("status")) if not problems else "; ".join(problems),
    }


def gate_budget_reached(ctx: dict) -> dict:
    """Every committed seed must have reached the frozen budget, or stopped early INSIDE it.

    Early stopping inside the budget is not a shortfall: the bar records the budget as a MAXIMUM,
    so treating a converged seed as a failure would fail every well-behaved run.
    """
    run, frozen = ctx["run"], ctx["frozen_epochs"]
    per_seed = run.get("per_seed") or {}
    missing, short = [], []
    rows = {}
    for s in ctx["seeds"]:
        key = str(s)
        row = per_seed.get(key)
        if row is None:
            missing.append(key)
            continue
        fb = row.get("frozen_budget") or {}
        read = int(fb.get("n_epochs_read") or 0)
        rows[key] = {"n_epochs_read": read, "value": fb.get("value"),
                     "stopped_early": bool(fb.get("stopped_early"))}
        if fb.get("value") is None:
            short.append("%s: no frozen-budget value" % key)
        elif read < frozen and not fb.get("stopped_early"):
            short.append("%s: %d of %d epochs and did not stop early" % (key, read, frozen))
    ok = not missing and not short
    return {
        "passed": ok, "per_seed": rows,
        "n_seeds_expected": len(ctx["seeds"]), "n_seeds_present": len(rows),
        "missing_seeds": missing, "short_seeds": short,
        "rule": ("the frozen budget %d is a MAXIMUM; a seed that stopped early inside it PASSES"
                 % frozen),
        "detail": ("all %d seeds reached the frozen budget of %d epochs (or stopped early inside "
                   "it)" % (len(rows), frozen)) if ok
        else "budget not reached: missing %r; short %r" % (missing, short),
    }


def gate_split_equality(ctx: dict) -> dict:
    """The run's per-seed validation sizes must match the COMMITTED reference exactly.

    The reference is the reconcile artifact, not the heads bar — see the module note. Byte-identical
    splits are what make this arm's number comparable with the incumbent at all.
    """
    want, run = ctx["committed_n_test"], ctx["run"]
    per_seed = run.get("per_seed") or {}
    got, bad = {}, []
    for s in ctx["seeds"]:
        key = str(s)
        row = per_seed.get(key) or {}
        n_val = row.get("n_val")
        got[key] = n_val
        if n_val is None:
            bad.append((key, int(want[key]), None))
        elif int(want[key]) != int(n_val):
            bad.append((key, int(want[key]), int(n_val)))
    shared = [str(s) for s in ctx["seeds"]
              if (per_seed.get(str(s)) or {}).get("peaks_shared_between_splits")]
    if shared:
        bad.append(("peaks shared between splits on seeds", shared, None))
    return {
        "passed": not bad, "n_val_per_seed": got,
        "committed_reference": "data/promoteratlas/local_reconcile_2026-11-01.json",
        "why_this_reference": (
            "the fine-tune bar carries no per-seed sizes and the heads bar's copies are "
            "byte-identical to this artifact, which its own comparator.source names. Reading the "
            "PRIMARY avoids coupling two frozen bars."),
        "detail": ("all %d seeds reproduce the committed test-set sizes with no shared peaks"
                   % len(ctx["seeds"])) if not bad
        else "split mismatch: %r" % bad,
    }


def gate_shuffled_label_null(ctx: dict) -> dict:
    """The null must be small RELATIVE TO THE ARM'S OWN SIGNAL, per the bar's ratio rail.

    Reports the untrained baseline beside it, because a high null has two causes that look identical
    in the number and only that baseline separates them.
    """
    run = ctx["run"]
    null = run.get("null") or {}
    values = [float(v) for v in (null.get("values_at_frozen_budget") or [])]
    signal = ctx.get("arm_median")
    if not values:
        return {"passed": False, "n_replicates": 0,
                "detail": ("the run records no null replicates. A missing null is not a passing "
                           "null -- the bar requires one.")}
    if signal is None:
        return {"passed": False, "detail": "no arm median to scale the null against"}

    mx = float(np.max(np.abs(values)))
    med = float(np.median(values))
    ratio = (mx / abs(signal)) if signal else float("inf")
    ok = ratio < NULL_SIGNAL_RATIO
    untrained = null.get("untrained_baseline_per_replicate") or []
    return {
        "passed": bool(ok),
        "kind": null.get("kind"),
        "per_replicate": [round(v, 4) for v in values],
        "n_replicates": len(values),
        "median": round(med, 4), "max_abs": round(mx, 4),
        "arm_median_signal": round(float(signal), 4),
        "max_abs_over_signal": round(ratio, 4),
        "ratio_ceiling": NULL_SIGNAL_RATIO,
        "untrained_baseline_per_replicate": untrained,
        "how_to_read_a_high_null": null.get("how_to_read_a_HIGH_null"),
        "scope_limit": null.get("scope_limit"),
        "detail": ("max|null| %.4f = %.3fx the arm's %.4f (bar: <%.2fx)"
                   % (mx, ratio, signal, NULL_SIGNAL_RATIO)) if ok
        else ("null FAILED: max|null| %.4f = %.3fx the arm's %.4f. Compare against the untrained "
              "baseline %r before reading this as the fine-tune learning unjustified structure."
              % (mx, ratio, signal, untrained)),
    }


GATE_FUNCS = {
    "run_completed": gate_run_completed,
    "budget_reached": gate_budget_reached,
    "split_equality": gate_split_equality,
    "shuffled_label_null": gate_shuffled_label_null,
}


# --------------------------------------------------------------------------- scoring


def score_arm(ctx: dict) -> dict:
    """Paired bootstrap of per-seed frozen-budget values against the committed atlas per-seed."""
    run = ctx["run"]
    per_seed = run["per_seed"]
    rows, gains, wins = {}, [], 0
    for s in ctx["seeds"]:
        key = str(s)
        value = float((per_seed[key]["frozen_budget"])["value"])
        want = float(ctx["committed_atlas"][key])
        rows[key] = {"frozen_budget_spearman": round(value, 4), "atlas_committed": want,
                     "gain": round(value - want, 4),
                     "best_epoch": (per_seed[key]["frozen_budget"]).get("best_epoch"),
                     "stopped_early": (per_seed[key]["frozen_budget"]).get("stopped_early")}
        gains.append(value - want)
        wins += int(value > want)
    boot = paired_bootstrap(gains, percentiles=ctx["pct"], n_boot=N_BOOT)
    return {
        "per_seed": rows,
        "median_spearman": median_of(r["frozen_budget_spearman"] for r in rows.values()),
        "seed_wins_over_atlas": wins,
        "n_seeds": len(rows),
        "gain_bootstrap": boot,
        "frozen_budget_rule_source": (
            "frozen_budget.value is READ from the run artifact, which computed it with "
            "finetune_arm_run.frozen_budget_value. Recomputing the early-stopping rule here would "
            "put two implementations of one rule in the repo."),
    }


def verdict_from_artifact(rec: dict, bar: dict) -> tuple:
    """Re-derive the verdict from the emitted numbers plus the frozen bar, and nothing else."""
    gates = rec.get("gates", {})
    margin = bar["inherited_not_chosen"]["margin"]["value"]
    for gate, outcome in GATE_SEQUENCE:
        g = gates.get(gate)
        if g is None:
            return (outcome, "gate %r was never evaluated" % gate)
        if not g.get("passed"):
            return (outcome, g.get("detail", "gate %r failed" % gate))
    scored = rec.get("scored")
    if not scored:
        return ("FINETUNE_TIES", "no scored result was produced")
    b = scored["gain_bootstrap"]
    return classify(b["ci_lo"], b["ci_hi"], margin=margin)


# --------------------------------------------------------------------------- run


def verify(*, bar_path: Path = BAR_PATH, run_path: Path = RUN_PATH,
           reconcile_path: Path = RECONCILE) -> dict:
    for p in (bar_path, run_path, reconcile_path):
        if not Path(p).is_file():
            raise ScoreError("anchor source absent: %s" % p)

    bar = json.loads(Path(bar_path).read_text(encoding="utf-8"))
    run = json.loads(Path(run_path).read_text(encoding="utf-8"))
    rec = json.loads(Path(reconcile_path).read_text(encoding="utf-8"))

    if bar.get("schema") != "finetune-arm-acceptance-bar-v1":
        raise ScoreError("bar schema is %r" % bar.get("schema"))
    if rec.get("verdict") != "RECONCILED":
        raise ScoreError("the comparator is %r, not RECONCILED" % rec.get("verdict"))
    if rec.get("mismatched"):
        raise ScoreError("the comparator records mismatches: %r" % rec["mismatched"])

    comp = bar["inherited_not_chosen"]["comparator"]
    if abs(float(rec["atlas_median"]) - float(comp["median_spearman"])) > 1e-9:
        raise ScoreError("the bar inherits %r but the comparator now reads %r"
                         % (comp["median_spearman"], rec["atlas_median"]))
    # The reference and the bar must describe the SAME design before per-seed values are compared.
    if abs(float(rec["frac"]) - float(comp["split_frac"])) > 1e-12:
        raise ScoreError("reconcile frac %r != bar split_frac %r; these are different designs and "
                         "their per-seed numbers are not comparable"
                         % (rec["frac"], comp["split_frac"]))
    if int(rec["n_seeds"]) != int(comp["n_seeds"]):
        raise ScoreError("reconcile n_seeds %r != bar n_seeds %r"
                         % (rec["n_seeds"], comp["n_seeds"]))

    reach = reachable_outcomes()
    if reach != set(OUTCOMES):
        raise ScoreError("reachable outcomes %r != declared %r"
                         % (sorted(reach), sorted(OUTCOMES)))
    for gate, _ in GATE_SEQUENCE:
        if gate not in GATE_FUNCS:
            raise ScoreError("gate %r is in the bar's order with no implementation" % gate)
    if set(GATE_FUNCS) != {g for g, _ in GATE_SEQUENCE}:
        raise ScoreError("GATE_FUNCS %r does not match the bar's gate set %r"
                         % (sorted(GATE_FUNCS), sorted(g for g, _ in GATE_SEQUENCE)))

    return {"bar": bar, "run": run, "reconcile": rec}


def run_scoring(v: dict, *, job_state: str | None = None) -> dict:
    bar, run, rec = v["bar"], v["run"], v["reconcile"]
    frozen = int(bar["frozen_hyperparameters"]["n_epochs_budget"])
    seeds = list(range(int(rec["n_seeds"])))

    ctx = {
        "bar": bar, "run": run, "reconcile": rec, "seeds": seeds,
        "frozen_epochs": frozen,
        "committed_atlas": rec["atlas_per_seed"],
        "committed_n_test": rec["n_test_per_seed"],
        "pct": ci_percentiles(bar),
        "job_state": job_state,
    }

    out: dict = {
        "schema": "finetune-arm-result-v1",
        "arm": ARM,
        "bar": "wiki/finetune_arm_acceptance_bar.json",
        "run_artifact": "wiki/finetune_arm_run_2026-11-13.json",
        "criterion": bar["criterion"],
        "frozen_budget_epochs": frozen,
        "horizon_epochs": run.get("horizon_epochs"),
        "run_status": run.get("status"),
        "ci_percentiles": ctx["pct"],
        "ci_percentiles_derivation": (
            "DERIVED from the bar's frozen multiplicity block (correction=none over 1 declared "
            "comparison -> two-sided alpha 0.05 -> [2.5, 97.5]). The bar states the correction and "
            "the comparison count but omits the percentiles; deriving from frozen fields is "
            "checkable, inventing a number is not."),
        "comparator": {
            "source": "data/promoteratlas/local_reconcile_2026-11-01.json",
            "atlas_median": rec["atlas_median"],
            "atlas_per_seed": rec["atlas_per_seed"],
            "why_the_primary_not_the_heads_bar": (
                "the heads bar's comparator.source names this artifact and its per-seed copies are "
                "byte-identical to it, so reading the primary avoids coupling two frozen bars."),
        },
        "margin": bar["inherited_not_chosen"]["margin"]["value"],
        "gates": {},
    }

    # The null gate needs the arm's own signal to scale against, so the median is computed first --
    # but the SCORE is only published once every gate has passed, in the bar's declared order.
    prelim = None
    try:
        prelim = score_arm(ctx)
        ctx["arm_median"] = prelim["median_spearman"]
    except (KeyError, TypeError, ValueError):
        ctx["arm_median"] = None

    for gate, outcome in GATE_SEQUENCE:
        g = GATE_FUNCS[gate](ctx)
        out["gates"][gate] = g
        if not g["passed"]:
            out["verdict"], out["verdict_detail"] = outcome, g["detail"]
            out["verdict_is_rederivable_from_this_artifact"] = True
            return out

    out["scored"] = prelim
    out["arm_metric"] = prelim["median_spearman"]
    out["diagnostic_tail_per_seed"] = {
        s: (run["per_seed"][s].get("diagnostic_tail")) for s in run["per_seed"]}
    out["ineligible_diagnostics_note"] = run.get("ineligible_diagnostics_note")
    out["budget_versus_representation"] = bar["budget_versus_representation"]
    out["honest_limits"] = {
        "what_a_null_here_bounds": bar["what_a_null_here_bounds"],
        "multiplicity": bar["multiplicity"],
    }
    verdict, detail = verdict_from_artifact(out, bar)
    out["verdict"], out["verdict_detail"] = verdict, detail
    out["verdict_is_rederivable_from_this_artifact"] = True
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-artifact", default=str(RUN_PATH))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--job", default="finetune_arm")
    a = ap.parse_args(argv)

    job_state = None
    try:
        from scripts.detached_launch import classify as job_classify

        job_state = job_classify(a.job).state
    except Exception:                                          # noqa: BLE001
        pass

    try:
        v = verify(run_path=Path(a.run_artifact))
        rec = run_scoring(v, job_state=job_state)
    except ScoreError as e:
        print("REFUSED: %s" % e)
        return 2

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    for name, g in rec["gates"].items():
        print("GATE %-22s %s  %s" % (name, "PASS" if g["passed"] else "FAIL", g["detail"]))
    if rec.get("scored"):
        b = rec["scored"]["gain_bootstrap"]
        print()
        print("ARM  median %.4f vs atlas %.4f  gain %+.4f  CI [%+.4f, %+.4f]  wins %d/%d"
              % (rec["arm_metric"], rec["comparator"]["atlas_median"], b["median"],
                 b["ci_lo"], b["ci_hi"], rec["scored"]["seed_wins_over_atlas"],
                 rec["scored"]["n_seeds"]))
    print()
    print("VERDICT  %s" % rec["verdict"])
    print("         %s" % rec["verdict_detail"])
    print("-> %s" % out)
    return 0 if rec["verdict"] in OUTCOMES else 1


if __name__ == "__main__":
    raise SystemExit(main())
