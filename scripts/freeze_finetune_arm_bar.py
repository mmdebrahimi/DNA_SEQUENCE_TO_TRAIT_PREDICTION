"""Freeze the acceptance bar for the fine-tuned-encoder arm, BEFORE any score exists.

The arm asks one question: does fine-tuning the PromoterAtlas encoder on this rung beat the
FROZEN encoder's median Spearman of 0.5076 by the inherited 0.02 margin, on byte-identical
splits?

**Four things are decided here and must not be decided later.**

1. **THE WALL-CLOCK CEILING, AND THEREFORE THE ROUTE.** The probe recorded
   `LOCAL_WALL_CEILING_SECONDS = 6 * 3600` and called it, in its own words, *operational, not
   statistical*: the stated reason is that an overnight run is recoverable while a multi-day
   local run "tends to be abandoned part-way". That is a claim about ATTENDANCE, not about
   statistics -- and an authorised ten-day unattended window falsifies its premise directly.
   So the ceiling here is DERIVED FROM THE DECLARED WINDOW (`derive_ceiling`), and the route
   follows from the ceiling. The window, the derived ceiling and the resulting route are all
   recorded, so a reader can check that the route followed from the window rather than from a
   preference. The rule routes BOTH WAYS and is tested in both directions: a short window sends
   this arm to Kaggle exactly as the 6-hour constant did.

2. **THE FROZEN EPOCH BUDGET IS THE ARM'S RESULT, AND LONGER HORIZONS ARE INELIGIBLE.** The run
   is carried to a long horizon because training is deterministic given a seed, so one long run
   contains every shorter budget as a prefix. The verdict is read at the frozen budget. Epochs
   beyond it are PRE-DECLARED diagnostics and can never become the arm's result -- the same shape
   as max-over-windows on the heads arm, where having declared the secondary ineligible in
   advance was the only reason the highest-scoring rung could not be re-read as a win.

3. **THE TARGET TRANSFORM, INHERITED AND NOT CHOSEN.** Trained on `log1p(expression)`, scored
   against RAW expression. This is `atlas_embed.ridge_spearman`'s behaviour, and its docstring
   states that any deviation makes the number incomparable with the committed incumbent tables.
   Training on raw expression instead would confound the OBJECTIVE TRANSFORM with the
   REPRESENTATION, which is the one thing this arm exists to isolate.

4. **WHICH NULL, AND WHY IT DIFFERS FROM THE HEADS ARM'S.** Here a SHUFFLED-LABEL null is
   informative, because the fine-tune genuinely fits: destroy the labels and you destroy the fit.
   The heads arm needed a PREDICTOR-SIDE null instead, and the reason is measured rather than
   stylistic -- a one-feature zero-shot predictor scored through `StandardScaler` + `Ridge` can
   only scale and SIGN-FLIP, so shuffling its training labels returns plus-or-minus the arm's
   entire signal (measured: +-0.8814 against a real 0.8814). A null that returns the signal is
   vacuous. With 1.8M trainable parameters that degeneracy is gone. The two arms' nulls therefore
   differ for a measured reason, and this bar records it so the difference cannot later read as
   one arm having been held to a weaker standard.

**The two readings that must not be conflated.** Beating 0.5076 shows the frozen REPRESENTATION
was the limit. Failing to beat it shows either that the representation was NOT the limit, or that
the substrate or the BUDGET was -- and only a convergence curve separates those. Today that
separation rests on a SUBSAMPLE curve (2,160 tiles). A full-cohort metric that is flat from the
frozen budget out to several times it closes the budget-limited reading EMPIRICALLY instead of by
inference. That is the entire purpose of the ineligible diagnostics.

**No multiplicity correction, and that is a difference from the heads arm, not an omission.**
There is exactly ONE declared comparison here: the median over seeds of the fine-tuned encoder
against the committed 0.5076. The heads arm is corrected at alpha/5 because it declares five
heads and takes a maximum over them. One comparison needs no divisor; applying one anyway would
widen the interval for no reason and make a real gain harder to see.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

from scripts.freeze_positional_arm_bar import (
    MARGIN,
    BarRefused,
    ceiling_arms_and_criterion,
)
from scripts.freeze_expression_heads_bar import NULL_SIGNAL_RATIO

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/finetune_arm_acceptance_bar.json"
PROBE = REPO / "wiki/finetune_arm_feasibility_2026-11-02.json"
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"

ARM = "finetuned-promoteratlas-encoder-150bp"
INCUMBENT_ARM = "promoteratlas-base-embeddings-frozen-150bp-padded-to-200"
EXPECTED_ATLAS_MEDIAN = 0.5076
PROBE_SCHEMA = "finetune-arm-feasibility-v1"

ROUTE_LOCAL = "local_cpu"
ROUTE_KAGGLE = "kaggle_t4"

#: The unattended window this bar's ceiling is derived from. Declared as DATES, and the hours are
#: computed from them -- a restated "about ten days" is a number nobody can check, and the plan's
#: round "240 hours" is nine days plus rounding.
WINDOW = {
    "name": "vacation-2026-11",
    "start": "2026-11-04T00:00:00",
    "end": "2026-11-13T00:00:00",
    "authorised_by": "the user, who stated the absence and asked for the longest run that fits it",
    "attendance": ("NONE. I execute only when invoked, so nothing inside this window can depend "
                   "on me mid-flight -- which is why a Kaggle chain is unavailable here "
                   "regardless of its compute, and why the ceiling question is live at all."),
}

#: The fraction of an unattended window that may be COMMITTED to compute.
#:
#: DERIVED FROM THE MEASURED THROUGHPUT VARIANCE, not chosen for comfort. This host's training
#: rate was measured twice and disagreed by 34% (40.7 versus 27.0 train seq/s, with the fastest
#: batch size moving between the two), so a horizon sized on today's rate can take 1.34x as long.
#: Surviving that alone requires committing no more than 1/1.34 = 0.746 of the window. The rest of
#: the headroom pays for the other named hazard: a Windows Update reboot inside a ten-day window is
#: plausible, nothing restarts the run automatically, and recovery costs whole hours. 0.60 leaves
#: both covered with the slack landing on the side where overrun costs every seed half-done.
WINDOW_COMMIT_FRACTION = 0.60

#: Measured throughput swing between two timings on this host, used in the derivation above.
MEASURED_THROUGHPUT_SWING = 0.34

OUTCOMES = (
    "FINETUNE_ADDS",
    "FINETUNE_TIES",
    "FINETUNE_NO_GAIN",
    "INDETERMINATE_RUN_DIED",
    "INDETERMINATE_BUDGET_EXHAUSTED",
    "INDETERMINATE_SPLIT_MISMATCH",
    "INDETERMINATE_NULL_FAIL",
)

#: Gate -> the outcome a failure returns, in the order the scorer must run them.
#:
#: LIVENESS FIRST, which inverts the heads arm's ordering on purpose. There, split equality came
#: first because every later number is meaningless if the test sets differ. Here a run can fail to
#: produce numbers AT ALL -- this host has six recorded long-run deaths, two of which left no trace
#: anywhere -- so the first question is whether there is anything to read. Only then does budget
#: completion matter, and only then the numbers.
CHECK_ORDER = (
    ("run_completed", "INDETERMINATE_RUN_DIED"),
    ("budget_reached", "INDETERMINATE_BUDGET_EXHAUSTED"),
    ("split_equality", "INDETERMINATE_SPLIT_MISMATCH"),
    ("shuffled_label_null", "INDETERMINATE_NULL_FAIL"),
    ("primary_decision", ("FINETUNE_ADDS", "FINETUNE_TIES", "FINETUNE_NO_GAIN")),
)


def reachable_outcomes(check_order=CHECK_ORDER) -> set[str]:
    """Every outcome the check order can actually produce.

    Asserted EQUAL to `OUTCOMES`, never merely contained in it: containment would let a declared
    outcome be unreachable, which is a bar promising a verdict it cannot deliver.
    """
    out: set[str] = set()
    for _gate, o in check_order:
        if isinstance(o, str):
            out.add(o)
        else:
            out.update(o)
    return out


def classify(ci_lo: float, ci_hi: float, *, margin: float = MARGIN) -> tuple[str, str]:
    """Total on the interval. The `(0, margin]` band is absorbed into TIES, decided in advance.

    A gain that is real but smaller than the margin is a TIE here. Saying so now is what stops it
    being re-read as a win later, which is the move a bar exists to prevent.
    """
    if ci_lo > margin:
        return ("FINETUNE_ADDS",
                "the whole interval [%+.4f, %+.4f] lies above the margin %+.4f, so a gain of at "
                "least the margin is demonstrated and the frozen representation WAS the limit"
                % (ci_lo, ci_hi, margin))
    if ci_hi < 0.0:
        return ("FINETUNE_NO_GAIN",
                "the whole interval [%+.4f, %+.4f] lies below zero, so a LOSS is demonstrated"
                % (ci_lo, ci_hi))
    return ("FINETUNE_TIES",
            "the interval [%+.4f, %+.4f] neither clears the margin %+.4f nor lies below zero, so "
            "no gain of at least the margin is demonstrated and no loss is either"
            % (ci_lo, ci_hi, margin))


def window_hours(window: dict = WINDOW) -> float:
    """Hours between the declared window's endpoints. Computed, never restated."""
    start = datetime.fromisoformat(window["start"])
    end = datetime.fromisoformat(window["end"])
    hours = (end - start).total_seconds() / 3600.0
    if hours <= 0:
        raise BarRefused("window %r ends at or before it starts" % window.get("name"))
    return hours


def derive_ceiling(hours_available: float,
                   *, fraction: float = WINDOW_COMMIT_FRACTION) -> dict:
    """The wall-clock ceiling, DERIVED from an unattended window rather than fixed.

    This is the function that replaces `LOCAL_WALL_CEILING_SECONDS = 6 * 3600`. The constant's own
    stated justification was operational -- that a multi-day local run tends to be abandoned
    part-way -- so it is a claim about attendance, and an authorised window of known length is
    direct evidence about exactly that. Replacing it with a derivation keeps the operational logic
    and removes the stale assumption about how long a window can be.

    It is NOT a disguised constant: pass tonight's ~11 hours and it returns ~6.6 h, which rejects
    this arm's 23.82 h local projection just as the old constant did. Both directions are tested.
    """
    if hours_available <= 0:
        raise BarRefused("an unattended window of %r hours cannot derive a ceiling"
                         % hours_available)
    if not 0.0 < fraction <= 1.0:
        raise BarRefused("commit fraction %r must lie in (0, 1]" % fraction)
    ceiling = hours_available * fraction
    return {
        "hours_available": round(hours_available, 2),
        "commit_fraction": fraction,
        "ceiling_hours": round(ceiling, 2),
        "derivation": (
            "%.2f h of window x %.2f committable = %.2f h. The fraction is derived from this "
            "host's MEASURED throughput variance: two timings disagreed by %.0f%%, so a horizon "
            "sized on today's rate can take %.2fx as long, which alone permits at most %.3f of "
            "the window. The remaining headroom pays for an unattended reboot, which is plausible "
            "over days and which nothing restarts automatically."
            % (hours_available, fraction, ceiling, MEASURED_THROUGHPUT_SWING * 100,
               1 + MEASURED_THROUGHPUT_SWING, 1 / (1 + MEASURED_THROUGHPUT_SWING))),
        "replaces": {
            "constant": "LOCAL_WALL_CEILING_SECONDS = 6 * 3600",
            "its_own_stated_justification": (
                "'six hours is chosen because an overnight run is recoverable while a multi-day "
                "local run often is not' -- recorded in the probe as operational, NOT statistical"),
            "why_it_no_longer_binds": (
                "the justification is a claim about ATTENDANCE, and an authorised unattended "
                "window of known length is direct evidence about attendance. The premise that a "
                "multi-day run gets abandoned part-way is false for a window the user declared in "
                "advance and asked to have filled."),
            "what_is_kept": ("the operational logic itself -- a run must fit the time actually "
                             "available with margin. Only the fixed number is replaced."),
        },
    }


def choose_route(probe: dict, ceiling: dict) -> dict:
    """Route the arm by comparing its MEASURED local projection against the derived ceiling.

    The comparison is the probe's own number, read live, never retyped.
    """
    proj = probe.get("local_projection") or {}
    try:
        local_hours = float(proj["hours_total"])
    except (KeyError, TypeError, ValueError) as e:
        raise BarRefused("the probe records no local projection to route against (%s)" % e) from e

    ceiling_hours = float(ceiling["ceiling_hours"])
    fits = local_hours <= ceiling_hours

    if fits:
        chosen = ROUTE_LOCAL
        why = ("the measured local projection %.2f h FITS the window-derived ceiling %.2f h "
               "(%.0f%% of it), so the arm runs locally and needs no mid-flight attendance. This "
               "is the decisive property, not the compute: a Kaggle chain requires an invocation "
               "to retrieve its output and nothing will invoke me inside this window."
               % (local_hours, ceiling_hours, 100.0 * local_hours / ceiling_hours))
    else:
        chosen = ROUTE_KAGGLE
        why = ("the measured local projection %.2f h EXCEEDS the window-derived ceiling %.2f h by "
               "%.1fx, so the arm cannot run locally inside this window."
               % (local_hours, ceiling_hours, local_hours / ceiling_hours))

    return {
        "chosen": chosen,
        "why": why,
        "rejected_alternative": ROUTE_KAGGLE if chosen == ROUTE_LOCAL else ROUTE_LOCAL,
        "local_projection_hours": local_hours,
        "ceiling_hours": ceiling_hours,
        "fits": fits,
        "probe_had_chosen": (probe.get("route") or {}).get("chosen"),
        "why_the_probe_chose_differently": (
            "the probe routed to %r against the FIXED 6.00 h ceiling, which this bar replaces. "
            "Nothing measured changed: the local projection is the same %.2f h in both. Only the "
            "ceiling moved, and it moved because a window was authorised."
            % ((probe.get("route") or {}).get("chosen"), local_hours)),
    }


def size_horizon(ceiling: dict, probe: dict, frozen_epochs: int) -> dict:
    """How many epochs per seed the derived ceiling actually buys.

    The frozen budget is the arm's result; this is how far past it the run can be carried for the
    ineligible diagnostics. Derived from the ceiling and the measured per-epoch cost, so a slower
    host yields a shorter horizon rather than an overrun.
    """
    proj = probe.get("local_projection") or {}
    try:
        n_seeds = int(proj["n_seeds"])
        seconds_per_epoch = float(proj["seconds_per_epoch"])
    except (KeyError, TypeError, ValueError) as e:
        raise BarRefused("the probe records no per-epoch cost to size a horizon from (%s)" % e)

    budget_s = float(ceiling["ceiling_hours"]) * 3600.0
    epochs = int(budget_s // (n_seeds * seconds_per_epoch))
    if epochs < frozen_epochs:
        raise BarRefused(
            "the derived ceiling %.2f h buys only %d epochs per seed across %d seeds, fewer than "
            "the frozen budget of %d. The arm's own verdict would not fit the window, so the bar "
            "must not be frozen against it." % (ceiling["ceiling_hours"], epochs, n_seeds,
                                                frozen_epochs))
    return {
        "frozen_budget_epochs": frozen_epochs,
        "diagnostic_horizon_epochs": epochs,
        "multiple_of_frozen_budget": round(epochs / frozen_epochs, 1),
        "n_seeds": n_seeds,
        "minutes_per_epoch": round(seconds_per_epoch / 60.0, 1),
        "projected_hours_at_horizon": round(
            n_seeds * epochs * seconds_per_epoch / 3600.0, 1),
        "why_one_long_run_not_three": (
            "training is sequential and deterministic given a seed, so a single run carried to "
            "epoch %d with the validation metric recorded at EVERY epoch contains the frozen "
            "budget's verdict as a prefix -- same seeds, same data order. Three separate runs at "
            "E, 2E and 4E epochs would cost 7E epochs per seed for strictly less information."
            % epochs),
        "what_the_extra_epochs_are_for": (
            "closing the budget-limited reading EMPIRICALLY. A metric flat from epoch %d to epoch "
            "%d (%sx the budget) rules out 'more budget would have helped' by measurement instead "
            "of by inference from a 2,160-tile subsample curve."
            % (frozen_epochs, epochs, round(epochs / frozen_epochs, 1))),
    }


# ---------------------------------------------------------------------------------------------------
# verification -- a bar anchored to retyped numbers is anchored to nothing
# ---------------------------------------------------------------------------------------------------
def verify(*, reconcile_path: Path = RECONCILE, ledger_path: Path = LEDGER,
           probe_path: Path = PROBE) -> dict:
    """Re-read every anchor live. Raises `BarRefused` on any drift. Returns what it verified."""
    for p in (reconcile_path, ledger_path, probe_path):
        if not Path(p).is_file():
            raise BarRefused("anchor source absent: %s" % p)

    rec = json.loads(Path(reconcile_path).read_text(encoding="utf-8"))
    led = json.loads(Path(ledger_path).read_text(encoding="utf-8"))
    probe = json.loads(Path(probe_path).read_text(encoding="utf-8"))

    # --- comparator -------------------------------------------------------------------------
    if rec.get("verdict") != "RECONCILED":
        raise BarRefused("the comparator is %r, not RECONCILED; every inherited number is void"
                         % rec.get("verdict"))
    if abs(float(rec["atlas_median"]) - EXPECTED_ATLAS_MEDIAN) > 1e-9:
        raise BarRefused("atlas median drifted to %r; the bar inherits %r verbatim"
                         % (rec.get("atlas_median"), EXPECTED_ATLAS_MEDIAN))
    if rec.get("mismatched"):
        raise BarRefused("the comparator records mismatches: %r" % rec["mismatched"])

    # --- arm --------------------------------------------------------------------------------
    arms, criterion = ceiling_arms_and_criterion(led)
    if ARM not in arms:
        raise BarRefused("%r is NOT declared on criterion %r; declare it under R6 first"
                         % (ARM, criterion))
    arm = arms[ARM]
    if arm.get("status") != "planned":
        raise BarRefused("arm status is %r, not 'planned'; a bar frozen after a run is not a bar"
                         % arm.get("status"))
    if arm.get("metric") is not None:
        raise BarRefused("arm already carries metric %r; freezing a bar around an existing number "
                         "is choosing the bar that the number passes" % arm.get("metric"))
    if INCUMBENT_ARM not in arms:
        raise BarRefused("the incumbent %r is absent from the ledger, so there is nothing to "
                         "compare against" % INCUMBENT_ARM)

    # --- probe ------------------------------------------------------------------------------
    if probe.get("schema") != PROBE_SCHEMA:
        raise BarRefused("probe artifact schema is %r, not %r"
                         % (probe.get("schema"), PROBE_SCHEMA))
    if not str(probe.get("verdict", "")).startswith("FEASIBLE"):
        raise BarRefused("the probe verdict is %r; the bar must not be frozen against a probe that "
                         "did not establish feasibility" % probe.get("verdict"))

    early = (probe.get("convergence") or {}).get("early_stop") or {}
    if not early.get("turned"):
        raise BarRefused(
            "the probe's convergence curve did NOT turn (%r). The budget-versus-representation "
            "distinction this bar pre-commits to rests on a curve that peaks and then declines; "
            "a monotonically rising curve means the budget was the binding constraint and a null "
            "from this arm would be uninterpretable." % early)

    hp = probe.get("hyperparameters_used") or {}
    for key in ("learning_rate", "batch_size", "patience"):
        if hp.get(key) is None:
            raise BarRefused("the probe does not record %r, so it cannot be frozen" % key)

    tt = probe.get("target_transform") or {}
    if tt.get("trained_on") != "log1p(expression)":
        raise BarRefused("the probe trained on %r, not log1p(expression); the number would not be "
                         "comparable with the committed incumbent tables"
                         % tt.get("trained_on"))

    devs = probe.get("deviations_from_pretraining_recipe")
    if devs is None:
        raise BarRefused("the probe records no deviation list from the pretraining recipe; an "
                         "UNRECORDED deviation is the thing its own deviation rule forbids")

    proj = probe.get("local_projection") or {}
    try:
        frozen_epochs = int(proj["n_epochs"])
    except (KeyError, TypeError, ValueError) as e:
        raise BarRefused("the probe records no epoch count to freeze as the budget (%s)" % e)
    if frozen_epochs <= 0:
        raise BarRefused("frozen epoch budget %r is not positive" % frozen_epochs)

    return {"reconcile": rec, "arms": arms, "criterion": criterion, "probe": probe,
            "frozen_epochs": frozen_epochs}


def build_bar(v: dict, *, window: dict = WINDOW) -> dict:
    """Pure function of the verified anchors plus the declared window. Same inputs, same bar."""
    rec, probe = v["reconcile"], v["probe"]
    frozen_epochs = v["frozen_epochs"]
    hp = probe["hyperparameters_used"]
    early = probe["convergence"]["early_stop"]

    ceiling = derive_ceiling(window_hours(window))
    route = choose_route(probe, ceiling)

    # The horizon is a property of the LOCAL route only, so it is sized only when the route is
    # local. Sizing it unconditionally conflated two different things: "this window cannot host
    # the local run" -- which is TRUE whenever the route is Kaggle, and is the very reason it is
    # Kaggle -- with "the bar cannot be frozen", which is false. Found by the short-window test.
    if route["fits"]:
        horizon = size_horizon(ceiling, probe, frozen_epochs)
    else:
        horizon = {
            "frozen_budget_epochs": frozen_epochs,
            "diagnostic_horizon_epochs": None,
            "applicable": False,
            "why_not": (
                "the arm routes to %s, so there is no local horizon to size. The window-derived "
                "ceiling %.2f h does not fit even the frozen budget locally, which is what sent "
                "the route off-host in the first place."
                % (route["chosen"], ceiling["ceiling_hours"])),
            "what_this_costs": (
                "the ineligible budget diagnostics are UNAVAILABLE on this route. A Kaggle chain "
                "needs an invocation to retrieve each segment's output, so it cannot carry a "
                "many-times-the-budget horizon through an unattended window either. On this route "
                "the budget-limited reading stays closed by INFERENCE from the subsample curve "
                "rather than by measurement -- a real weakening of the arm, recorded here so it "
                "is not discovered at reporting time."),
            "what_the_extra_epochs_are_for": (
                "nothing on this route: there are none. See what_this_costs."),
        }

    reach = reachable_outcomes()
    if reach != set(OUTCOMES):
        raise BarRefused("reachable outcomes %r != declared %r" % (sorted(reach), sorted(OUTCOMES)))

    return {
        "schema": "finetune-arm-acceptance-bar-v1",
        "frozen": "2026-11-03",
        "arm": ARM,
        "criterion": v["criterion"],
        "question": ("does fine-tuning the PromoterAtlas encoder on this rung beat the FROZEN "
                     "encoder's median Spearman %.4f by the inherited margin %.2f, on "
                     "byte-identical splits?" % (float(rec["atlas_median"]), MARGIN)),

        "inherited_not_chosen": {
            "margin": {
                "value": MARGIN,
                "derivation": ("INHERITED verbatim from the tile-chain rung's incumbent bar. Never "
                               "re-chosen here; a margin picked per-arm is a margin picked to be "
                               "cleared."),
            },
            "comparator": {
                "arm": INCUMBENT_ARM,
                "median_spearman": float(rec["atlas_median"]),
                "n_seeds": (probe.get("comparator") or {}).get("n_seeds"),
                "split_frac": (probe.get("comparator") or {}).get("split_frac"),
                "comparison": ("the fine-tuned encoder's median Spearman over seeds against the "
                               "frozen encoder's %.4f on BYTE-IDENTICAL splits. Identical splits "
                               "are what make the difference attributable to the representation "
                               "rather than to the data partition."
                               % float(rec["atlas_median"])),
                "verified_live": "RECONCILED, zero mismatched seeds, median re-read not retyped",
            },
            "target_transform": {
                **(probe.get("target_transform") or {}),
                "inherited_from": "dna_decode.glm.atlas_embed.ridge_spearman",
                "why_it_is_not_a_free_parameter": (
                    "ridge_spearman fits on log1p(expression) and scores Spearman against RAW "
                    "expression, and its docstring states that any deviation makes the number "
                    "incomparable with the committed incumbent tables. Training on raw expression "
                    "instead would confound the OBJECTIVE TRANSFORM with the REPRESENTATION, which "
                    "is the single thing this arm exists to isolate."),
            },
            "null_signal_ratio": {
                "value": NULL_SIGNAL_RATIO,
                "imported_from": "scripts.freeze_expression_heads_bar.NULL_SIGNAL_RATIO",
                "why_imported_not_restated": (
                    "the rail is scale-free -- a ceiling expressed as a fraction of the arm's own "
                    "signal -- so it transfers without re-derivation, and importing it means the "
                    "two arms cannot silently drift apart. A restated copy is a copy that goes "
                    "stale."),
            },
        },

        "ceiling": ceiling,
        "window": dict(window, hours=round(window_hours(window), 2)),
        "route": route,

        "frozen_hyperparameters": {
            "learning_rate": hp["learning_rate"],
            "batch_size": hp["batch_size"],
            "patience": hp["patience"],
            "n_epochs_budget": frozen_epochs,
            "early_stopping": {
                "monitor": "validation Spearman against RAW expression",
                "patience": hp["patience"],
                "rule": ("stop after %d non-improving epochs and report the BEST epoch, which is "
                         "what the probe's curve measured. The frozen budget %d is therefore a "
                         "MAXIMUM; a seed may stop earlier and that is not a budget failure."
                         % (hp["patience"], frozen_epochs)),
            },
            "read_from": "the probe artifact's hyperparameters_used and local_projection blocks",
            "deviations_from_pretraining_recipe": probe["deviations_from_pretraining_recipe"],
            "deviation_rule": probe.get("deviation_rule"),
        },

        "budget_versus_representation": {
            "the_distinction": (
                "beating %.4f shows the frozen REPRESENTATION was the limit. Failing to beat it "
                "shows either that the representation was NOT the limit, or that the substrate or "
                "the BUDGET was. Only a convergence curve separates those readings."
                % float(rec["atlas_median"])),
            "evidence_today": {
                "probe_field": "convergence.early_stop",
                "curve_field": "convergence.val_spearman_by_epoch",
                "turned": early["turned"],
                "best_epoch": early["best_epoch"],
                "epochs_after_best": early["epochs_after_best"],
                "n_epochs_scored": early["n_epochs_scored"],
                "scope_limit": (probe.get("convergence") or {}).get("scope"),
                "why_this_is_not_enough": (
                    "it is a SUBSAMPLE curve over %s tiles. It bounds the SHAPE of convergence, "
                    "not the full-cohort epoch count, so the budget-limited reading is closed by "
                    "INFERENCE rather than by measurement."
                    % (probe.get("convergence_sizing") or {}).get("n_subsample")),
            },
            "how_it_gets_closed_empirically": horizon["what_the_extra_epochs_are_for"],
            "ineligible_diagnostics": {
                "what": (
                    ("the per-epoch validation metric at every epoch from %d to %d"
                     % (frozen_epochs + 1, horizon["diagnostic_horizon_epochs"]))
                    if horizon.get("diagnostic_horizon_epochs")
                    else "UNAVAILABLE on this route -- see horizon.what_this_costs"),
                "status": "INELIGIBLE as the arm's result, pre-declared",
                "purpose": "answering 'would more budget have helped?' and nothing else",
                "why_pre_declared": (
                    "the same shape as max-over-windows on the heads arm, where the "
                    "highest-scoring rung was an ineligible secondary and having declared it "
                    "ineligible IN ADVANCE is the only reason that run could not be re-read as a "
                    "win. A diagnostic that may become the result is not a diagnostic."),
            },
        },

        "horizon": horizon,

        "null": {
            "kind": "shuffled_training_labels",
            "procedure": ("re-run the fine-tune with the TRAINING labels permuted, leaving the "
                          "architecture, the split, the data order and every hyperparameter "
                          "untouched, then score against the true test labels."),
            "gate": ("max|null| must fall below %.2f x the arm's own measured signal"
                     % NULL_SIGNAL_RATIO),
            "is_it_vacuous_here": False,
            "why_not": ("the fine-tune genuinely FITS -- 1.8M trainable parameters over 33,080 "
                        "training tiles -- so destroying the labels destroys the fit and the null "
                        "collapses toward zero. A null that can move is a null that tests "
                        "something."),
            "why_the_heads_arm_needed_a_DIFFERENT_null": (
                "a one-feature zero-shot predictor scored through StandardScaler + Ridge can only "
                "SCALE and SIGN-FLIP, and it takes that sign from the training labels. Shuffling "
                "them therefore returns plus-or-minus the arm's entire signal -- MEASURED at "
                "+-0.8814 against a real 0.8814, a ratio of 1.00. That null is vacuous by "
                "construction, which is why the heads arm uses a PREDICTOR-SIDE null (fixed "
                "predictor, permuted TEST labels, no refit) measured at max|null| 0.0371 against "
                "the same signal, a ratio of 0.04."),
            "why_this_is_recorded_here": (
                "so the difference between the two arms' nulls reads as a measured consequence of "
                "their predictors rather than as one arm having been held to a weaker standard. "
                "With 1.8M trainable parameters the one-feature degeneracy is simply absent."),
            "sampling_diagnostic_not_a_gate": (
                "a permutation Spearman over the test split has sampling spread near 1/sqrt(n), "
                "which on the full cohort is an order of magnitude tighter than this ratio "
                "ceiling. The tighter figure is REPORTED beside the null and is not the gate -- "
                "re-deriving a gate against the number it is about to judge is the failure this "
                "project has logged twice."),
        },

        "multiplicity": {
            "correction": "none",
            "n_declared_comparisons": 1,
            "why": ("exactly ONE comparison is declared: the median over seeds against the "
                    "committed %.4f. The heads arm is corrected at alpha/5 because it declares "
                    "five heads and takes a MAXIMUM over them. Applying a divisor to a single "
                    "comparison would widen the interval for no reason and make a real gain "
                    "harder to see." % float(rec["atlas_median"])),
            "contrast_is_deliberate": ("recorded because an uncorrected arm beside a corrected one "
                                       "invites the reading that this one was given an easier "
                                       "test. It was given a smaller hypothesis space."),
        },

        "outcomes": list(OUTCOMES),
        "check_order": [[g, o if isinstance(o, str) else list(o)] for g, o in CHECK_ORDER],
        "reachable_equals_declared": True,
        "why_liveness_is_checked_first": (
            "this inverts the heads arm's ordering on purpose. There, split equality came first "
            "because every later number is meaningless if the test sets differ. Here the run can "
            "fail to produce numbers AT ALL -- this host has six recorded long-run deaths, two of "
            "which left no event on any of 118 logs -- so the first question is whether there is "
            "anything to read at all. INDETERMINATE_RUN_DIED is the DEAD_NO_SENTINEL case from "
            "scripts/detached_launch.py, and INDETERMINATE_BUDGET_EXHAUSTED is a run that stopped "
            "cleanly but short of its frozen budget. Neither is a verdict about the arm."),

        "what_a_null_here_bounds": (
            "the fine-tuned ENCODER on THIS rung at THIS budget. It does not bound fine-tuning in "
            "general, it does not bound the architecture, and at %d epochs it cannot bound what a "
            "much larger budget would do -- which is precisely why the ineligible diagnostics "
            "exist." % frozen_epochs),

        "anchors_verified_live": {
            "comparator": str(RECONCILE.relative_to(REPO)),
            "ledger": str(LEDGER.relative_to(REPO)),
            "probe": str(PROBE.relative_to(REPO)),
            "rule": ("every number in this bar is re-read from an anchor at freeze time. A bar "
                     "anchored to retyped numbers is anchored to nothing."),
        },
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=str(OUT))
    p.add_argument("--print", action="store_true", help="print the bar without writing it")
    a = p.parse_args(argv)

    try:
        v = verify()
    except BarRefused as e:
        print("REFUSED: %s" % e)
        return 2

    bar = build_bar(v)
    text = json.dumps(bar, indent=2, sort_keys=True)
    if a.print:
        print(text)
        return 0

    out = Path(a.out)
    if out.is_file():
        print("REFUSED: %s already exists; a frozen bar is not overwritten" % out)
        return 3
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text + "\n", encoding="utf-8")
    print("froze %s" % out)
    print("  route:   %s (%s)" % (bar["route"]["chosen"], bar["route"]["why"][:80]))
    print("  ceiling: %.2f h derived from a %.2f h window"
          % (bar["ceiling"]["ceiling_hours"], bar["window"]["hours"]))
    print("  budget:  %d epochs frozen, %d-epoch horizon (%sx) as ineligible diagnostics"
          % (bar["horizon"]["frozen_budget_epochs"],
             bar["horizon"]["diagnostic_horizon_epochs"],
             bar["horizon"]["multiple_of_frozen_budget"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
