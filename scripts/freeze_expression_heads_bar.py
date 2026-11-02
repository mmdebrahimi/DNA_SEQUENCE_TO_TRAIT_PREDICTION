"""Freeze the acceptance bar for the shipped-expression-heads arm, BEFORE any head is scored.

The arm asks one question: do the five trained prediction heads that ship beside the frozen
PromoterAtlas encoder transfer, zero-shot, to the GSE144621 peak-tile rung, beating the frozen
encoder's median Spearman of 0.5076 by the inherited 0.02 margin?

**Three things are decided here and must not be decided later.**

1. **WHICH WINDOWING ROUTE DEFINES THE ARM.** The heads are geometry-locked to 86 positions against
   a 150-position substrate, so a head must be run over 65 windows, and there are two ways to do
   that which are not the same function. The choice is made here, from Step 2's measurement, with
   the tolerance DERIVED FROM THE MARGIN rather than from the observed agreement -- see
   `choose_route` for the derivation and for why the measured numbers reject the cheap route.

2. **WHICH REDUCTION IS THE PRIMARY.** Mean over the 65 windows, one primary per head. Max over
   windows is a pre-declared SECONDARY and is INELIGIBLE as the arm's result. That shape is not
   decoration: on the positional arm the highest-scoring rung was an ineligible secondary, and
   having declared it ineligible in advance is the only reason that run could not be re-read as a
   win.

3. **HOW MANY HEADS, AND THEREFORE THE CORRECTION.** Five heads were MEASURED, not the four an
   earlier note recorded, so the arm is corrected at alpha/5. The divisor is read from the probe
   artifact rather than restated, so it cannot go stale if the head set ever changes. The arm's
   metric is the MAXIMUM OVER THE DECLARED HEADS -- never over heads crossed with reductions, which
   would be ten comparisons wearing the costume of five.

**What a null here bounds, pre-committed.** All five heads were trained on OTHER libraries
(hossain2020, lafleur2022, urtecho2018, yu2021, kosuri2013). This arm therefore measures TRANSFER
to this rung. A null bounds that transfer; it does NOT bound the heads, and it must not be reported
as showing the heads are poor predictors of promoter activity in general.

**Two scoped limits, recorded rather than resolved.** The `urtecho2018` head IS circular on the
sigma70 rung (GSE108535) and must never be a comparator there; it is not circular on this substrate.
The `tslpred-kosuri2013` head predicts TRANSLATION, a different quantity from the promoter activity
this rung measures, and is included under the standing instruction to leave no stones unturned --
which costs interval width rather than buying it, since five heads mean a wider interval than four.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.freeze_positional_arm_bar import (
    MARGIN,
    BarRefused,
    bonferroni_percentiles,
    ceiling_arms_and_criterion,
)

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/expression_heads_acceptance_bar.json"
PROBE = REPO / "wiki/expression_heads_feasibility_2026-11-02.json"
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"

ARM = "shipped-expression-heads-zeroshot-150bp"

PRIMARY_REDUCTION = "mean"
SECONDARY_REDUCTION = "max"

ROUTE_TRUE = "per_window_recompute"
ROUTE_CHEAP = "full_latent_then_slice"

#: Rank agreement the cheap route must reach to be allowed to DEFINE the arm.
#:
#: DERIVED FROM THE MARGIN, not from the measurement. The verdict turns on whether a gain's
#: interval clears 0.02, and measured gains on this rung sit in the 0.00-0.05 band. Two predictors
#: agreeing at rho ~0.99 can differ in their Spearman against the labels by order 0.01, which is
#: HALF the margin -- enough for the route choice alone to move a verdict. Requiring 0.999 keeps
#: route-induced movement an order of magnitude below the margin, around 5% of it, so the arm's
#: answer is a property of the heads rather than of an optimisation.
#:
#: This threshold was written down before `choose_route` was pointed at the probe, and the probe's
#: measured minimum (0.9888) FAILS it. The cheap route is therefore not used, and the arm pays
#: 4.64 h instead of 0.24 h. Recorded because the opposite ordering -- reading the numbers, then
#: picking a tolerance they clear -- is the exact failure this project has logged twice.
ROUTE_RANK_TOLERANCE = 0.999

#: How the arm's scalar becomes a Spearman. DECLARED, because the two candidates are different
#: claims and leaving it unstated would be a free parameter inside a frozen bar.
#:
#: DIRECT is the arm's definition: `spearman(head_output[test], y[test])`, sign AS GIVEN. An
#: anti-correlated head reports a negative Spearman and that IS the answer.
#:
#: The one-feature-ridge alternative was REJECTED, and not on taste. `StandardScaler` + `Ridge` on a
#: single column can only scale and SIGN-FLIP, and it fits that sign from the TRAINING LABELS -- so
#: scoring this arm through a ridge would borrow one bit of supervision from this substrate and the
#: arm would no longer be zero-shot, which is the whole claim its own name makes. It is reported as
#: an INELIGIBLE secondary with that caveat attached.
SCORING_DIRECT = "direct_spearman_sign_as_given"
SCORING_RIDGE_INELIGIBLE = "one_feature_ridge_borrows_sign_from_train"

#: The null's magnitude ceiling, expressed as a fraction of the arm's OWN real signal.
#:
#: DERIVED FROM THE TESTED QUANTITY, not from the margin -- and the first version of this rail got
#: that wrong. A permutation Spearman over ~1,000 test tiles has sampling spread near
#: 1/sqrt(n) ~ 0.032, so holding every null below the 0.02 margin would fail on ordinary scatter and
#: say nothing. Measured: the fixed null gives max|null| 0.0371 against a real 0.8772 (ratio 0.04),
#: while the REJECTED one-feature-refit null gives max|null| 0.8814 against the same signal
#: (ratio 1.00). A ratio ceiling separates those two by an order of magnitude and is scale-free.
NULL_SIGNAL_RATIO = 0.5

OUTCOMES = (
    "HEADS_ADDS",
    "HEADS_TIES",
    "HEADS_NO_GAIN",
    "INDETERMINATE_SPLIT_MISMATCH",
    "INDETERMINATE_GEOMETRY_FAIL",
    "INDETERMINATE_COLLAPSED_HEAD",
    "INDETERMINATE_NULL_FAIL",
)

#: Gate -> the outcome a failure returns, in the order the scorer must run them. Validity before
#: decision, and split equality first of all, because every later number is meaningless if the test
#: sets differ.
CHECK_ORDER = (
    ("split_equality", "INDETERMINATE_SPLIT_MISMATCH"),
    ("geometry", "INDETERMINATE_GEOMETRY_FAIL"),
    ("head_validity", "INDETERMINATE_COLLAPSED_HEAD"),
    ("permutation_null", "INDETERMINATE_NULL_FAIL"),
    ("primary_decision", ("HEADS_ADDS", "HEADS_TIES", "HEADS_NO_GAIN")),
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
        return ("HEADS_ADDS",
                "the whole interval [%+.4f, %+.4f] lies above the margin %+.4f, so a gain of at "
                "least the margin is demonstrated" % (ci_lo, ci_hi, margin))
    if ci_hi < 0.0:
        return ("HEADS_NO_GAIN",
                "the whole interval [%+.4f, %+.4f] lies below zero, so a LOSS is demonstrated"
                % (ci_lo, ci_hi))
    return ("HEADS_TIES",
            "the interval [%+.4f, %+.4f] neither clears the margin %+.4f nor lies below zero, so "
            "no gain of at least the margin is demonstrated and no loss is either"
            % (ci_lo, ci_hi, margin))


def choose_route(probe: dict, *, tolerance: float = ROUTE_RANK_TOLERANCE) -> dict:
    """Pick the windowing route from the probe's measurement against a margin-derived tolerance.

    The rule, fixed before the numbers were read: the cheap route may DEFINE the arm only if its
    rank agreement with the faithful route reaches `tolerance` on EVERY head. One head below it is
    enough to reject, because the arm's metric is a maximum over heads and a route that distorts
    any head's ranking can distort the maximum.
    """
    routes = probe.get("routes") or {}
    if not routes:
        raise BarRefused("the probe artifact records no route measurements")

    rhos = {}
    for name, row in routes.items():
        try:
            rhos[name] = float(row["per_reduction"][PRIMARY_REDUCTION]["spearman_cheap_vs_true"])
        except (KeyError, TypeError) as e:
            raise BarRefused("head %r has no %s-reduction rank agreement to read (%s)"
                             % (name, PRIMARY_REDUCTION, e)) from e

    worst_head = min(rhos, key=lambda k: rhos[k])
    worst = rhos[worst_head]
    cost = probe.get("cost_projection") or {}
    cheap_h = (cost.get("cheap_route") or {}).get("hours_all_heads")
    true_h = (cost.get("true_route") or {}).get("hours_all_heads")

    if worst >= tolerance:
        chosen, why = ROUTE_CHEAP, (
            "every head's cheap-versus-faithful rank agreement reaches the margin-derived "
            "tolerance %.4f (worst: %s at %.4f), so the cheap route cannot move the verdict and "
            "the arm takes the %sx saving." % (tolerance, worst_head, worst,
                                               round((true_h / cheap_h), 1) if cheap_h else "?"))
    else:
        chosen, why = ROUTE_TRUE, (
            "head %s agrees with the faithful route at only %.4f, below the margin-derived "
            "tolerance %.4f, so the cheap route is NOT allowed to define the arm. The arm pays %s h "
            "instead of %s h. The cheap route is probably fine -- but 'probably fine' is not a "
            "frozen bar, and a %sx saving does not justify an unquantified risk to the verdict."
            % (worst_head, worst, tolerance, true_h, cheap_h,
               round((true_h / cheap_h), 1) if cheap_h else "?"))

    return {
        "chosen": chosen,
        "why": why,
        "tolerance": tolerance,
        "tolerance_derivation": (
            "DERIVED FROM THE MARGIN, not from these numbers. Gains on this rung sit in the "
            "0.00-0.05 band against a margin of %.2f; two predictors agreeing at rho ~0.99 can "
            "differ in Spearman-against-labels by order 0.01, half the margin. Requiring %.3f "
            "holds route-induced movement near 5%% of the margin." % (MARGIN, tolerance)),
        "measured_rank_agreement": {k: round(v, 4) for k, v in sorted(rhos.items())},
        "worst_head": worst_head,
        "worst_rank_agreement": round(worst, 4),
        "cheap_route_hours": cheap_h,
        "true_route_hours": true_h,
        "rejected_alternative": ROUTE_CHEAP if chosen == ROUTE_TRUE else ROUTE_TRUE,
        "faithful_route_is": ROUTE_TRUE,
        "why_the_true_route_is_the_faithful_one": (
            "a head is a function from an 86-position latent to a scalar, so handing it 86-position "
            "latents is what APPLYING the head means. The cheap route feeds it a pooled "
            "representation computed with context the head never saw at training time, which is an "
            "approximation introduced for speed and nothing else."),
    }


# ---------------------------------------------------------------------------------------------------
# verification -- a bar anchored to retyped numbers is anchored to nothing
# ---------------------------------------------------------------------------------------------------
def verify(*, reconcile_path: Path = RECONCILE, ledger_path: Path = LEDGER,
           probe_path: Path = PROBE) -> dict:
    for p in (reconcile_path, ledger_path, probe_path):
        if not Path(p).is_file():
            raise BarRefused("anchor source absent: %s" % p)

    rec = json.loads(Path(reconcile_path).read_text(encoding="utf-8"))
    led = json.loads(Path(ledger_path).read_text(encoding="utf-8"))
    probe = json.loads(Path(probe_path).read_text(encoding="utf-8"))

    if rec.get("verdict") != "RECONCILED":
        raise BarRefused("the comparator is %r, not RECONCILED; every inherited number is void"
                         % rec.get("verdict"))
    if abs(float(rec["atlas_median"]) - 0.5076) > 1e-9:
        raise BarRefused("atlas median drifted to %r; the bar inherits 0.5076 verbatim"
                         % rec.get("atlas_median"))
    if rec.get("mismatched"):
        raise BarRefused("the comparator records mismatches: %r" % rec["mismatched"])

    arms, criterion = ceiling_arms_and_criterion(led)
    if ARM not in arms:
        raise BarRefused("%r is NOT declared on criterion %r; declare it under R6 first"
                         % (ARM, criterion))
    arm = arms[ARM]
    if arm.get("status") != "planned":
        raise BarRefused("arm status is %r, not 'planned'; a bar frozen after a run is not a bar"
                         % arm.get("status"))
    if arm.get("metric") is not None:
        raise BarRefused("arm already carries metric %r" % arm.get("metric"))

    if probe.get("schema") != "expression-heads-feasibility-v1":
        raise BarRefused("probe artifact schema is %r" % probe.get("schema"))
    if not str(probe.get("verdict", "")).startswith("FEASIBLE"):
        raise BarRefused("the probe verdict is %r; the bar must not be frozen against a probe that "
                         "did not establish feasibility" % probe.get("verdict"))
    heads = probe.get("heads") or {}
    if not heads:
        raise BarRefused("the probe artifact names no heads")
    if set(heads) != set(probe.get("routes") or {}):
        raise BarRefused("the probe measured %d heads but loaded %d; a head cannot be allowed to "
                         "go missing between loading and measurement"
                         % (len(probe.get("routes") or {}), len(heads)))
    for name, row in heads.items():
        if row.get("base_tensors_identical") != row.get("base_tensors_checked"):
            raise BarRefused("head %r has drifted base weights; this would be a FINE-TUNED-encoder "
                             "arm, not this one" % name)
        if list(row.get("fc1_shape") or []) != [64, 11008]:
            raise BarRefused("head %r has fc1 %r, not [64, 11008]" % (name, row.get("fc1_shape")))

    return {"reconcile": rec, "arms": arms, "criterion": criterion, "probe": probe,
            "n_heads": len(heads), "head_names": sorted(heads)}


def build_bar(v: dict) -> dict:
    """Pure function of the verified anchors. Same inputs, same bar."""
    rec, probe = v["reconcile"], v["probe"]
    n_heads = v["n_heads"]
    route = choose_route(probe)
    pct = bonferroni_percentiles(n_heads)

    return {
        "schema": "expression-heads-acceptance-bar-v1",
        "frozen": "2026-11-02",
        "arm": ARM,
        "criterion": v["criterion"],
        "question": ("do the five shipped expression heads, applied zero-shot, beat the frozen "
                     "encoder's median Spearman %.4f by the inherited margin %.2f on byte-identical "
                     "splits?" % (float(rec["atlas_median"]), MARGIN)),
        "inherited_not_chosen": {
            "margin": {"value": MARGIN,
                       "derivation": ("INHERITED verbatim from the tile-chain rung's incumbent bar; "
                                      "not re-chosen for this arm, because a per-arm margin is a "
                                      "free parameter")},
            "comparator": {"atlas_median": float(rec["atlas_median"]),
                           "atlas_per_seed": rec["atlas_per_seed"],
                           "source": "data/promoteratlas/local_reconcile_2026-11-01.json",
                           "must_be_recomputed_in_run": False},
            "split": {"n_seeds": len(rec["atlas_per_seed"]), "frac": rec["frac"],
                      "n_test_per_seed": rec["n_test_per_seed"],
                      "rule": "leave_peak_out_indices, byte-identical to the incumbent's"},
        },
        "windowing": dict(route, n_windows=probe["geometry"]["n_windows"],
                          geometry=probe["geometry"]["arithmetic"]),
        "reductions": {
            "primary": PRIMARY_REDUCTION,
            "secondary": SECONDARY_REDUCTION,
            "secondary_is_eligible_as_the_arms_result": False,
            "why": ("one primary per head. The max-over-windows reduction is declared now so it "
                    "can be reported without being promotable: on the positional arm the "
                    "highest-scoring rung was an ineligible secondary, and having said so in "
                    "advance is the only reason that result could not be re-read as a win."),
        },
        "multiplicity": {
            "n_primaries": n_heads,
            "primaries": v["head_names"],
            "correction": "BONFERRONI over the %d declared heads (alpha/%d)" % (n_heads, n_heads),
            "ci_percentiles": pct,
            "divisor_is_read_from_the_probe": True,
            "why_not_heads_times_reductions": (
                "the arm's metric is the MAXIMUM OVER HEADS at the primary reduction only. "
                "Correcting over heads x reductions would be ten comparisons wearing the costume "
                "of five, and would also make the ineligible secondary count toward the verdict."),
            "divisor_moved_from_four_to_five": (
                "an earlier note recorded four heads; five were measured. The correction is "
                "therefore alpha/5, which makes the arm HARDER to pass than alpha/4 -- the honest "
                "direction for an inherited miscount."),
        },
        "metric_rule": ("the arm's metric is max over the %d declared heads of the median Spearman "
                        "at the %s reduction; never over a ladder, never across arms"
                        % (n_heads, PRIMARY_REDUCTION)),
        "outcomes": list(OUTCOMES),
        "check_order": [[g, o if isinstance(o, str) else list(o)] for g, o in CHECK_ORDER],
        "reachable_equals_declared": sorted(reachable_outcomes()) == sorted(OUTCOMES),
        "validity_controls": {
            "split_equality": "the test sets must match the committed n_test_per_seed exactly",
            "geometry": "every head's fc1 must be (64, 11008) and the window count must be 65",
            "head_validity": ("each head must produce finite, NON-CONSTANT output across the "
                              "cohort. A collapsed head would score as a clean null rather than as "
                              "a fault, which is the failure mode that most resembles a finding."),
            "permutation_null": {
                "form": ("FIXED predictor, PERMUTED TEST labels, NO refit: "
                         "spearman(head_output[test], permute(y[test])) over the seed set."),
                "condition_median": "|median null| < %.2f" % MARGIN,
                "condition_magnitude": ("max|null| < %.2f x |real median Spearman| -- this is the "
                                        "condition that catches a null which cannot destroy the "
                                        "signal it is meant to destroy" % NULL_SIGNAL_RATIO),
                "why_not_the_inherited_shuffled_TRAINING_label_null": (
                    "The sibling arms shuffle TRAINING labels and refit, which works because their "
                    "primaries are 1,280-feature matrices that a shuffled fit genuinely destroys. "
                    "This arm's primary is ONE feature, where StandardScaler + Ridge can only scale "
                    "and SIGN-FLIP, so shuffling the training labels cannot destroy the ranking -- "
                    "only its sign. MEASURED: every such null returns +/-0.8814 against a real "
                    "0.8814, i.e. the arm's entire signal with a random sign, so the gate either "
                    "FAILS for a reason unrelated to the heads or PASSES by sign cancellation. "
                    "Neither is a control. The 128-feature case measured 0.01-0.12 (median 0.0088) "
                    "on the same data, which is why the sibling arms were unaffected."),
                "measured_replacement": ("the permuted-test-label null measured |median| 0.0015, "
                                         "max|null| 0.0371 against a real 0.8772 -- noise-sized, "
                                         "and consistent with the 1/sqrt(n_test) ~ 0.032 sampling "
                                         "scale a permutation Spearman has by construction"),
            },
        },
        "scoring": {
            "definition": SCORING_DIRECT,
            "formula": "spearman(head_output[test], y[test]) with the sign AS GIVEN",
            "ineligible_alternative": SCORING_RIDGE_INELIGIBLE,
            "why_the_ridge_form_is_ineligible": (
                "a one-feature ridge fits the SIGN from the training labels, so the arm would "
                "borrow supervision from this substrate and stop being zero-shot -- which is the "
                "claim its own name makes. Reported beside the primary, never as it."),
            "an_anticorrelated_head_reports_a_NEGATIVE_spearman": True,
            "comparability_note": (
                "the incumbent 0.5076 comes from ridge_spearman, which fits on log1p(expression) "
                "and scores Spearman against RAW expression. This arm fits nothing, so only the "
                "SCORING half applies: Spearman against RAW expression, identically."),
        },
        "zero_shot_framing": probe["zero_shot_framing"],
        "forecast": {
            "predicted_outcome": "HEADS_NO_GAIN",
            "p_heads_no_gain": 0.70,
            "p_heads_ties": 0.22,
            "p_heads_adds": 0.08,
            "reasoning": ("these heads were fitted on different libraries at 86 bp and are applied "
                          "here with no fitting on this substrate at all, while the incumbent is a "
                          "ridge FITTED on this rung. A zero-shot predictor losing to a fitted one "
                          "is the expected direction."),
            "my_forecast_record_in_this_project": ("5 of 9. Do not defer to this prior; it is "
                                                   "recorded so the bar can be scored on it."),
        },
        "pre_commitments": [
            "a sub-margin gain returns HEADS_TIES, decided here and NOT re-readable as a win",
            "the %s reduction is reported but is NOT eligible as the arm's result"
            % SECONDARY_REDUCTION,
            "the arm's metric is the max over heads at the primary reduction ONLY",
            "a null bounds TRANSFER to this rung, NOT the heads",
            "the route is %s and is NOT to be switched after seeing a score" % route["chosen"],
            "a collapsed head yields INDETERMINATE_COLLAPSED_HEAD, never a scored null",
            "the arm is scored by DIRECT Spearman with the sign as given; the ridge form is "
            "reported but is NOT eligible as the arm's result",
            "the null is PREDICTOR-side-valid (fixed predictor, permuted test labels) and must "
            "satisfy BOTH its median and its magnitude condition",
        ],
        "amendments_before_first_score": [
            {
                "amended": "2026-11-02",
                "what": ("replaced the inherited shuffled-TRAINING-label null with a permuted-TEST-"
                         "label null, and declared the scoring definition and sign handling"),
                "why": ("the inherited null is INFORMATIVE AT 1,280 FEATURES AND VACUOUS AT ONE. "
                        "This arm's primary is one feature, so the inherited gate would have "
                        "returned INDETERMINATE_NULL_FAIL after a 4.64 h run for a reason that has "
                        "nothing to do with the heads. The scoring definition was simply absent, "
                        "and the two candidates differ on whether the arm is zero-shot at all."),
                "outcome_neutrality": ("TRIVIALLY SATISFIED AND STATED RATHER THAN ASSUMED: no "
                                       "head has been scored on this arm, so there is no measured "
                                       "verdict that this amendment could move. A bar may be "
                                       "corrected before its first use; after it, this same edit "
                                       "would be a post-hoc change and would not be permissible."),
                "new_threshold_introduced": ("YES -- NULL_SIGNAL_RATIO = %.1f. Declared openly "
                                             "because the usual rule for amending a frozen bar is "
                                             "that no new threshold may be added. That rule governs "
                                             "edits made AFTER results exist; here none do."
                                             % NULL_SIGNAL_RATIO),
                "revert_pre_commitment": ("if the permuted-test-label null is later shown to be "
                                          "insensitive to a fault the shuffled-training-label null "
                                          "would have caught, this amendment is reverted and the "
                                          "arm is re-scored, not re-interpreted."),
            },
        ],
        "honest_limits_fixed_in_advance": [
            ("ZERO-SHOT TRANSFER, not head quality: all five heads were trained on other "
             "libraries, so a null bounds transfer to the GSE144621 tile rung and says nothing "
             "about the heads as predictors of the libraries they were fitted on."),
            ("GEOMETRY MISMATCH: the heads expect 86 positions and the substrate is 150, so every "
             "number here is a windowed application and not a native one. %d windows per tile."
             % probe["geometry"]["n_windows"]),
            ("THE tslpred-kosuri2013 HEAD'S ENDPOINT DOES NOT MATCH THIS RUNG. It predicts "
             "TRANSLATION; promoter activity is a transcription quantity. It is included under the "
             "standing instruction to leave no stones unturned, and its inclusion WIDENS every "
             "head's interval via the alpha/5 correction rather than helping the arm."),
            ("CIRCULARITY IS SCOPED, NOT ABSENT: the urtecho2018 head is circular on the sigma70 "
             "rung (GSE108535) and must never be a comparator there. It is not circular on this "
             "substrate, and that scoping is recorded so nobody carries it to the wrong rung."),
            ("THE INCUMBENT IS A FITTED RIDGE AND THIS ARM IS NOT FITTED AT ALL. The comparison is "
             "therefore asymmetric in tuning budget, and it is asymmetric AGAINST this arm. A loss "
             "bounds the pair, not the heads alone."),
            ("0.941 IS THE REPLICATE CEILING, so the incumbent 0.5076 is 54%% of what is "
             "achievable. An arm that ties the incumbent has not approached the ceiling."),
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    try:
        bar = build_bar(verify())
    except BarRefused as e:
        print("REFUSED: %s" % e, flush=True)
        return 2
    if not bar["reachable_equals_declared"]:
        print("INTERNAL: reachable outcomes != declared outcomes")
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(bar, indent=2, sort_keys=True), encoding="utf-8")

    w = bar["windowing"]
    print("ARM      %s" % bar["arm"])
    print("ROUTE    %s" % w["chosen"])
    print("         %s" % w["why"])
    print("PRIMARY  %s over %d windows, per head; %s declared INELIGIBLE"
          % (PRIMARY_REDUCTION, w["n_windows"], SECONDARY_REDUCTION))
    print("CORRECT  %s -> CI percentiles %s"
          % (bar["multiplicity"]["correction"], bar["multiplicity"]["ci_percentiles"]))
    print("FORECAST %s (p=%.2f)" % (bar["forecast"]["predicted_outcome"],
                                    bar["forecast"]["p_heads_no_gain"]))
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
