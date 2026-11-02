"""Freeze the as-given arm's bar. It is NOT a blind pre-registration, and that is stated, not hidden.

`asgiven-position-resolved-promoteratlas-150bp` was declared 2026-11-02 under R6 because the
positional arm's orientation diagnostic resolved AGAINST orienting: as-given scored **0.5086**
against oriented **0.4893**, so reversing the bin axis on minus-strand tiles destroys about 0.0193
of signal rather than aligning the cohorts. The positional arm's declared co-primary was the
ORIENTED rung, so its eligible primary was handicapped by the transform that same run showed to be
wrong. Changing the DECLARED PRIMARY to the as-given reduction is what makes this a distinct arm
under R6 rather than a re-read.

**THE HAZARD, FIRST, BECAUSE IT IS THE WHOLE WEAKNESS OF THIS ARM.** The number this bar grades is
ALREADY VISIBLE. 0.5086 was measured and published as the positional run's INELIGIBLE secondary
before this arm existed, so no amount of care makes this a blind pre-registration and it must never
be presented as one. Three things follow, all pre-committed here rather than discovered later:

1. **The predicted outcome is a TIE, not a win.** 0.5086 sits above the incumbent 0.5076, but its
   paired CI is [-0.0079, +0.0017] and CONTAINS ZERO. The arm's job is NOT to re-report 0.5086 as a
   gain, and `ASGIVEN_TIES` is named as the expected terminal.
2. **A single primary, and therefore NO multiplicity correction.** One comparison at the plain 95%
   level. Declaring a second rung here would re-open exactly the forking path the positional bar
   closed by making its secondaries ineligible, and would be a far worse offence given the number
   is already known.
3. **What the arm can still establish despite the hazard** is narrow and real: whether the as-given
   reduction, as a DECLARED PRIMARY scored through the full validity chain, clears the rung's
   inherited margin. The answer being foreseeable does not make the gates optional, and the gates
   are what the positional run's secondaries never had to pass as a primary.

Everything else is INHERITED and re-read with a `verify()` that refuses on drift: the margin 0.02,
the committed atlas per-seed comparator, the 10-seed `leave_peak_out(frac=0.25)` split, and
`ridge_spearman` unchanged.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from scripts.freeze_positional_arm_bar import (
    EXPECTED_ATLAS_MEDIAN,
    EXPECTED_GC_MEDIAN,
    EXPECTED_N_SEEDS,
    EXPECTED_N_TILES,
    FAMILY_ALPHA,
    MARGIN,
    BarRefused,
    bonferroni_percentiles,
    ceiling_arms_and_criterion,
    derive_parity_tolerance,
)

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/asgiven_arm_acceptance_bar.json"
FEASIBILITY = REPO / "wiki/positional_arm_feasibility_2026-11-01.json"
POSITIONAL_RESULT = REPO / "wiki/positional_arm_result_2026-11-01.json"
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"

ARM = "asgiven-position-resolved-promoteratlas-150bp"
PRIMARY = "bins10_mean_asgiven"

#: The already-visible number, read from the positional arm's committed artifact rather than typed.
#: Recorded so the hazard is a FIGURE in the bar, not an adjective.
OBSERVED_BEFORE_DECLARED_RUNG = PRIMARY

OUTCOMES = (
    "ASGIVEN_ADDS",
    "ASGIVEN_TIES",
    "ASGIVEN_NO_GAIN",
    "INDETERMINATE_SPLIT_MISMATCH",
    "INDETERMINATE_PARITY_FAIL",
    "INDETERMINATE_NULL_FAIL",
)

CHECK_ORDER = (
    ("split_equality", "INDETERMINATE_SPLIT_MISMATCH"),
    ("parity", "INDETERMINATE_PARITY_FAIL"),
    ("shuffled_label_null", "INDETERMINATE_NULL_FAIL"),
    ("primary_decision", ("ASGIVEN_ADDS", "ASGIVEN_TIES", "ASGIVEN_NO_GAIN")),
)


def reachable_outcomes(check_order=CHECK_ORDER) -> set[str]:
    out: set[str] = set()
    for _gate, result in check_order:
        out.update((result,) if isinstance(result, str) else result)
    return out


def classify(ci_lo: float, ci_hi: float, *, margin: float = MARGIN) -> tuple[str, str]:
    """Pure, ordered, first match wins. Exhaustive over the real line by construction.

    `ASGIVEN_TIES` deliberately absorbs an interval lying wholly inside `(0, margin]`. That is a
    demonstrated gain which does NOT clear the margin the rung grades against, and calling it an
    ADD would be grading this arm against a softer bar than its siblings. Pre-committed, so it
    cannot be re-read after a number lands between 0 and 0.02.
    """
    if ci_lo > margin:
        return ("ASGIVEN_ADDS",
                "the 95%% interval lies entirely above the inherited margin %+.4f" % margin)
    if ci_hi < 0.0:
        return ("ASGIVEN_NO_GAIN",
                "the interval lies entirely BELOW zero, so the as-given reduction is strictly "
                "worse than the mean-pool incumbent")
    return ("ASGIVEN_TIES",
            "the interval neither clears the margin %+.4f nor lies below zero, so no gain of at "
            "least the margin is demonstrated and no loss is either" % margin)


def verify(*, reconcile_path: Path = RECONCILE, ledger_path: Path = LEDGER,
           positional_result: Path = POSITIONAL_RESULT,
           feasibility_path: Path = FEASIBILITY) -> dict:
    """Re-read every anchor live. Raises `BarRefused` on any drift."""
    for p in (reconcile_path, ledger_path, positional_result, feasibility_path):
        if not Path(p).is_file():
            raise BarRefused("anchor source absent: %s" % p)

    rec = json.loads(Path(reconcile_path).read_text(encoding="utf-8"))
    if rec.get("verdict") != "RECONCILED" or rec.get("mismatched"):
        raise BarRefused("the comparator is %r with mismatched %r"
                         % (rec.get("verdict"), rec.get("mismatched")))
    if rec.get("atlas_median") != EXPECTED_ATLAS_MEDIAN:
        raise BarRefused("atlas median drifted: %r vs %r"
                         % (rec.get("atlas_median"), EXPECTED_ATLAS_MEDIAN))
    if rec.get("gc_median") != EXPECTED_GC_MEDIAN:
        raise BarRefused("GC median drifted: %r vs %r" % (rec.get("gc_median"),
                                                          EXPECTED_GC_MEDIAN))
    if int(rec.get("n_seeds", -1)) != EXPECTED_N_SEEDS:
        raise BarRefused("comparator has %r seeds" % rec.get("n_seeds"))
    if int(rec.get("n_tiles", -1)) != EXPECTED_N_TILES:
        raise BarRefused("comparator has %r tiles" % rec.get("n_tiles"))

    ledger = json.loads(Path(ledger_path).read_text(encoding="utf-8"))
    arms, criterion = ceiling_arms_and_criterion(ledger)
    if ARM not in arms:
        raise BarRefused(
            "the arm %r is NOT declared on %r; a bar for an undeclared arm would grade work the "
            "ledger does not know exists" % (ARM, criterion)
        )
    if arms[ARM].get("status") != "planned":
        raise BarRefused(
            "the arm's status is %r, not 'planned'. A bar may only be frozen while NO number "
            "exists for the arm." % arms[ARM].get("status")
        )
    if arms[ARM].get("metric") is not None:
        raise BarRefused("the arm already carries metric %r" % arms[ARM].get("metric"))

    # the positional arm's committed artifact is where the already-visible number comes from
    pos = json.loads(Path(positional_result).read_text(encoding="utf-8"))
    rung = pos.get("rungs", {}).get(PRIMARY)
    if rung is None or rung.get("skipped"):
        raise BarRefused(
            "the positional arm's artifact has no scored %r rung, so the observed-before-declared "
            "hazard cannot be stated with a figure and this bar must not be frozen without it"
            % PRIMARY
        )
    if rung.get("eligible_as_the_arms_result") is not False:
        raise BarRefused(
            "%r is not recorded as INELIGIBLE in the positional artifact; the premise of this arm "
            "is that it was a secondary there and becomes a PRIMARY here" % PRIMARY
        )

    feas = json.loads(Path(feasibility_path).read_text(encoding="utf-8"))
    if feas.get("verdict") != "FEASIBLE" or not feas["parity"]["passed"]:
        raise BarRefused("feasibility verdict %r / parity passed %r"
                         % (feas.get("verdict"), feas["parity"]["passed"]))

    extra = reachable_outcomes() ^ set(OUTCOMES)
    if extra:
        raise BarRefused(
            "the reachable outcome set and the closed set differ by %r; containment is not enough, "
            "because an unreachable outcome is a dead validity branch" % sorted(extra)
        )
    return {"reconcile": rec, "arms": arms, "criterion": criterion, "positional": pos,
            "feasibility": feas}


def build_bar(v: dict) -> dict:
    rec, pos, feas, criterion = v["reconcile"], v["positional"], v["feasibility"], v["criterion"]
    rung = pos["rungs"][PRIMARY]
    oriented = pos["rungs"]["bins10_mean_oriented"]
    lo, hi = bonferroni_percentiles(1)          # one comparison -> the plain 95% interval
    parity_tol = derive_parity_tolerance(float(feas["parity"]["max_abs_deviation"]))

    return {
        "schema": "asgiven-arm-acceptance-bar-v1",
        "frozen": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "frozen_before_any_score_for_THIS_arm_exists": True,
        "arm": ARM,
        "criterion": criterion,
        "question": (
            "Does the AS-GIVEN position-resolved reduction, as a DECLARED PRIMARY scored through "
            "the full validity chain, beat the mean-pool incumbent (median Spearman %.4f) by the "
            "rung's inherited margin of %.2f?" % (rec["atlas_median"], MARGIN)
        ),
        "why_this_is_a_distinct_arm_under_R6": (
            "It changes the DECLARED PRIMARY. In the positional arm `bins10_mean_asgiven` was a "
            "SECONDARY and INELIGIBLE as that arm's result; here it is the primary and the only "
            "rung. The positional arm's verdict is untouched and is never retro-edited."
        ),
        "NOT_A_BLIND_PRE_REGISTRATION": {
            "statement": (
                "The number this bar grades is ALREADY VISIBLE and was published before this arm "
                "existed. This is not a blind pre-registration and must never be presented as one."
            ),
            "already_observed_rung": OBSERVED_BEFORE_DECLARED_RUNG,
            "already_observed_median": rung["median_spearman"],
            "already_observed_gain": rung["gain_bootstrap"]["median"],
            "already_observed_ci": [rung["gain_bootstrap"]["ci_lo"], rung["gain_bootstrap"]["ci_hi"]],
            "already_observed_ci_level": rung["gain_bootstrap"]["ci_percentiles"],
            "source": "wiki/positional_arm_result_2026-11-01.json (read, not retyped)",
            "what_it_costs": (
                "A pre-registration's value is that the bar could not have been shaped by the "
                "result. That guarantee is absent here, so this arm's evidential weight is LOWER "
                "than the positional arm's, whatever it returns."
            ),
            "what_survives": (
                "Whether the as-given reduction clears the inherited margin as a PRIMARY, scored "
                "through split-equality, parity and a shuffled-label null -- gates the positional "
                "run's secondaries never had to pass in that role. A foreseeable answer does not "
                "make the gates optional."
            ),
        },
        "inherited_not_chosen": {
            "margin": {"value": MARGIN,
                       "derivation": ("INHERITED verbatim from the tile-rung arms, as in the "
                                      "positional bar, so the bars stay commensurable")},
            "comparator": {
                "source": "data/promoteratlas/local_reconcile_2026-11-01.json",
                "verdict": rec["verdict"],
                "atlas_median": rec["atlas_median"],
                "atlas_per_seed": rec["atlas_per_seed"],
                "must_be_recomputed_in_run": False,
                "why": "the incumbent is never re-run, so this arm cannot move its own comparator",
            },
            "split": {
                "spec": "genomewide.leave_peak_out(tiles, frac=0.25, seed=s) for s in 0..9",
                "frac": rec["frac"], "n_seeds": rec["n_seeds"],
                "n_test_per_seed": rec["n_test_per_seed"],
                "leakage_unit": "the called expression PEAK, not the tile",
            },
            "scorer": {"spec": "dna_decode.glm.atlas_embed.ridge_spearman, unchanged"},
            "features": {
                "source": "the COMMITTED stored block; no second forward pass",
                "block": "D:/dna_decode_cache/promoteratlas_positional/positional_block.npy",
                "reduction": ("bins10_mean sliced from the stored block and flattened WITHOUT "
                              "`orient` -- which is the entire difference from the positional "
                              "arm's primary"),
            },
        },
        "multiplicity": {
            "n_primaries": 1,
            "primary": PRIMARY,
            "family_alpha": FAMILY_ALPHA,
            "correction": "NONE -- one comparison, judged at the plain 95% level",
            "ci_percentiles": [lo, hi],
            "no_secondaries_declared": True,
            "why_no_secondaries": (
                "A second rung here would re-open the forking path the positional bar closed by "
                "making its secondaries ineligible, and would be a worse offence given the "
                "primary's number is already known. One rung, one interval, one verdict."
            ),
        },
        "validity_controls": {
            "split_equality": {
                "spec": "recompute the split for seeds 0-9 and assert n_test_per_seed matches",
                "on_failure": "INDETERMINATE_SPLIT_MISMATCH, and stop",
            },
            "parity": {
                "spec": ("derive the global mean from the stored 10-bin means and assert it "
                         "reproduces the committed embedding and its per-seed Spearmans"),
                "tolerance": parity_tol,
                "on_failure": "INDETERMINATE_PARITY_FAIL, and stop",
                "why": "it is what makes this arm's number comparable to %.4f" % rec["atlas_median"],
            },
            "shuffled_label_null": {
                "spec": ("refit the ridge on shuffled TRAIN labels, same splits and seeds, scoring "
                         "against REAL test labels -- the incumbent's IMPLEMENTATION, not its prose"),
                "bar": "|median| < %.2f" % MARGIN,
                "on_failure": "INDETERMINATE_NULL_FAIL, and SPEND MORE SEEDS -- never re-specify",
            },
        },
        "outcomes_defined_before_running": {
            "ASGIVEN_ADDS": "the 95%% interval lies entirely above +%.2f" % MARGIN,
            "ASGIVEN_TIES": ("the interval neither clears +%.2f nor lies below zero. This "
                             "DELIBERATELY absorbs an interval inside (0, %.2f]: a sub-margin gain "
                             "is not a gain by the bar this rung grades against, and calling it an "
                             "ADD would grade this arm against a softer bar than its siblings"
                             % (MARGIN, MARGIN)),
            "ASGIVEN_NO_GAIN": "the interval lies entirely below zero",
            "INDETERMINATE_SPLIT_MISMATCH": "the recomputed splits differ from the committed ones",
            "INDETERMINATE_PARITY_FAIL": "the derived global mean does not reproduce the committed "
                                         "embedding",
            "INDETERMINATE_NULL_FAIL": "the shuffled-label control cleared the bar; nothing graded",
        },
        "verdict_check_order": [g for g, _ in CHECK_ORDER],
        "check_order_rationale": "Validity FIRST, split equality before everything else.",
        "closed_outcome_set": list(OUTCOMES),
        "reachable_outcome_set": sorted(reachable_outcomes()),
        "reachable_equals_closed": sorted(reachable_outcomes()) == sorted(OUTCOMES),
        "pre_commitments": [
            "The predicted outcome is ASGIVEN_TIES. The arm's job is NOT to re-report %.4f as a "
            "gain." % rung["median_spearman"],
            "The margin %.2f is INHERITED and will not be re-specified after a number exists."
            % MARGIN,
            "A sub-margin gain returns ASGIVEN_TIES, decided here and not after seeing where the "
            "interval lands.",
            "No secondary rung will be added to this arm. If another reduction is worth scoring it "
            "is a NEW arm under R6.",
            "The positional arm's 0.4893, its verdict and its bar are never retro-edited.",
            "If the null fails, the response is MORE SEEDS, never a weaker control.",
        ],
        "honest_limits_fixed_in_advance": [
            "NOT a blind pre-registration -- see NOT_A_BLIND_PRE_REGISTRATION. This arm's "
            "evidential weight is lower than the positional arm's whatever it returns.",
            "The orientation convention remains a 1-BIT INFERENCE from a performance difference "
            "between two rungs, not established biology. This arm tests the better-performing side "
            "of that inference; it does not establish which orientation `seq` is given in.",
            "No promoter anchor exists on this substrate, so even an ADDS verdict would say that "
            "position bins relative to an ARBITRARY 150 bp window help -- not that a promoter "
            "register was recovered.",
            "The encoder is FROZEN; this bounds nothing about a fine-tuned one.",
            "689 duplicate sequences are unchanged and inherited by BOTH sides of the paired "
            "comparison.",
            "The rung ceiling is 0.941 and replicate-derived. No tile-rung number may be read as a "
            "fraction of it without saying so.",
        ],
        "forecast": {
            "recorded_before_any_score_for_this_arm": True,
            "predicted_outcome": "ASGIVEN_TIES",
            "p_asgiven_ties": 0.80,
            "p_asgiven_adds": 0.05,
            "p_asgiven_no_gain": 0.15,
            "rationale": (
                "The interval is already known to contain zero at the corrected level the "
                "positional run used, and widening to a plain 95%% interval makes containing zero "
                "MORE likely, not less. The 0.15 on NO_GAIN is not symmetry: the positional run's "
                "value came from a secondary scored at alpha/6, and re-scoring at 95%% could in "
                "principle resolve a small negative, though the point estimate %+.4f argues "
                "against it." % rung["gain_bootstrap"]["median"]
            ),
            "comparison_point": {
                "oriented_primary_median": oriented["median_spearman"],
                "asgiven_median": rung["median_spearman"],
                "difference": round(rung["median_spearman"] - oriented["median_spearman"], 4),
            },
            "my_forecast_record_in_this_project": (
                "4 of 8. The positional arm's bar put p=0.45 on POSITIONAL_NO_GAIN and it landed "
                "there, which is the fourth hit. DO NOT defer to these priors -- and note this "
                "particular forecast is nearly free, because the number is already visible."
            ),
        },
        "reproduce": [
            "uv run python scripts/freeze_asgiven_arm_bar.py",
            "uv run python -m scripts.asgiven_arm_score",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--print-only", action="store_true")
    args = ap.parse_args(argv)

    try:
        v = verify()
    except BarRefused as e:
        print("REFUSED TO FREEZE: %s" % e, flush=True)
        return 2

    bar = build_bar(v)
    if args.print_only:
        print(json.dumps(bar, indent=2, sort_keys=True))
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(bar, indent=2, sort_keys=True), encoding="utf-8")

    h = bar["NOT_A_BLIND_PRE_REGISTRATION"]
    print("arm          %s (status planned, no metric)" % ARM)
    print("primary      %s -- ONE comparison, plain 95%% interval %r"
          % (PRIMARY, bar["multiplicity"]["ci_percentiles"]))
    print("margin       %+.2f INHERITED" % MARGIN)
    print("HAZARD       the number is ALREADY VISIBLE: %s = %.4f, gain %+.4f, CI %r"
          % (h["already_observed_rung"], h["already_observed_median"], h["already_observed_gain"],
             h["already_observed_ci"]))
    print("             -> NOT a blind pre-registration, stated in the bar")
    print("parity tol   %.0e" % bar["validity_controls"]["parity"]["tolerance"])
    print("outcomes     %d closed, reachable EQUALS closed: %s"
          % (len(OUTCOMES), bar["reachable_equals_closed"]))
    print("forecast     %s (p=%.2f) -- record 4 of 8, and this forecast is nearly free"
          % (bar["forecast"]["predicted_outcome"], bar["forecast"]["p_asgiven_ties"]))
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
