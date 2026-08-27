"""PRE-REGISTERED predictions for the v3 staleness prompt, committed BEFORE the fix is written.

WHY PRE-REGISTER. I have adjudicated all 110 items by hand, which makes them a real test set -- and that
is exactly what makes tuning against them dangerous. The v2 memo refused a third benchmark round for this
reason ("tuning a prompt against a 10-item set a third time is fitting"). The 110-item set is bigger but
the trap is identical: run, look, adjust, re-run, and a PASS means "I adjusted until it passed".

The escape is to commit the predictions FIRST, in their own commit, so the fix cannot be quietly reshaped
to match whatever came back. Each prediction below is falsifiable and specific: which flags must flip,
which must NOT, and -- the part that makes it a real test -- which currently-`supported` items must stay
supported. A fix that flips things I did not predict is informative even when the aggregate improves.

THE DIAGNOSED MECHANISM (from 3 independent instances, all the same root cause):
the v2 `ALSO CRITICAL` instruction "implemented code exists -> stale" was written UNSCOPED and appended
AFTER the exceptions, so it silently outranks them. It fires on:
  1. FINDING claims  (bvbrc census: the script is the INSTRUMENT that produced the infeasibility finding,
     so its existence SUPPORTS the claim rather than refuting it)
  2. CORRECTION text (browser.py: the only "deferred" text IS the correction recording that it shipped --
     a case the prompt explicitly names as `supported` and then overrides)
  3. HISTORICAL claims (negative_results_map: "the gates this family needed and did not have" is
     past-tense about what was missing, not refuted by the gates now existing)

Run: uv run python scripts/staleness_v3_preregistration.py
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Prediction:
    artifact: str
    now: str            # the v2 verdict
    predict: str        # what v3 must return
    because: str
    kind: str           # which arm of the fix this exercises


# Committed before `staleness_auditor_v3` exists. Scored by `score_v3.py` afterwards.
PREDICTIONS: tuple[Prediction, ...] = (
    # --- arm 1: FINDING claims must stop being flagged --------------------------------------------
    Prediction(
        "scripts/bvbrc_strict_mic_4drug_census.py", "stale", "supported",
        "'the gent substrate is infeasible' is a FINDING. The census script is the instrument that "
        "produced that finding, so implemented code supports rather than refutes it.",
        "finding"),
    # --- arm 2: CORRECTION text must stop being flagged -------------------------------------------
    Prediction(
        "dna_decode/genome_map/browser.py", "stale", "supported",
        "the only deferred+browser text in CLAUDE.md IS the correction recording that the browser "
        "shipped -- the exception the prompt already states and v2 overrode.",
        "correction"),
    # --- arm 3: HISTORICAL/past-tense claims must stop being flagged -------------------------------
    Prediction(
        "wiki/negative_results_map_2026-06-13.md", "stale", "supported",
        "'the decoder-side gates this family needed and did not have' is past-tense about what was "
        "missing; the gates now existing does not refute it. NOTE: this file also carries the run's "
        "true positive on a DIFFERENT bullet, so the two must be distinguished by claim, not by file.",
        "historical"),

    # --- the part that makes this a real test: what must NOT change --------------------------------
    Prediction(
        "tests/test_models_cache.py", "stale", "stale",
        "a COUNT claim ('8 unit tests', file has 42) is neither capability nor finding. If the fix "
        "silences this, it has over-corrected and traded a false-positive problem for a false-negative "
        "one.",
        "must-hold-TP"),
    Prediction(
        "wiki/decoder_v0_ux_and_success_criterion.md", "stale", "supported",
        "an ordinary misreading (the heading 'Explicit non-criteria for v0' AGREES with the claim). "
        "The fix does not target this, so if it flips, the improvement is not coming from the mechanism "
        "I diagnosed.",
        "not-targeted"),
    Prediction(
        "dna_decode/forward/prosst_scorer.py", "stale", "supported",
        "the model asserted the opposite of the text ('Local, not Kaggle'). Also NOT targeted by the "
        "fix -- listed so an unexplained flip is visible rather than counted as success.",
        "not-targeted"),
)

# The aggregate bar, also pre-registered. v2 scored 11 flags / 2 TP / 9 FP on 110 items.
BAR = {
    "max_flags": 6,          # at least 5 of the 9 false positives must go
    "min_true_positives": 2, # both real catches must survive -- non-negotiable
    "max_false_positives": 4,
}


def targeted() -> tuple[Prediction, ...]:
    return tuple(p for p in PREDICTIONS if p.kind in {"finding", "correction", "historical"})


def must_hold() -> tuple[Prediction, ...]:
    return tuple(p for p in PREDICTIONS if p.kind == "must-hold-TP")


def main() -> int:
    print("PRE-REGISTERED v3 predictions (committed before the fix is written)\n")
    for p in PREDICTIONS:
        print(f"  [{p.kind:14}] {p.artifact}")
        print(f"       {p.now} -> {p.predict}")
    print(f"\naggregate bar: {BAR}")
    print(f"\n{len(targeted())} targeted flips, {len(must_hold())} must-hold true positive(s).")
    print("A fix that clears the aggregate bar while missing the targeted flips has improved by "
          "some other mechanism, and the diagnosis is wrong.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
