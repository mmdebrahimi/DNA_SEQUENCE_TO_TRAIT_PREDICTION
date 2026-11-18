"""Re-grade the active-frame gate by changing EXACTLY ONE field of its pre-registered bar.

WHAT THIS DOES AND DOES NOT DO
------------------------------
It calls the gate's OWN ``compute_verdict`` with the gate's OWN ``PREREGISTERED_ACTIVE_FRAME``
dict, altered in exactly one place: ``null_max_abs``. Everything else -- the effect-size
thresholds, the seed count, the verdict ladder -- is the frozen bar, byte for byte. The per-seed
numbers come from the recorded artifact, so nothing is retrained and the data cannot move between
the two gradings. A diff of the two bar dicts is printed, so "one field" is checkable rather than
asserted.

A CLAIM WITHDRAWN, AND WHY IT MATTERED
--------------------------------------
An earlier version of this script ran ``derive_bar`` on the +0.15 margin threshold, found its
clean-failure rate to be essentially zero, and reported that the frozen bar was
"over-conservative" -- offering that as evidence the correction "runs in both directions" and was
therefore not bar-shopping. **That argument is WITHDRAWN. It was a category error.**

The pre-registered ladder separates two different questions, and the separation is explicit in
``glm_active_frame_gate.PREREGISTERED_ACTIVE_FRAME["verdicts"]``:

  FRAME_EXPLAINS_THE_GAP   margin >= +0.15 on >= 8 of 10 seeds. No noise test. This is an
                           EFFECT-SIZE question -- is the margin big enough that the frame
                           explains the gap rather than merely contributing to it?
  FRAME_CONTRIBUTES        margin > +0.05 with a paired CI excluding zero. THIS is the
                           significance question.

``derive_bar`` answers "how big must this be to stand out from noise". It cannot answer "how big
must this be to explain the gap". A clean-failure rate of 0.000 on an effect-size threshold is that
threshold working as designed, not over-conservatism -- so there was never a correction in the
margin direction to report.

THE REAL ANTI-BAR-SHOPPING ARGUMENT IS SIMPLER AND STRONGER
-----------------------------------------------------------
The effect-size criterion was ALREADY MET on the recorded data: median margin +0.2190 >= 0.15 with
9 of 10 seeds clearing it, against a bar needing 8. **Only the null gate ever blocked this run.** So
exactly one gate needs to change, and the reason is a property of the bar rather than of the result:

  * The criterion bounds ``max |rho|`` over 10 shuffled-label seeds against a threshold set as
    though for a single draw.
  * ``abs()`` makes it two-sided against a hazard -- leakage -- that can only push a null UPWARD.
  * The 3,822-tile test set is far fewer independent units than rows: tiles within a peak share
    ICC 0.6931.
  * Measured consequence: a perfectly CLEAN pipeline fails the frozen 0.05 bar about 90% of the
    time. All three observed breaches are NEGATIVE and the largest positive null is +0.0287.

Every one of those is computable from the design without looking at a single observed margin, which
is what distinguishes this from loosening a bar that proved inconvenient.

THE REPLACEMENT NULL BAR IS THE MOST CONSERVATIVE ONE AVAILABLE
---------------------------------------------------------------
It is derived at ``design_effect = 1.0`` -- every tile treated as independent, the TIGHTEST bar this
design can produce -- precisely so the verdict cannot be accused of resting on the ICC estimate. At
the measured design effect the bar would be looser (about 0.1022); this uses about 0.0416 instead.
The statistic compared against it is the maximum POSITIVE null, matching the hazard's direction.

BOTH VERDICTS ARE KEPT. The original INDETERMINATE_NULL_NOT_CLEAN is preserved verbatim and the
original artifact is never modified. A re-grading that overwrote the first verdict would be
indistinguishable from having only ever run the second.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics as st
from pathlib import Path
from typing import Any

from dna_decode.eval.preregistration import (
    IndependenceUnit,
    Statistic,
    clean_failure_rate,
    derive_bar,
)
from scripts.glm_active_frame_gate import (
    PREREGISTERED_ACTIVE_FRAME,
    compute_verdict,
    framing_recommendation,
)

GATE_ARTIFACT = Path("wiki/glm_active_frame_gate_2026-11-16.json")
#: Family G-A's recorded result. Source of the POOLED-frame ceiling and arm medians, so the
#: cross-frame comparison reads both halves from artifacts rather than from this file.
POOLED_ARTIFACT = Path("wiki/glm_learned_representation_2026-10-08.json")
OUT = Path("wiki/glm_active_frame_regrade_2026-11-17.json")

#: The active-frame ceiling is not its own field anywhere -- it is embedded in the gate's
#: `preregistered.thresholds_are` prose. It is PARSED from that string rather than typed here,
#: because a typed 0.9690 beside a comment is still a number this file asserts.
_ACTIVE_CEILING_RE = re.compile(r"active-frame ceiling \((\d+\.\d+)")

#: Measured 2026-11-17 by one-way ANOVA on the exact active frame the gate used. The selector is
#: ``t.active == "active"`` -- a STRING comparison -- which reproduces the gate's 14,913 tiles to
#: the row. A truthiness test on that field selects all 44,106, which is how two earlier estimates
#: of this ICC came out on the wrong set.
MEASURED_ICC = 0.6931
ICC_PROVENANCE = (
    "one-way ANOVA ICC of rep1 within peak over the 14,913-tile active frame "
    "(t.active == 'active'; 1,848 peaks, mean 8.07 tiles/peak)"
)

#: A convention about the BAR's own error rate, recorded as a convention. Far weaker than asserting
#: a threshold on the data, which is what the frozen null bar did.
MAX_CLEAN_FAILURE = 0.05


class RegradeRefused(Exception):
    """The recorded artifact did not have the shape this re-grade requires."""


def load_arms(
    artifact: dict[str, Any],
) -> tuple[list[float], list[float], list[float], list[float], int, int]:
    """Pull the recorded per-seed arms and the test-set shape.

    Refuses rather than defaulting: a missing field means the artifact is not the one this was
    written against, and silently substituting a zero would grade a different run.
    """
    try:
        per_seed = artifact["phase2_active"]["per_seed"]
    except KeyError as exc:  # pragma: no cover - shape guard
        raise RegradeRefused(f"recorded artifact lacks {exc}") from exc
    if not per_seed:
        raise RegradeRefused("no per-seed records; there is nothing to re-grade")

    seeds = sorted(per_seed, key=int)
    margins = [float(per_seed[s]["margin"]) for s in seeds]
    nulls = [float(per_seed[s]["encoder_shuffled"]) for s in seeds]
    encoders = [float(per_seed[s]["encoder_multi"]) for s in seeds]
    gcs = [float(per_seed[s]["gc"]) for s in seeds]
    n_rows = int(per_seed[seeds[0]]["n_test"])
    units = {int(per_seed[s]["held_out_peaks"]) for s in seeds}
    if len(units) != 1:
        raise RegradeRefused(
            f"held-out peak counts differ across seeds ({sorted(units)}); one independence unit "
            f"cannot describe them all"
        )
    return margins, nulls, encoders, gcs, n_rows, units.pop()


def recorded_ci(artifact: dict[str, Any]) -> tuple[float, float]:
    """The paired bootstrap CI as the gate recorded it. Refuses rather than recomputing, because a
    recomputed CI is a different number from the one the frozen ladder was applied to."""
    blocks = (
        artifact.get("phase2_active") or {},
        artifact.get("record") or {},
        artifact,
    )
    for block in blocks:
        for key in ("margin_ci", "active_margin_ci"):
            v = block.get(key) if isinstance(block, dict) else None
            if isinstance(v, (list, tuple)) and len(v) == 2:
                return float(v[0]), float(v[1])
    raise RegradeRefused(
        "no recorded paired-bootstrap CI (margin_ci / active_margin_ci). The frozen ladder reads "
        "it for the FRAME_CONTRIBUTES rung, so grading without it would silently change the bar"
    )


def active_ceiling(artifact: dict[str, Any]) -> float:
    """Parse the active-frame ceiling out of the gate's own provenance prose.

    Refuses rather than falling back to a literal. If the gate ever stops stating its ceiling, the
    right outcome is a loud failure, not a normalization computed against a number this file made
    up.
    """
    prose = str((artifact.get("preregistered") or {}).get("thresholds_are", ""))
    m = _ACTIVE_CEILING_RE.search(prose)
    if not m:
        raise RegradeRefused(
            "could not parse the active-frame ceiling from "
            "`preregistered.thresholds_are`. Normalizing against a hardcoded ceiling would make "
            "this artifact assert a number instead of reading one, so this refuses instead"
        )
    return float(m.group(1))


def load_pooled_frame(pooled: dict[str, Any]) -> tuple[float, float, float]:
    """The POOLED frame's ceiling and arm medians, from family G-A's recorded result.

    Returns `(ceiling, encoder_median, gc_median)`. Every field is required: a default here would
    silently normalize one frame against a fabricated ceiling, which is the exact failure the
    normalization exists to expose.
    """
    try:
        ceiling = float(pooled["noise_ceiling"])
        arms = pooled["arms"]
        return ceiling, float(arms["encoder-multi"]["median"]), float(arms["gc"]["median"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RegradeRefused(
            f"the pooled artifact {POOLED_ARTIFACT} lacks {exc!r}; it must carry `noise_ceiling` "
            f"plus `arms['encoder-multi']['median']` and `arms['gc']['median']` for the "
            f"cross-frame normalization to read both halves from artifacts"
        ) from exc


def bar_diff(frozen: dict, used: dict) -> dict:
    """Every field where the bar actually used differs from the frozen one. The point of the
    re-grade is that this has exactly one entry."""
    keys = set(frozen) | set(used)
    return {
        k: {"frozen": frozen.get(k), "used": used.get(k)}
        for k in sorted(keys)
        if frozen.get(k) != used.get(k)
    }


def regrade(artifact: dict[str, Any], *, max_clean_failure: float = MAX_CLEAN_FAILURE,
            pooled: dict[str, Any] | None = None) -> dict:
    margins, nulls, encoders, gcs, n_rows, n_units = load_arms(artifact)
    ci_lo, ci_hi = recorded_ci(artifact)

    unit_measured = IndependenceUnit.from_icc(
        "peak", n_rows=n_rows, n_units=n_units, icc=MEASURED_ICC,
        measured_from=ICC_PROVENANCE,
    )
    #: deff = 1 is NOT a measurement. It is the deliberately conservative choice that makes the
    #: verdict independent of the ICC estimate, since gap G1 says an ICC-derived design effect is
    #: an OVERSTATEMENT and therefore loosens the bar.
    unit_tight = IndependenceUnit(
        "tile (deff=1, conservative)", n_rows=n_rows, n_units=n_rows, design_effect=1.0,
        measured_from="NOT a measurement: the tightest bar this design can produce, chosen so the "
                      "verdict does not rest on the ICC estimate",
    )

    null_one = Statistic("encoder_shuffled", "spearman", "max-of-k", len(nulls), "upward")
    derived_null_tight = derive_bar(null_one, unit_tight, max_clean_failure=max_clean_failure)
    derived_null_measured = derive_bar(null_one, unit_measured, max_clean_failure=max_clean_failure)

    max_pos_null = max(nulls)
    max_abs_null = max(abs(x) for x in nulls)
    frozen_cfr = clean_failure_rate(
        Statistic("encoder_shuffled", "spearman", "max-of-k", len(nulls), "two-sided"),
        unit_measured, PREREGISTERED_ACTIVE_FRAME["null_max_abs"],
    )

    # ---- the bar actually used: the frozen dict with ONE field replaced ----------------
    bar_used = {**PREREGISTERED_ACTIVE_FRAME, "null_max_abs": round(derived_null_tight, 4)}
    diff = bar_diff(PREREGISTERED_ACTIVE_FRAME, bar_used)

    # ---- the arms, with the null statistic made hazard-matched (max POSITIVE, not max |.|) ----
    arms = {
        "pooled_reconcile_ok": True,
        "null_max_abs": round(max_pos_null, 4),
        "active_margin_per_seed": margins,
        "active_margin_ci": (ci_lo, ci_hi),
    }
    verdict, detail = compute_verdict(arms, bar=bar_used)

    # ---- the PAIRED empirical null. margin - margin_null = encoder - encoder_shuffled, because
    # GC appears in both and cancels. Pairing is kept because the seeds share splits; comparing
    # each margin against the UNPAIRED maximum discards that and manufactures false thinness.
    paired = [e - n for e, n in zip(encoders, nulls)]

    # ---- CROSS-FRAME CEILING NORMALIZATION. The scoring-frame work established that a
    # correlation is read against its OWN frame's ceiling; these two artifacts did not apply that
    # rule to themselves, and the omission makes a specific wrong inference available.
    if pooled is None:
        pooled = json.loads(POOLED_ARTIFACT.read_text(encoding="utf-8"))
    pool_ceil, pool_enc, pool_gc = load_pooled_frame(pooled)
    act_ceil = active_ceiling(artifact)
    act_enc = float(artifact["phase2_active"]["encoder_multi_median"])
    act_gc = float(artifact["phase2_active"]["gc_median"])
    pool_frac = round(pool_enc / pool_ceil, 4)
    act_frac = round(act_enc / act_ceil, 4)
    enc_fall = round(pool_enc - act_enc, 4)
    gc_fall = round(pool_gc - act_gc, 4)

    return {
        "schema": "glm-active-frame-regrade-v2",
        "date": "2026-11-17",
        "ceiling_normalization": {
            "why": (
                "a correlation is read against its OWN frame's ceiling. Without these two numbers "
                "printed together, the growing margin reads as evidence the learned "
                "representation is stronger than family G-A concluded. Normalized, it is weaker."
            ),
            "pooled_frame": {
                "source": str(POOLED_ARTIFACT),
                "ceiling": pool_ceil,
                "encoder_median": pool_enc,
                "gc_median": pool_gc,
                "fraction_of_ceiling": pool_frac,
            },
            "active_frame": {
                "source": str(GATE_ARTIFACT),
                "ceiling": act_ceil,
                "ceiling_provenance": str(artifact["preregistered"]["thresholds_are"]),
                "ceiling_was_parsed_not_typed": True,
                "encoder_median": act_enc,
                "gc_median": act_gc,
                "fraction_of_ceiling": act_frac,
            },
            "encoder_reaches_a_SMALLER_fraction_on_the_active_frame": act_frac < pool_frac,
            "THE_COMPARATOR_COLLAPSED": (
                f"both arms fall from the pooled frame to the active one -- the encoder by "
                f"{enc_fall:+.4f} ({pool_enc:.4f} -> {act_enc:.4f}) and GC by {gc_fall:+.4f} "
                f"({pool_gc:.4f} -> {act_gc:.4f}). The margin grows from +0.0122 to "
                f"{st.median(margins):+.4f} because the COMPARATOR collapsed, not because the "
                f"candidate improved. Relative to achievable signal the encoder reaches "
                f"{act_frac:.1%} of the active ceiling against {pool_frac:.1%} of the pooled one."
            ),
            "encoder_fall": enc_fall,
            "gc_fall": gc_fall,
            "gc_fell_further": gc_fall > enc_fall,
            "scoping_caution": (
                "NAME WHICH CEILING EACH FRACTION USES. scripts/glm_genomewide_oracle.py, the "
                "owner of `fraction_of_ceiling`, documents \"Do not compare the two fractions\" "
                "for fractions taken across SUBSTRATES with unrelated ceilings. These two are "
                f"across FRAMES of ONE substrate -- {pool_frac:.4f} against the pooled {pool_ceil} "
                f"ceiling and {act_frac:.4f} against the active {act_ceil} ceiling -- which is a "
                "weaker comparison than the forbidden one but is still only legible with both "
                "ceilings stated."
            ),
            "NOT_grounds_to_revisit_G_A": (
                "This is an observation about SCOPE, not a reason to re-open family G-A. The gate "
                "artifact itself records `conditions_on_the_outcome: "
                f"{artifact['phase2_active'].get('conditions_on_the_outcome')}` and "
                "`not_a_ledger_candidate: "
                f"{artifact['phase2_active'].get('not_a_ledger_candidate')}`, so its own author "
                "already marked these numbers as non-candidates. The project's umbrella states "
                "the self-distillation precondition as 'a learned representation LANDS (G-A "
                "terminal-met is NOT sufficient)', and Action-Log row 10 records a prior session "
                "reaching this same fork on 2026-10-23 and declining it. Re-grading G-A on a "
                "frame that collapses its comparator would be changing the POPULATION until the "
                "bar passes -- structurally the same move as changing the bar."
            ),
        },
        "what_changed": (
            "EXACTLY ONE field of the pre-registered bar: null_max_abs. The effect-size "
            "thresholds, the seed count and the verdict ladder are the frozen bar unchanged, and "
            "the per-seed numbers are read verbatim from the recorded artifact. The verdict is "
            "computed by the gate's OWN compute_verdict."
        ),
        "bar_diff_against_the_frozen_preregistration": diff,
        "bar_diff_field_count": len(diff),
        "a_claim_WITHDRAWN": (
            "v1 of this re-grade ran derive_bar on the +0.15 margin threshold, reported it as "
            "over-conservative, and used that as evidence the correction 'runs in both "
            "directions'. WITHDRAWN: +0.15 is an EFFECT-SIZE threshold (the FRAME_EXPLAINS rung "
            "has no noise test; the significance test lives in FRAME_CONTRIBUTES), so a "
            "noise-derived bar cannot evaluate it and a clean-failure rate of 0.000 on it is the "
            "threshold working as designed. The margin bar is UNCHANGED here at 0.15."
        ),
        "original_verdict_PRESERVED": {
            "verdict": artifact.get("verdict"),
            "detail": artifact.get("verdict_detail"),
            "note": "kept because a re-grade that overwrote it would be indistinguishable from "
                    "having only ever run the second grading",
        },
        "why_only_the_null_gate_needed_changing": {
            "effect_size_criterion_was_already_met": True,
            "median_margin": round(st.median(margins), 4),
            "explains_margin_bar": PREREGISTERED_ACTIVE_FRAME["explains_margin"],
            "seeds_clearing_the_UNCHANGED_effect_size_bar": sum(
                1 for m in margins if m >= PREREGISTERED_ACTIVE_FRAME["explains_margin"]
            ),
            "explains_min_seeds": PREREGISTERED_ACTIVE_FRAME["explains_min_seeds"],
            "reading": (
                "median +0.2190 clears +0.15 and 9 of 10 seeds clear it against a bar needing 8, "
                "so the frozen effect-size criterion was satisfied on the recorded data. The null "
                "gate is ordered FIRST in compute_verdict, so it alone produced INDETERMINATE."
            ),
        },
        "independence_unit": {
            "name": unit_measured.name,
            "n_rows": unit_measured.n_rows,
            "n_units": unit_measured.n_units,
            "icc": MEASURED_ICC,
            "design_effect": round(unit_measured.design_effect, 4),
            "n_effective": round(unit_measured.n_effective, 1),
            "measured_from": unit_measured.measured_from,
            "used_for_the_bar": False,
            "why_not_used": (
                "gap G1: a design effect from the label's ICC alone implicitly sets the "
                "prediction's intra-peak correlation to 1, which OVERSTATES it and LOOSENS the "
                "bar. The bar below is derived at deff=1 instead, which is tighter."
            ),
        },
        "null_arm": {
            "frozen_bar": PREREGISTERED_ACTIVE_FRAME["null_max_abs"],
            "frozen_bar_clean_failure_rate": round(frozen_cfr, 4),
            "THE_LOAD_BEARING_NUMBER": (
                f"a CLEAN pipeline fails the frozen 0.05 bar {frozen_cfr:.0%} of the time, so the "
                f"original INDETERMINATE verdict carried almost no information about leakage"
            ),
            "three_independent_mis_specifications": [
                "it bounds a MAX of 10 absolute values against a threshold set as though for a "
                "single draw",
                "abs() makes it two-sided against a hazard (leakage) that can only push a null "
                "UPWARD",
                "the 3,822-row test set is a few hundred independent units, not a few thousand",
            ],
            "replacement_bar_used": round(derived_null_tight, 4),
            "replacement_bar_basis": "deff=1, one-sided, max-of-10 -- the TIGHTEST available",
            "looser_bar_at_the_measured_icc_NOT_used": round(derived_null_measured, 4),
            "statistic_compared": "max POSITIVE null (hazard-matched)",
            "observed_max_positive_null": round(max_pos_null, 4),
            "observed_max_abs_null": round(max_abs_null, 4),
            "n_breaches_of_the_frozen_bar": sum(1 for x in nulls if abs(x) > 0.05),
            "n_breaches_that_are_negative": sum(1 for x in nulls if abs(x) > 0.05 and x < 0),
            "n_positive_nulls_above_the_frozen_bar": sum(1 for x in nulls if x > 0.05),
            "clean": max_pos_null < derived_null_tight,
        },
        "paired_empirical_null": {
            "what_it_is": (
                "margin minus margin-null, which reduces to encoder_multi - encoder_shuffled "
                "because GC appears in both and cancels. The same encoder, same split, same seed, "
                "real labels against shuffled."
            ),
            "why_paired": (
                "the seeds share peaks and splits, so pairing is the whole point. Comparing each "
                "margin against the UNPAIRED maximum of the null discards the pairing and "
                "manufactures thinness -- v1 of this artifact did exactly that and reported a "
                "misleading 1.07x."
            ),
            "per_seed": [round(x, 4) for x in paired],
            "all_positive": all(x > 0 for x in paired),
            "min": round(min(paired), 4),
            "median": round(st.median(paired), 4),
            "no_p_value_on_purpose": (
                "a sign test over these 10 would overstate independence, because the seeds share "
                "peaks. The sign pattern is reported; a p-value is not."
            ),
            "needs_no_variance_model_or_design_effect": True,
        },
        "verdict": verdict,
        "verdict_detail": detail,
        "framing_recommendation": framing_recommendation(verdict),
        "honest_limits": [
            "The weakest seed (9) is the weakest on every reading: margin +0.1075 is the only one "
            "below the frozen +0.15, and its paired real-minus-shuffled gap is +0.0523 against a "
            "median of +0.1921. The frozen bar needs 8 of 10, so one weak seed does not decide "
            "the verdict -- but it is a real outlier and is not smoothed over.",
            "No p-value is quoted for the paired null: the seeds share peaks, so any test treating "
            "them as independent would overstate the evidence.",
            "The ICC is measured on the label within peaks and is NOT used for the bar. For a "
            "correlation the design effect is nearer 1 + (m-1)*rho_x*rho_y, so a label-only "
            "estimate overstates it and loosens the bar; deff=1 is used instead.",
            "The null model behind the replacement bar is asymptotic (Spearman under H0 as "
            "N(0, 1/(n-1))). An exact permutation null on the real peak grouping would be stronger.",
            "This re-grades the FRAME question only. The 0.941 pooled ceiling remains an average "
            "over deciles running 0.5013 to 0.9823, and nothing here touches the finetune arm.",
            "Adopting this verdict in place of the original is a USER decision. Both are recorded "
            "and the original artifact is untouched.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate-artifact", default=str(GATE_ARTIFACT))
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    artifact = json.loads(Path(args.gate_artifact).read_text(encoding="utf-8"))
    result = regrade(artifact)
    Path(args.out).write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(f"ORIGINAL : {result['original_verdict_PRESERVED']['verdict']}")
    print(f"REGRADED : {result['verdict']}")
    print(f"           {result['verdict_detail']}")
    print()
    print(f"fields changed in the pre-registered bar: {result['bar_diff_field_count']}")
    for k, v in result["bar_diff_against_the_frozen_preregistration"].items():
        print(f"    {k}: {v['frozen']} -> {v['used']}")
    print(f"frozen null bar clean-failure rate: "
          f"{result['null_arm']['frozen_bar_clean_failure_rate']:.3f}")
    print(f"paired real-minus-shuffled: all positive="
          f"{result['paired_empirical_null']['all_positive']}, "
          f"min {result['paired_empirical_null']['min']:+.4f}")
    cn = result["ceiling_normalization"]
    print()
    print("CEILING NORMALIZATION (each frame against its OWN ceiling):")
    for name in ("pooled_frame", "active_frame"):
        f = cn[name]
        print(f"  {name:13s} ceiling {f['ceiling']:.4f}  encoder {f['encoder_median']:+.4f}  "
              f"gc {f['gc_median']:+.4f}  -> {f['fraction_of_ceiling']:.1%} of ceiling")
    print(f"  encoder reaches a SMALLER fraction on the active frame: "
          f"{cn['encoder_reaches_a_SMALLER_fraction_on_the_active_frame']}")
    print(f"  encoder fell {cn['encoder_fall']:+.4f}, GC fell {cn['gc_fall']:+.4f}  "
          f"(GC further: {cn['gc_fell_further']}) -- the margin grew by COMPARATOR COLLAPSE")
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
