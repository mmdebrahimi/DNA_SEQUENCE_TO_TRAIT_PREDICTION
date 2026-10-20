"""Pre-run contract for the GLM sampling knobs.

WHY THIS EXISTS (2026-10-20). A temperature sweep ran three grid points, 21 minutes of
compute, and produced three IDENTICAL numbers -- because `temperature` never reached the
sampler. The invariance was then reported as a FINDING ("temperature does not rescue the
generator") before anyone checked that the swept variable varied at all. One comparison at
n=1 would have caught it.

The lesson is not "be careful". A sweep's entire validity rests on the swept variable
actually moving the output, and that is a CHECKABLE PRECONDITION. This module makes it one:

  * `validate_sampling` runs at CONSTRUCTION, before the ~80 s model load, so a bad value
    fails in milliseconds instead of deep inside generation.
  * `assert_knob_wired` refuses to let a run proceed -- and so refuses to let an artifact
    STAMP a non-default sampling value -- when the requested knob produced no warper. That
    is the exact regression that made three artifacts claim `is_recorded_default: false`
    over runs that all used default sampling.

Neither is enforcement of honesty; both are cheap preconditions on a run that would
otherwise be silently void.
"""
from __future__ import annotations

# Derived, not asserted (2026-10-20). `top_k` on this model truncates the 4096-entry 6-mer
# distribution, NOT the 4 bases -- the remote `_BPLogitsProcessor` marginalises k-mer logits
# to per-base probabilities AFTER any warper runs. Truncating too hard makes a base
# UNREACHABLE at a position (its marginal becomes exactly 0), i.e. silently greedy there.
#
# Measured over 15 conditions (5 seeds x logit sharpness 1.0/2.0/4.0), counting (position,
# base) cells driven to zero out of 24:
#     top_k=4  -> 7-9 cells zero   (this is the value someone types meaning "4 bases")
#     top_k=10 -> up to 5
#     top_k=14..21 -> up to 1 (flat across this whole band)
#     top_k>=22 -> 0 in every condition tested, and at every larger value
# So 22 is the smallest value with no measured reachability loss.
#
# IT WAS FIRST SET TO 24, AND THAT WAS WRONG -- the original scan went 16, 24, 32, so 17-23
# were never measured and the constant was PADDED by 2 on the assumption that a rounder
# number was safer. `test_the_floor_is_the_SMALLEST_value_with_no_reachability_loss` caught
# it by re-deriving the boundary, which is the only reason it is 22 now. Padding a threshold
# is still asserting it; the test refuses both a floor that is too low AND one that is too
# high, so the constant cannot drift back to a remembered number.
#
# CAVEAT: measured on synthetic Gaussian logits over 5 seeds x sharpness 1.0/2.0/4.0, not on
# the real model's logits -- indicative, not exhaustive. The real-model end-to-end check is
# the `--verify-knob` diff.
TOP_K_KMER_FLOOR = 22


class SamplingContractError(ValueError):
    """A sampling setting that would be inert, crash late, or silently distort."""


def validate_sampling(temperature: float | None, top_k: int | None) -> None:
    """Refuse sampling settings that cannot mean what the caller intends.

    Called from `HFGenerator.__init__` so it fires BEFORE the model load. Previously
    `--temperature 0` was silently ignored; once the knob was wired it became a ValueError
    raised ~80 s in, after the load. Both are bad; this is the fix for both.
    """
    if temperature is not None:
        if temperature != temperature:  # NaN
            raise SamplingContractError("temperature is NaN")
        if temperature <= 0:
            raise SamplingContractError(
                f"temperature must be > 0, got {temperature}. Temperature divides the "
                f"logits, so 0 is undefined and negatives invert the distribution. For "
                f"greedy decoding do not use this sampler."
            )

    if top_k:
        if top_k < 0:
            raise SamplingContractError(f"top_k must be >= 0 (0 disables), got {top_k}")
        if top_k < TOP_K_KMER_FLOOR:
            raise SamplingContractError(
                f"top_k={top_k} is below the measured floor of {TOP_K_KMER_FLOOR}.\n"
                f"On this model top_k truncates the 4096 SIX-MERS, not the 4 bases -- the "
                f"per-base distribution is marginalised from k-mer logits afterwards. "
                f"Measured: below {TOP_K_KMER_FLOOR}, at least one base becomes UNREACHABLE "
                f"at some position (marginal exactly 0), making generation silently greedy "
                f"there; top_k=4 zeroes 7-9 of 24 (position, base) cells.\n"
                f"If you meant 'sample from the 4 bases', that is what temperature alone "
                f"already does. Use 0 to disable, or >= {TOP_K_KMER_FLOOR}."
            )


def requested_non_default(temperature: float | None, top_k: int | None) -> bool:
    """True when the caller asked for sampling that differs from the recorded default."""
    return bool((temperature is not None and temperature != 1.0) or top_k)


def assert_knob_wired(temperature: float | None, top_k: int | None, processors) -> None:
    """Refuse a run whose requested knob produced no warper.

    This is the guard the temperature sweep lacked. A non-default knob that installs nothing
    is INERT, so the run would measure default sampling while its artifact stamped the
    requested value -- a false provenance claim, which is exactly what happened across three
    sweep artifacts. Cheap: a list-length check, no model needed.
    """
    if requested_non_default(temperature, top_k) and not list(processors):
        raise SamplingContractError(
            f"sampling was requested (temperature={temperature}, top_k={top_k}) but NO "
            f"warper was installed, so the setting is INERT and any artifact stamping it "
            f"would be false. This is the 2026-10-20 regression: passing `temperature=` to "
            f"generate() lands AFTER the remote one-hot collapse and does nothing. Warpers "
            f"must go through `logits_processor=`."
        )


def verify_knob_live(generate_at, *, low: float = 0.7, high: float = 1.5) -> dict:
    """End-to-end proof that the knob moves the OUTPUT. The check the sweep lacked.

    `generate_at(temperature) -> list[str]` must generate with a FIXED seed, so the only
    difference between the two calls is temperature. (`HFGenerator.generate` reseeds to
    `self.seed` on entry, which is what makes one loaded model a valid control.)

    Returns a stampable record. `live` is True only when the two settings produced different
    sequences -- identical output under a doubled temperature is the dead-knob signature,
    and it is the ONLY thing that distinguishes "temperature does not help this generator"
    from "temperature was never applied".

    Deliberately NOT cheap to fake: it compares generated strings, not the processor list.
    The processor-list check (`assert_knob_wired`) catches our half of the plumbing; this
    catches the half that lives in third-party code and in framework ordering.
    """
    a = list(generate_at(low))
    b = list(generate_at(high))
    n = min(len(a), len(b))
    identical = sum(1 for x, y in zip(a, b) if x == y)
    if n == 0:
        return {"live": False, "verdict": "INDETERMINATE_NO_SEQUENCES",
                "n": 0, "identical": 0, "low": low, "high": high}
    return {
        "live": identical < n,
        "verdict": ("TEMPERATURE_IS_IGNORED" if identical == n
                    else "TEMPERATURE_IS_LIVE" if identical == 0
                    else f"PARTIAL_{identical}_of_{n}"),
        "n": n, "identical": identical, "low": low, "high": high,
        "note": "same seed both calls; only temperature differs",
    }


# --------------------------------------------------------------------------- #
# Comparison-mode contract (Step 1 of the sampling-sweep plan, 2026-10-20)
# --------------------------------------------------------------------------- #
# WHY: a sweep was run with `--natural-set promoter` and NO `--conditional`, which
# scores an UNCONDITIONAL sample (drawn from the whole-genome prior, mostly coding)
# against ALIGNED promoter windows, while the Markov null is fitted on the promoter
# FIT half. The null holds the target distribution and the model does not, so the
# comparison is close to tautological. The repo documented this in two places
# BEFORE the run and nothing enforced it:
#   corpus.py:168-175  "An unconditional sample cannot be scored against an aligned
#                       set fairly, so the meaningful task is *given the genomic
#                       context, continue into the promoter*."
#   glm_generator_falsifier.py --conditional help: "required for a meaningful
#                       comparison against the ALIGNED promoter set"
VALID_COMPARISON = ("promoter", True)


class ComparisonModeError(ValueError):
    """A (natural_set, conditional) pairing that cannot produce a meaningful number."""


def validate_comparison_mode(natural_set: str, conditional: bool) -> None:
    """Refuse every (reference, conditioning) pairing but the one that is meaningful.

    Called before any corpus build or model load, so a mis-specified run costs
    milliseconds instead of ~30 minutes of generation plus a number nobody can use.
    """
    if natural_set == "both":
        raise ComparisonModeError(
            "natural_set='both' is refused: the runner builds conditional prompts from "
            "`naturals[kinds[0]]['eval_windows']` only, so the SECOND set's candidates "
            "would be prompted from the FIRST set's windows -- a silent cross-reference "
            "mismatch. Run the two references as separate invocations."
        )

    if natural_set == "promoter":
        if not conditional:
            raise ComparisonModeError(
                "natural_set='promoter' without --conditional is refused.\n"
                "An UNCONDITIONAL sample is drawn from the whole-genome prior (mostly "
                "coding) while the promoter reference is ALIGNED to the CDS start and the "
                "Markov null is fitted on the promoter FIT half -- so the null holds the "
                "target distribution and the model does not, making the comparison close "
                "to tautological.\n"
                "This is not a new judgement: corpus.py:168-175 says an unconditional "
                "sample cannot be scored against an aligned set fairly, and the "
                "--conditional help calls it REQUIRED for this reference. A sweep run this "
                "way on 2026-10-19 produced distinguishability 0.8141 against a null of "
                "0.3876 and the gap was almost entirely the conditioning mismatch.\n"
                "Pass --conditional."
            )
        return  # the one valid pairing

    if natural_set == "genomic":
        if conditional:
            raise ComparisonModeError(
                "natural_set='genomic' with --conditional is refused: a prompted "
                "continuation is positionally anchored, and the genomic reference is "
                "UNALIGNED uniform windows, so there is no shared offset structure to "
                "score against."
            )
        raise ComparisonModeError(
            "natural_set='genomic' without --conditional is refused on MEASURED grounds, "
            "not on preference: corpus.py:171-173 records that against unaligned random "
            "genomic windows an order-3 Markov null scores 0.5086 against a 0.5255 "
            "ceiling -- BELOW its own noise floor -- so that reference has no "
            "discriminating power at all. A run against it cannot distinguish a good "
            "generator from a bad one."
        )

    raise ComparisonModeError(
        f"unknown natural_set {natural_set!r}; expected 'promoter', 'genomic' or 'both'"
    )
