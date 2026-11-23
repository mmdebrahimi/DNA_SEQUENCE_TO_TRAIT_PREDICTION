"""Equivalence against a MEASURED margin -- never "free" from a non-significant difference.

Discipline 5 of P2's final measurement-discipline document (2026-11-22), which is the generalised
fix for the one place a peer was told to act on an underpowered null:

    "annotation-anchored re-patching is free" rested on +0.0098 whose CI reached +0.0316 -- a
    possible cost LARGER than the entire arm spread (0.0226) it was measured against.

A confidence interval containing zero means the data could not resolve a difference. It does NOT
mean there is no difference, and it certainly does not mean the difference is too small to matter.
Those are three distinct statements and only the third is "free". This module refuses to let them
be confused: `equivalence()` reports non-significance and equivalence as SEPARATE fields, and a
reader who wants "free" has to supply a margin and watch the interval fit inside it.

WHY A MARGIN CANNOT DEFAULT. A margin is the smallest difference anyone would act on, so it is a
claim about consequences and cannot be derived from the data being tested. Supplying zero turns an
equivalence test back into a test against zero, which is the defect. So `margin` is required, must
be positive, and must arrive with a `source` saying what measured quantity it came from.

MARGINS DO NOT TRANSFER ACROSS SUBSTRATES, AND THIS IS THE ONE PLACE TO SAY SO. P2's
MARGIN_LOOSE = 0.0226 (its arm spread) and MARGIN_STRICT = 0.0075 (its leak family's per-candidate
minimum detectable difference, family-wise 0.0101, superseding an earlier 0.0130 that came from a
pooled-40 file) are measured on E. coli with a logistic-regression discriminator over AUROC. The
human lentiMPRA bar carries an INHERITED_PROHIBITED_COMPARISON forbidding exactly that kind of
cross-substrate transfer: different assay, different noise ceiling, different kingdom. So the
FUNCTION transfers and the NUMBERS do not. A human-arm margin must come from the human arm's own
spread, which is why none is hard-coded here and why `source` is mandatory rather than decorative.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "EquivalenceRefusal",
    "EquivalenceResult",
    "equivalence",
    "P2_REFERENCE_MARGINS",
]

#: Recorded for provenance and for the regression fixture ONLY -- see the module docstring. These
#: are E. coli figures from another substrate and another instrument; importing them as a human-arm
#: margin would be the cross-substrate comparison the frozen bar prohibits.
P2_REFERENCE_MARGINS = {
    "MARGIN_LOOSE": 0.0226,
    "MARGIN_STRICT": 0.0075,
    "family_wise_mdd": 0.0101,
    "superseded": {"value": 0.0130, "why": "came from the superseded pooled-40 file"},
    "substrate": "E. coli tokenization axis; logistic-regression discriminator over AUROC",
    "DO_NOT_TRANSFER": ("the human lentiMPRA bar forbids cross-substrate comparison of this kind; "
                        "derive a human margin from the human arm's own spread"),
}


class EquivalenceRefusal(ValueError):
    """The inputs cannot support an equivalence claim."""


@dataclass(frozen=True)
class EquivalenceResult:
    lo: float
    hi: float
    margin: float
    source: str
    #: True only when the WHOLE interval lies strictly inside (-margin, +margin).
    equivalent: bool
    #: True when the interval contains zero. Independent of `equivalent` -- that is the point.
    contains_zero: bool
    #: The largest absolute effect the interval still permits.
    worst_case_abs: float
    verdict: str

    def as_dict(self) -> dict:
        return {
            "test": "equivalence-against-a-measured-margin-v1",
            "ci": [self.lo, self.hi],
            "margin": self.margin,
            "margin_source": self.source,
            "equivalent": self.equivalent,
            "contains_zero": self.contains_zero,
            "worst_case_abs_effect_still_permitted": self.worst_case_abs,
            "verdict": self.verdict,
            "reading": (
                "`contains_zero` and `equivalent` are SEPARATE. An interval containing zero means "
                "the data could not resolve a difference; only `equivalent` means the difference "
                "is bounded below something that matters. Reporting the first as the second is how "
                "an underpowered null becomes the word 'free'."),
        }


def equivalence(ci, margin: float, *, source: str) -> EquivalenceResult:
    """Is the whole interval inside +/- margin?

    `ci` is a two-element (lo, hi). `margin` is the smallest difference worth acting on, and
    `source` must name the measured quantity it came from.
    """
    try:
        lo, hi = (float(ci[0]), float(ci[1]))
    except (TypeError, IndexError, ValueError) as exc:
        raise EquivalenceRefusal(f"ci must be a two-element (lo, hi); got {ci!r}") from exc
    if math.isnan(lo) or math.isnan(hi) or math.isinf(lo) or math.isinf(hi):
        raise EquivalenceRefusal(f"ci ({lo}, {hi}) is not finite")
    if hi < lo:
        raise EquivalenceRefusal(
            f"ci ({lo}, {hi}) is inverted. A silently swapped interval would widen or narrow the "
            f"comparison below without any other symptom.")
    if not isinstance(source, str) or not source.strip():
        raise EquivalenceRefusal(
            "margin `source` is required and must name the MEASURED quantity the margin came from. "
            "A margin is a claim about consequences, so an unsourced one cannot be checked by a "
            "later reader -- and an unsourced margin is how a convenient one gets chosen.")
    m = float(margin)
    if math.isnan(m):
        raise EquivalenceRefusal("margin is NaN")
    if m <= 0.0:
        raise EquivalenceRefusal(
            f"margin must be > 0; got {m}. A zero margin turns this back into a test against "
            f"ZERO, which is the exact defect this function exists to prevent: no finite interval "
            f"is ever inside +/-0, so every claim would read 'not equivalent' for the wrong reason.")

    worst = max(abs(lo), abs(hi))
    equivalent = (lo > -m) and (hi < m)
    contains_zero = (lo <= 0.0 <= hi)

    if equivalent:
        verdict = "EQUIVALENT_WITHIN_MARGIN"
    elif contains_zero:
        # The dangerous state, and the one that produced the withdrawn claim.
        verdict = "UNDERPOWERED_NOT_EQUIVALENT"
    else:
        verdict = "DIFFERENT_BEYOND_MARGIN"

    return EquivalenceResult(lo=lo, hi=hi, margin=m, source=source.strip(), equivalent=equivalent,
                             contains_zero=contains_zero, worst_case_abs=worst, verdict=verdict)
