"""Grade the corpus-selection-exposure experiment against its frozen bar. Pure, ordered, testable.

**Validity first, and every outcome reachable.** Run 1 of the parent PromoterAtlas arm emitted its
STRONGEST verdict while the shuffled-label control was failing, because that gate's clause read
`abs(med(null) - med(gc)) < MARGIN` and could therefore never fire -- the most consequential line in
a validity-first order was dead code. The fix is two things together: the indeterminates are
evaluated BEFORE the substantive outcomes, and `tests/test_exposure_verdict.py` constructs an input
that reaches each outcome, so no branch can be dead.

The order is the frozen bar's `verdict_check_order` and is not re-derived here.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Mirrors wiki/promoteratlas_exposure_acceptance_bar.json verdict_check_order. First match wins.
CHECK_ORDER: tuple[str, ...] = (
    "INDETERMINATE_HARNESS_DRIFT",
    "INDETERMINATE_CONTROLS_UNSEPARABLE",
    "INDETERMINATE_REMOVAL_COST_UNBOUNDED",
    "EXPOSURE_EXPLAINS",
    "EXPOSURE_PARTIAL",
    "EXPOSURE_EXCLUDED_LINEAR",
)

#: The cheap-feature separability floor below which the control cohort, not the encoder, is the
#: thing being measured. An AUROC of 0.5 is chance; 0.60 is the frozen bar's floor.
CHEAP_SEPARABILITY_FLOOR = 0.60


@dataclass(frozen=True)
class ExposureInput:
    """Everything the grader reads. All medians over the committed seeds."""

    harness_reconciled: bool
    sep_auroc_gc_only: float
    atlas_baseline: float
    atlas_minus_separation_dir: float
    atlas_minus_random_dir: float
    gc_baseline: float
    margin: float = 0.02

    @property
    def drop(self) -> float:
        """How much removing the control/tile separation direction costs the atlas channel."""
        return self.atlas_baseline - self.atlas_minus_separation_dir

    @property
    def drop_random(self) -> float:
        """How much removing an ARBITRARY direction costs. Bounds the removal's own price."""
        return self.atlas_baseline - self.atlas_minus_random_dir

    @property
    def surviving_margin_over_gc(self) -> float:
        """What the encoder still buys over the dumb baseline after the detection axis is gone."""
        return self.atlas_minus_separation_dir - self.gc_baseline


@dataclass(frozen=True)
class ExposureVerdict:
    verdict: str
    detail: str
    is_valid: bool


def grade(x: ExposureInput) -> ExposureVerdict:
    """Return the first matching outcome in the frozen check order."""
    if not x.harness_reconciled:
        return ExposureVerdict(
            "INDETERMINATE_HARNESS_DRIFT",
            "the local harness did not reproduce the committed figures; no exposure claim may rest "
            "on it. Diagnose the harness, do not widen the tolerance.",
            False,
        )

    if x.sep_auroc_gc_only < CHEAP_SEPARABILITY_FLOOR:
        return ExposureVerdict(
            "INDETERMINATE_CONTROLS_UNSEPARABLE",
            "GC-only control-vs-tile AUROC %.4f < %.2f, so a cheap feature cannot separate these "
            "controls either. 'The encoder cannot separate them' would be a statement about the "
            "CONTROLS. Pre-committed response: MORE controls, never a weaker control."
            % (x.sep_auroc_gc_only, CHEAP_SEPARABILITY_FLOOR),
            False,
        )

    if x.drop_random > x.margin:
        return ExposureVerdict(
            "INDETERMINATE_REMOVAL_COST_UNBOUNDED",
            "removing a RANDOM direction already costs %.4f > margin %.2f, so a drop cannot be "
            "attributed to the separation direction specifically -- the cell would be measuring the "
            "price of projection itself." % (x.drop_random, x.margin),
            False,
        )

    if x.surviving_margin_over_gc < x.margin:
        return ExposureVerdict(
            "EXPOSURE_EXPLAINS",
            "after removing the detection direction the encoder retains only %.4f over GC, below "
            "the inherited margin %.2f. The headline was substantially promoter DETECTION. The "
            "recorded metric keeps its measured value; the honest limit becomes a quantified "
            "discount." % (x.surviving_margin_over_gc, x.margin),
            True,
        )

    if x.drop > x.margin:
        return ExposureVerdict(
            "EXPOSURE_PARTIAL",
            "removing the detection direction costs %.4f (> margin %.2f) but %.4f still survives "
            "over GC. A material, quantified discount rather than a refutation."
            % (x.drop, x.margin, x.surviving_margin_over_gc),
            True,
        )

    return ExposureVerdict(
        "EXPOSURE_EXCLUDED_LINEAR",
        "removing the detection direction costs only %.4f (<= margin %.2f), so the gain survives "
        "LINEAR removal essentially intact. BOUNDED BY CONSTRUCTION: a nonlinear detection signal "
        "survives a linear projection, so this may NEVER be written up as 'exposure excluded'."
        % (x.drop, x.margin),
        True,
    )
