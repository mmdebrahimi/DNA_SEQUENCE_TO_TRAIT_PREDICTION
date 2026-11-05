"""Grade the frame-conditioning experiment against its frozen bar. Pure, ordered, testable.

**Validity first, and every outcome reachable.** Run 1 of the parent PromoterAtlas arm emitted its
STRONGEST verdict while a control was failing, because that gate's clause could never fire -- the
most consequential line in a validity-first order was dead code. The fix is two things together:
every indeterminate is evaluated BEFORE the substantive outcomes, and `tests/
test_frame_verdict.py` constructs an input reaching each outcome and asserts the reachable set
EQUALS `CHECK_ORDER` rather than merely being contained in it.

Mirrors `dna_decode/glm/exposure_verdict.py` in SHAPE rather than importing it: the inputs differ,
and the shared discipline is the ordering plus the reachability proof.

**DEVIATION FROM THE PLAN, recorded rather than quietly taken.** The technical plan names SEVEN
outcomes; this grader has EIGHT. The addition is `INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE`, and
it exists because the plan's own pre-execution review (that plan folder's `brainstorm.md`, under
critical issues) found that *"shrinkage is not absorption unless it exceeds what a
matched-informativeness, non-positional feature produces -- the subtractive experiment had a
random-direction control for exactly this; the additive one has no analogue"*, and that the bar
"cannot be frozen without" that control and the size-matched random stratification. A control that
cannot change the verdict is not a control: without this outcome the grader would be FORCED to
report `FRAME_EXPLAINS` over a measurement it knows is confounded with "the baseline simply got
better", which is the identical defect the ordering above exists to prevent. The direct ancestor is
`INDETERMINATE_REMOVAL_COST_UNBOUNDED` in the subtractive grader. Whether the plan's outcome set
should be amended to eight is recorded as a question for the user, not settled here.

**No threshold here is invented.** The only free constants are the inherited margin (0.02, from the
parent arm's bar) and the locate floor (0.95, which the measured 0.9918 clears with real headroom
rather than being fitted to it). Informativeness is tested against the MARGIN itself, on the
principle that a feature carrying less signal than the gap under test cannot explain that gap --
so no new number is conjured for it.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Mirrors wiki/promoteratlas_frame_acceptance_bar.json verdict_check_order. First match wins.
CHECK_ORDER: tuple[str, ...] = (
    "INDETERMINATE_HARNESS_DRIFT",
    "INDETERMINATE_LOCATE_RATE",
    # The plumbing fault is diagnosed BEFORE the frame's informativeness, and the order is
    # load-bearing rather than cosmetic: a baseline that never received the feature can MAKE the
    # frame look uninformative, so reporting "uninformative" first would attribute a harness fault
    # to the biology. The technical plan specifies this order ("a non-strengthened baseline beats
    # both"); an earlier draft of this module had the two reversed and its own ordering test caught
    # it.
    "INDETERMINATE_BASELINE_NOT_STRENGTHENED",
    "INDETERMINATE_FRAME_UNINFORMATIVE",
    "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE",
    "FRAME_EXPLAINS",
    "FRAME_PARTIAL",
    "FRAME_EXCLUDED_CONDITIONAL",
)

#: Below this unique-locate rate the cohort being measured is not the cohort that was scored.
#: Measured 0.9918, so the floor has headroom and is not fitted to the observation.
LOCATE_FLOOR = 0.95


@dataclass(frozen=True)
class FrameInput:
    """Everything the grader reads. Every figure a MEDIAN over the committed seeds.

    `gc_plus_matched` is the strongest NON-POSITIONAL strengthened baseline (the measured series is
    `gc + dinucleotide` and `gc + kmer4`); it is what makes absorption attributable to the frame
    rather than to the baseline having more columns. `random_stratum_min_margin` is the
    size-matched random stratification null for the C3 cell, which is what separates a genuinely
    collapsed within-stratum margin from the range restriction that stratifying causes by itself.
    """

    harness_reconciled: bool
    unique_locate_rate: float

    frame_alone_spearman: float          # the frame features' own expression signal
    train_fit_gc: float                  # training-set fit of the plain GC baseline
    train_fit_gc_frame: float            # ... and of the frame-strengthened one

    atlas: float
    gc_baseline: float
    gc_plus_frame: float
    gc_plus_matched: float               # strongest non-positional matched-strength baseline

    atlas_in_window: float
    gc_in_window: float
    atlas_out_window: float
    gc_out_window: float
    random_stratum_min_margin: float     # size-matched RANDOM stratification null

    margin: float = 0.02

    # ---- C2: the strengthened-baseline cell ------------------------------------------------
    @property
    def margin_over_gc(self) -> float:
        """What the encoder buys over the unstrengthened baseline."""
        return self.atlas - self.gc_baseline

    @property
    def margin_over_gc_frame(self) -> float:
        """What it still buys once the baseline is handed the corpus frame."""
        return self.atlas - self.gc_plus_frame

    @property
    def margin_over_gc_matched(self) -> float:
        """... and once handed an equally strong feature that carries no position."""
        return self.atlas - self.gc_plus_matched

    @property
    def frame_absorption(self) -> float:
        """How much of the margin the frame took."""
        return self.margin_over_gc - self.margin_over_gc_frame

    @property
    def matched_absorption(self) -> float:
        """How much ANY equally strong feature takes. The confound's own size."""
        return self.margin_over_gc - self.margin_over_gc_matched

    @property
    def attributable_absorption(self) -> float:
        """Absorption in EXCESS of the matched control. The frame-specific part."""
        return self.frame_absorption - self.matched_absorption

    # ---- C3: the stratified cell -----------------------------------------------------------
    @property
    def margin_in_window(self) -> float:
        return self.atlas_in_window - self.gc_in_window

    @property
    def margin_out_window(self) -> float:
        return self.atlas_out_window - self.gc_out_window

    @property
    def min_stratum_margin(self) -> float:
        """The weaker of the two strata. A collapse in EITHER is the cell firing."""
        return min(self.margin_in_window, self.margin_out_window)

    # ---- the two arms, each with its own attribution test ----------------------------------
    @property
    def c2_absorbs(self) -> bool:
        return self.margin_over_gc_frame < self.margin

    @property
    def c2_attributable(self) -> bool:
        """The frame must absorb strictly MORE than a matched non-positional feature."""
        return self.attributable_absorption > 0.0

    @property
    def c3_absorbs(self) -> bool:
        return self.min_stratum_margin < self.margin

    @property
    def c3_attributable(self) -> bool:
        """Stratifying at all must NOT already push a margin below the bar."""
        return self.random_stratum_min_margin >= self.margin

    @property
    def explains(self) -> bool:
        return (self.c2_absorbs and self.c2_attributable) or (
            self.c3_absorbs and self.c3_attributable
        )

    @property
    def absorbs_unattributably(self) -> bool:
        """A cell fired, but every cell that fired is confounded."""
        return (self.c2_absorbs or self.c3_absorbs) and not self.explains


@dataclass(frozen=True)
class FrameVerdict:
    verdict: str
    detail: str
    is_valid: bool


def grade(x: FrameInput) -> FrameVerdict:
    """Return the first matching outcome in the frozen check order."""
    if not x.harness_reconciled:
        return FrameVerdict(
            "INDETERMINATE_HARNESS_DRIFT",
            "the cached embedding pass no longer reproduces the committed arms.native150, so no "
            "frame claim may rest on it. Diagnose the harness; do not widen the tolerance.",
            False,
        )

    if x.unique_locate_rate < LOCATE_FLOOR:
        return FrameVerdict(
            "INDETERMINATE_LOCATE_RATE",
            "only %.4f of distinct tile sequences locate UNIQUELY in the genome (floor %.2f), so "
            "the frame is unknown for too much of the cohort and the scored subset is no longer "
            "the cohort the arm measured. Ambiguous and not-found tiles are excluded and counted, "
            "never resolved by taking the first hit."
            % (x.unique_locate_rate, LOCATE_FLOOR),
            False,
        )

    if x.train_fit_gc_frame < x.train_fit_gc:
        return FrameVerdict(
            "INDETERMINATE_BASELINE_NOT_STRENGTHENED",
            "the frame-strengthened baseline fits the TRAINING data worse (%.4f) than plain GC "
            "(%.4f). Adding a column cannot lower a least-squares training fit, so this is a "
            "plumbing fault, not a result: the baseline never received what it was given, and "
            "'no absorption' would be a statement about the harness. Checked BEFORE "
            "informativeness, because a feed this broken can also make the frame look "
            "uninformative." % (x.train_fit_gc_frame, x.train_fit_gc),
            False,
        )

    if abs(x.frame_alone_spearman) < x.margin:
        return FrameVerdict(
            "INDETERMINATE_FRAME_UNINFORMATIVE",
            "the frame features carry only %.4f of expression signal on their own, below the "
            "margin %.2f under test. 'Adding them does not absorb the gain' would be VACUOUSLY "
            "true -- a feature carrying less signal than the gap cannot explain the gap. "
            "Pre-committed response: a better positional feature, never a weaker test."
            % (x.frame_alone_spearman, x.margin),
            False,
        )

    if x.absorbs_unattributably:
        return FrameVerdict(
            "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE",
            "a cell fired but the absorption is not attributable to POSITION. The frame took "
            "%.4f of the margin while an equally strong NON-positional feature took %.4f "
            "(frame-specific excess %+.4f), and the size-matched random stratification null "
            "already leaves a minimum stratum margin of %.4f against the bar %.2f. So the "
            "shrinkage is confounded with 'the baseline simply got better' and with the range "
            "restriction stratifying causes by itself. Pre-committed response: a sharper "
            "positional feature or a better-matched control, never dropping the control."
            % (x.frame_absorption, x.matched_absorption, x.attributable_absorption,
               x.random_stratum_min_margin, x.margin),
            False,
        )

    if x.explains:
        return FrameVerdict(
            "FRAME_EXPLAINS",
            "corpus-window position carries the gain. Over the frame-strengthened baseline the "
            "encoder retains %.4f against the inherited margin %.2f, and the weaker stratum "
            "retains %.4f; the absorption exceeds the matched non-positional control by %+.4f, so "
            "it is attributable to POSITION rather than to baseline strength. arms.native150 "
            "keeps its measured value; the honest limit becomes a quantified positional discount."
            % (x.margin_over_gc_frame, x.margin, x.min_stratum_margin,
               x.attributable_absorption, ),
            True,
        )

    if x.frame_absorption > x.margin:
        return FrameVerdict(
            "FRAME_PARTIAL",
            "a material, quantified positional discount rather than a refutation. Handing the "
            "frame to the baseline absorbs %.4f of the margin (> %.2f) but %.4f still survives, "
            "and the weaker stratum retains %.4f. Of that absorption %+.4f is frame-specific, the "
            "rest being what any equally strong feature takes."
            % (x.frame_absorption, x.margin, x.margin_over_gc_frame,
               x.min_stratum_margin, x.attributable_absorption),
            True,
        )

    return FrameVerdict(
        "FRAME_EXCLUDED_CONDITIONAL",
        "the gain survives both the frame-strengthened baseline (%.4f retained) and BOTH strata "
        "(weaker stratum %.4f), absorbing only %.4f (<= margin %.2f). BOUNDED BY CONSTRUCTION: "
        "this excludes THE CORPUS-WINDOW RULE AS REPLICATED -- 200 bp upstream of an annotated "
        "gene start with a 5 bp intergenic gate, under sorted ordering, on ONE organism's "
        "annotation. A sharper or differently-granular positional feature is not excluded, and "
        "this may NEVER be written up as 'position excluded'."
        % (x.margin_over_gc_frame, x.min_stratum_margin, x.frame_absorption, x.margin),
        True,
    )
