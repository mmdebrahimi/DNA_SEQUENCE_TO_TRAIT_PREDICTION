"""Tests for the frame-conditioning grader.

**The reachability test is the load-bearing one.** Run 1 of the parent arm emitted its strongest
verdict over a gate whose predicate could never fire, so the test asserts the reachable set EQUALS
`CHECK_ORDER` rather than merely being contained in it -- containment would permit exactly that
defect. The per-outcome inputs are held as DATA in `CASES` so a future edit cannot silently cover
fewer outcomes than the order declares.
"""

from __future__ import annotations

import dataclasses

import pytest

from dna_decode.glm.frame_verdict import (
    CHECK_ORDER,
    LOCATE_FLOOR,
    FrameInput,
    grade,
)

#: A clean, valid, gain-survives input. Every other case is this one perturbed, so each test names
#: exactly one difference and nothing rides along.
BASE = FrameInput(
    harness_reconciled=True,
    unique_locate_rate=0.9918,
    frame_alone_spearman=0.10,
    train_fit_gc=0.30,
    train_fit_gc_frame=0.33,
    atlas=0.5076,
    gc_baseline=0.3112,
    gc_plus_frame=0.3200,
    gc_plus_matched=0.3180,
    atlas_in_window=0.50,
    gc_in_window=0.31,
    atlas_out_window=0.49,
    gc_out_window=0.30,
    random_stratum_min_margin=0.18,
)


def _with(**kw) -> FrameInput:
    return dataclasses.replace(BASE, **kw)


#: One input per outcome, as DATA. Keyed by the verdict each must produce.
CASES: dict[str, FrameInput] = {
    # --- the five indeterminates ---
    "INDETERMINATE_HARNESS_DRIFT": _with(harness_reconciled=False),
    "INDETERMINATE_LOCATE_RATE": _with(unique_locate_rate=0.80),
    "INDETERMINATE_FRAME_UNINFORMATIVE": _with(frame_alone_spearman=0.005),
    "INDETERMINATE_BASELINE_NOT_STRENGTHENED": _with(train_fit_gc_frame=0.29),
    # the frame absorbs the margin, but an equally strong non-positional feature absorbs just as
    # much -- so the shrinkage is "the baseline got better", not position
    "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE": _with(
        gc_plus_frame=0.50, gc_plus_matched=0.50
    ),
    # --- the three substantive outcomes ---
    "FRAME_EXPLAINS": _with(gc_plus_frame=0.50),
    "FRAME_PARTIAL": _with(gc_plus_frame=0.40),
    "FRAME_EXCLUDED_CONDITIONAL": BASE,
}


# ------------------------------------------------------------------------------- reachability


def test_every_outcome_in_the_check_order_is_reachable():
    for want, x in CASES.items():
        assert grade(x).verdict == want, "case for %s graded %s" % (want, grade(x).verdict)


def test_the_reachable_set_EQUALS_the_check_order_not_merely_contained_in_it():
    """Containment would permit a dead branch, which is the defect this grader exists to avoid."""
    reachable = {grade(x).verdict for x in CASES.values()}
    assert reachable == set(CHECK_ORDER), (
        "unreachable outcomes: %s / graded-but-undeclared: %s"
        % (sorted(set(CHECK_ORDER) - reachable), sorted(reachable - set(CHECK_ORDER)))
    )


def test_the_cases_table_covers_the_check_order_exactly():
    assert set(CASES) == set(CHECK_ORDER)
    assert len(CHECK_ORDER) == len(set(CHECK_ORDER)), "a duplicated outcome name"


def test_only_the_substantive_outcomes_are_valid():
    for name, x in CASES.items():
        v = grade(x)
        assert v.is_valid == (not name.startswith("INDETERMINATE_")), name


def test_every_detail_string_is_substantial_and_names_no_placeholder():
    for name, x in CASES.items():
        d = grade(x).detail
        assert len(d) > 80, name
        assert "%" not in d, "%s left an unformatted conversion in its detail" % name


# ------------------------------------------------------------------------------------ ordering


def test_harness_drift_beats_every_substantive_outcome():
    for name in ("FRAME_EXPLAINS", "FRAME_PARTIAL", "FRAME_EXCLUDED_CONDITIONAL"):
        x = dataclasses.replace(CASES[name], harness_reconciled=False)
        assert grade(x).verdict == "INDETERMINATE_HARNESS_DRIFT", name


def test_an_uninformative_frame_beats_frame_explains():
    """'Adding them did not absorb the gain' is vacuous when they carry no signal."""
    x = dataclasses.replace(CASES["FRAME_EXPLAINS"], frame_alone_spearman=0.001)
    assert grade(x).verdict == "INDETERMINATE_FRAME_UNINFORMATIVE"


def test_a_non_strengthened_baseline_beats_both_of_those():
    x = dataclasses.replace(
        CASES["FRAME_EXPLAINS"], frame_alone_spearman=0.001, train_fit_gc_frame=0.10
    )
    assert grade(x).verdict == "INDETERMINATE_BASELINE_NOT_STRENGTHENED"


def test_a_bad_locate_rate_beats_everything_except_harness_drift():
    x = dataclasses.replace(CASES["FRAME_EXPLAINS"], unique_locate_rate=0.5)
    assert grade(x).verdict == "INDETERMINATE_LOCATE_RATE"
    x2 = dataclasses.replace(x, harness_reconciled=False)
    assert grade(x2).verdict == "INDETERMINATE_HARNESS_DRIFT"


def test_non_attributable_absorption_beats_frame_explains():
    """The control must be able to DEMOTE the strongest verdict, or it is not a control."""
    explains = CASES["FRAME_EXPLAINS"]
    assert grade(explains).verdict == "FRAME_EXPLAINS"
    confounded = dataclasses.replace(explains, gc_plus_matched=explains.gc_plus_frame)
    assert grade(confounded).verdict == "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE"


# ------------------------------------------------------- the two arms, and each one's control


def test_the_c3_stratified_arm_can_fire_on_its_own():
    """A stratum margin collapsing is enough; C2 need not fire."""
    x = _with(atlas_in_window=0.32, gc_in_window=0.31)
    assert x.c3_absorbs and not x.c2_absorbs
    assert grade(x).verdict == "FRAME_EXPLAINS"


def test_a_collapsed_stratum_is_refused_when_random_stratification_alone_collapses_it():
    """Range restriction from stratifying is not absorption; the null is what separates them."""
    x = _with(atlas_in_window=0.32, gc_in_window=0.31, random_stratum_min_margin=0.01)
    assert grade(x).verdict == "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE"


def test_attributable_absorption_is_the_frame_minus_the_matched_control():
    x = _with(gc_plus_frame=0.40, gc_plus_matched=0.35)
    assert x.frame_absorption == pytest.approx(0.1964 - 0.1076, abs=1e-9)
    assert x.attributable_absorption == pytest.approx(
        x.frame_absorption - x.matched_absorption, abs=1e-12
    )
    assert x.attributable_absorption == pytest.approx(
        x.margin_over_gc_matched - x.margin_over_gc_frame, abs=1e-12
    )


def test_min_stratum_margin_takes_the_weaker_stratum_either_way_round():
    a = _with(atlas_in_window=0.32, gc_in_window=0.31)      # in-window is weaker
    b = _with(atlas_out_window=0.31, gc_out_window=0.30)    # out-of-window is weaker
    assert a.min_stratum_margin == pytest.approx(a.margin_in_window)
    assert b.min_stratum_margin == pytest.approx(b.margin_out_window)


# ------------------------------------------------------------------------------------- bounds


def test_the_excluded_verdict_carries_its_own_bound_so_a_quoter_cannot_drop_it():
    """The pattern that keeps EXPOSURE_EXCLUDED_LINEAR from being restated as 'excluded'."""
    d = grade(CASES["FRAME_EXCLUDED_CONDITIONAL"]).detail
    assert "AS REPLICATED" in d
    assert "200 bp" in d and "5 bp" in d
    assert "ONE organism" in d
    assert "never" in d.lower() and "position excluded" in d


def test_the_explains_verdict_does_not_claim_the_arm_metric_moves():
    d = grade(CASES["FRAME_EXPLAINS"]).detail
    assert "keeps its measured value" in d


def test_the_locate_floor_has_headroom_over_the_measured_rate():
    """0.95 was chosen above a measured 0.9918, not fitted down to it."""
    assert LOCATE_FLOOR == 0.95
    assert BASE.unique_locate_rate > LOCATE_FLOOR
