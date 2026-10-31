"""Verdict-logic pins for scripts/glm_selfdistill_learned_teacher.py.

Pure: no torch, no data, no network. `dna_decode.glm.encoder` lazy-imports torch inside
`fit_encoder`, so importing the probe module is safe to do while a torch job holds a model.

The load-bearing test here is `test_a_control_beating_its_teacher_does_not_manufacture_a_headline`.
It is not hypothetical -- the 1-seed smoke run measured `student_gc` at delta_rho **+0.0309** against
the ridge teacher, which EXCEEDS the 0.02 beat margin. A beats-check that did not exclude controls
would have reported DISTILLATION_HELPS_NOW, the single most consequential outcome in the bar, off a
control arm. The reason GC can do that is the same reason the v1 vacuity gate was mis-specified: on
the peak-tile rung GC alone is a strong scorer (0.3000), not an impoverished one.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from glm_selfdistill_learned_teacher import (  # noqa: E402
    AGREEMENT_FLOOR,
    BEAT_MARGIN,
    SHUF_MAX_ABS_AGREEMENT,
    verdict,
)

RIDGE = "ridge_positional_onehot"
CONV = "conv_encoder"


def arm(rho: float, agree: float, drho: float) -> dict:
    return {"rho_vs_measured": rho, "agreement_with_teacher": agree, "delta_rho_vs_teacher": drho}


def row(teacher: str, *, linear=(0.27, 1.0, 0.0), mlp=(0.07, 0.26, -0.19),
        shuffled=(-0.005, 0.013, -0.27), gc=(0.30, 0.75, 0.0), rate: float = 0.05) -> dict:
    return {
        "mutation_rate": rate, "teacher": teacher, "n_perturbed": 8000,
        "arms": {
            "student_linear": arm(*linear),
            "student_mlp": arm(*mlp),
            "student_shuffled": arm(*shuffled),
            "student_gc": arm(*gc),
        },
    }


def seed_rec(rows: list[dict], seed: int = 0) -> dict:
    return {
        "seed": seed, "held_out_peaks": 608, "n_train": 32908, "n_test": 11198,
        "n_sampled_for_perturbation": 8000,
        "teacher_rho": {RIDGE: 0.2691, CONV: 0.3235},
        "rows": rows,
    }


def both_teachers(*, conv_linear, ridge_linear=(0.2691, 1.0, 0.0), **kw) -> list[dict]:
    return [row(RIDGE, linear=ridge_linear, **kw), row(CONV, linear=conv_linear, **kw)]


# --------------------------------------------------------------------------------------------
# the headline outcomes
# --------------------------------------------------------------------------------------------

def test_dissolves_when_ridge_reproduces_and_conv_does_not():
    v, detail = verdict([seed_rec(both_teachers(conv_linear=(0.2677, 0.4279, -0.0558)))])
    assert v == "MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER"
    assert "0.4279" in detail


def test_persists_when_the_conv_teacher_is_also_reproduced():
    v, _ = verdict([seed_rec(both_teachers(conv_linear=(0.3235, 1.0, 0.0)))])
    assert v == "MECHANISM_PERSISTS_BAN_STANDS"


def test_agreement_exactly_at_the_floor_counts_as_reproduction():
    """The bar says 'at least 0.99'. A boundary that silently excluded its own threshold would
    turn a persist into a dissolve on an arithmetic accident."""
    v, _ = verdict([seed_rec(both_teachers(conv_linear=(0.3235, AGREEMENT_FLOOR, 0.0)))])
    assert v == "MECHANISM_PERSISTS_BAN_STANDS"


def test_high_agreement_but_drifted_rho_is_not_reproduction():
    """Both clauses are required: agreement AND |d rho| <= tol. Agreement alone is what a
    degenerate student also achieves, which is the whole reason the original probe logged both."""
    v, _ = verdict([seed_rec(both_teachers(conv_linear=(0.20, 0.999, -0.12)))])
    assert v == "MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER"


# --------------------------------------------------------------------------------------------
# the guards, in their frozen order
# --------------------------------------------------------------------------------------------

def test_vacuity_fires_when_the_shuffled_control_does_not_collapse():
    v, detail = verdict([seed_rec(both_teachers(
        conv_linear=(0.2677, 0.4279, -0.0558), shuffled=(0.31, 0.92, 0.0)))])
    assert v == "PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE"
    assert "NOT weakened a second time" in detail


def test_vacuity_fires_on_rho_alone_even_if_agreement_collapsed():
    v, _ = verdict([seed_rec(both_teachers(
        conv_linear=(0.2677, 0.4279, -0.0558), shuffled=(0.45, 0.01, 0.0)))])
    assert v == "PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE"


def test_vacuity_outranks_the_anchor_guard_and_the_headline():
    """Frozen check order. With vacuity, anchor failure and a dissolve all true at once, the
    earliest guard must win -- otherwise a broken probe reports a finding."""
    v, _ = verdict([seed_rec(both_teachers(
        ridge_linear=(0.10, 0.41, -0.17),            # anchor would fail
        conv_linear=(0.2677, 0.4279, -0.0558),       # conv would dissolve
        shuffled=(0.40, 0.95, 0.0)))])               # but the control did not collapse
    assert v == "PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE"


def test_anchor_guard_fires_when_the_ridge_teacher_is_not_degenerate_here():
    v, detail = verdict([seed_rec(both_teachers(
        ridge_linear=(0.11, 0.42, -0.16), conv_linear=(0.2677, 0.4279, -0.0558)))])
    assert v == "ANCHOR_FAILED_RIDGE_TEACHER_NOT_DEGENERATE_HERE"
    assert "UNINTERPRETABLE" in detail


def test_anchor_guard_outranks_the_headline():
    """A conv dissolve with a failed anchor must NOT be read as the headline: on a rung where the
    ridge teacher is not degenerate either, a conv dissolve is a property of the substrate."""
    v, _ = verdict([seed_rec(both_teachers(
        ridge_linear=(0.11, 0.42, -0.16), conv_linear=(0.3235, 1.0, 0.0)))])
    assert v == "ANCHOR_FAILED_RIDGE_TEACHER_NOT_DEGENERATE_HERE"


# --------------------------------------------------------------------------------------------
# the measured near-miss
# --------------------------------------------------------------------------------------------

def test_a_control_beating_its_teacher_does_not_manufacture_a_headline():
    """MEASURED, not hypothetical: smoke seed 0 had student_gc at delta_rho +0.0309 against the
    ridge teacher, past the 0.02 beat margin. DISTILLATION_HELPS_NOW is the outcome that would
    reopen a banned method; a control arm must never be able to trigger it."""
    assert 0.0309 > BEAT_MARGIN, "the premise of this test is that GC really did clear the margin"
    v, _ = verdict([seed_rec(both_teachers(
        conv_linear=(0.2677, 0.4279, -0.0558), gc=(0.30, 0.75, 0.0309)))])
    assert v == "MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER"


def test_the_shuffled_control_also_cannot_manufacture_a_headline():
    v, _ = verdict([seed_rec(both_teachers(
        conv_linear=(0.2677, 0.4279, -0.0558), shuffled=(-0.005, 0.013, 0.05)))])
    assert v == "MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER"


def test_the_beats_arm_is_not_vacuous_for_a_real_student():
    """Proving the exclusion above is a FILTER and not a disabled check: the same margin on a
    non-control student must still fire outcome 3."""
    v, detail = verdict([seed_rec(both_teachers(
        conv_linear=(0.3235, 0.80, 0.0309)))])
    assert v == "DISTILLATION_HELPS_NOW"
    assert "student_linear" in detail


def test_the_mlp_student_can_fire_the_beats_arm():
    """The MLP is the information-bound arm; a student exceeding its teacher there is the result
    that would contradict the measured bound, so it must not be filtered out either."""
    v, _ = verdict([seed_rec(both_teachers(
        conv_linear=(0.2677, 0.4279, -0.0558), mlp=(0.40, 0.70, 0.08)))])
    assert v == "DISTILLATION_HELPS_NOW"


# --------------------------------------------------------------------------------------------
# multi-seed / multi-rate aggregation
# --------------------------------------------------------------------------------------------

def test_one_non_reproducing_conv_cell_anywhere_is_enough_to_dissolve():
    """`all(reproduces)` over conv rows: a teacher that is degenerate at two rates and not at a
    third has not reproduced, and averaging would hide it."""
    rows = [row(RIDGE, rate=r) for r in (0.02, 0.05, 0.10)]
    rows += [row(CONV, linear=(0.3235, 1.0, 0.0), rate=0.02),
             row(CONV, linear=(0.3235, 1.0, 0.0), rate=0.05),
             row(CONV, linear=(0.2677, 0.43, -0.056), rate=0.10)]
    v, _ = verdict([seed_rec(rows)])
    assert v == "MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER"


def test_an_anchor_failure_in_any_seed_voids_the_whole_run():
    good = seed_rec(both_teachers(conv_linear=(0.2677, 0.4279, -0.0558)), seed=0)
    bad = seed_rec(both_teachers(ridge_linear=(0.1, 0.4, -0.17),
                                 conv_linear=(0.2677, 0.4279, -0.0558)), seed=1)
    v, detail = verdict([good, bad])
    assert v == "ANCHOR_FAILED_RIDGE_TEACHER_NOT_DEGENERATE_HERE"
    assert "(1," in detail, "the detail must name WHICH seed failed"


def test_thresholds_match_the_frozen_bar():
    """The bar is the authority; a drifted constant here would silently re-specify it."""
    import json
    bar = json.loads(Path("wiki/glm_selfdistill_learned_teacher_acceptance_bar.json")
                     .read_text(encoding="utf-8"))
    t = bar["thresholds_mirrored_from_2026_10_07"]
    assert AGREEMENT_FLOOR == t["agreement_floor"]
    assert BEAT_MARGIN == t["beat_margin"]
    assert SHUF_MAX_ABS_AGREEMENT == t["shuffled_control_max_abs_agreement"]


def test_the_bar_still_records_that_v1_was_mis_specified():
    """The amendment is only honest if its own record survives. A future tidy-up that deleted
    superseded_v1 would erase the fact that a control was replaced after a vacuous run."""
    import json
    bar = json.loads(Path("wiki/glm_selfdistill_learned_teacher_acceptance_bar.json")
                     .read_text(encoding="utf-8"))
    v1 = bar["superseded_v1"]
    assert "UNACHIEVABLE ON THIS RUNG BY CONSTRUCTION" in v1["what_changed"]
    assert v1["pre_commitment"]
    # the three headline thresholds must still be the mirrored ones, or the amendment was not
    # outcome-neutral as it claims
    t = bar["thresholds_mirrored_from_2026_10_07"]
    assert (t["agreement_floor"], t["delta_rho_tol"], t["beat_margin"]) == (0.99, 0.01, 0.02)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
