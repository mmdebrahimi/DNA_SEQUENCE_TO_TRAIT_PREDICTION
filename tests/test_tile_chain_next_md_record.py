"""`NEXT.md` records the tile chain's blocker, and HOW that blocker was discharged.

RE-PINNED 2026-10-27. This file was written while
`orthogonal-modality-stack-promcalc-plus-N-channels` was `planned` and SUBSTRATE-BLOCKED, and it
asserted the arm stays `planned` -- correctly, at the time: per R6 a declared arm cannot evaporate,
its blocker was substrate rather than merit, and demoting it would have destroyed the information
that it was still owed. **The substrate wall was then CLEARED (promcalc regenerated), the arm RAN,
and it is now `measured` at 0.4404 on its OWN executed evidence.** So these two guards went red
because the thing they guard legitimately moved -- the ONLY route off `planned` that R6 permits.

The invariant is unchanged and is what these tests now assert: an arm leaves `planned` by being
RUN, never by a memo, a sibling's result or a tidy-up. `tests/test_arms_ledger_stack_arm.py`
carries the metric pin and the empty-`evidence` rule for the three arms still planned.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
NEXT = REPO / "NEXT.md"
ARMS = REPO / "soraya_runs" / "arms-ledger.json"
ARM = "orthogonal-modality-stack-promcalc-plus-N-channels"


def test_next_md_names_the_promcalc_blocker():
    text = NEXT.read_text(encoding="utf-8")
    assert "promcalc" in text
    assert "data-absent" in text
    assert "code-closable" in text


def test_next_md_cites_the_memo_that_carries_the_derivation():
    assert "tile_chain_substrate_wall_2026-10-21" in NEXT.read_text(encoding="utf-8")


def test_next_md_records_the_correction_rather_than_only_the_finding():
    """A wrong belief the reader still holds is the thing worth correcting first."""
    text = NEXT.read_text(encoding="utf-8")
    assert "that was wrong" in text.lower()


def test_next_md_names_the_one_verified_runnable_item_this_leaves_standing():
    """A diagnosis that closes a path should say what is still open, or it reads as a dead end."""
    assert "Discharge_The_Queued_Sampling_Sweep" in NEXT.read_text(encoding="utf-8")


def test_next_md_says_the_frozen_bar_was_not_edited():
    text = NEXT.read_text(encoding="utf-8")
    assert "frozen bar was not edited" in text.lower() or "not edited either" in text.lower()


@pytest.mark.skipif(not ARMS.is_file(), reason="arms ledger absent")
def test_the_arm_left_planned_by_RUNNING_and_nothing_else():
    """R6: the one legitimate route off `planned` is the arm's own executed evidence.

    It is now `measured` at 0.4404, having been held `planned` through a sibling's `NO_GAIN`
    rather than demoted on it. Both halves matter: it was NOT killed on someone else's result,
    and it did NOT stay `planned` once it had a result of its own.

    The ledger's shape was READ, not guessed: each `criteria[<criterion text>]` maps arm NAME -> arm
    record. A first version of this test assumed `criteria[*]["arms"]` was a list and failed with the
    arm plainly present in the file -- a test asserting its author's model of the schema rather than
    the schema.
    """
    d = json.loads(ARMS.read_text(encoding="utf-8"))
    found = [arm for crit in d["criteria"].values()
             for name, arm in crit.items()
             if isinstance(arm, dict) and name == ARM]
    assert found, f"{ARM} not found under any criterion"
    assert all(a["status"] == "measured" for a in found)
    assert all(a["metric"] == 0.4404 for a in found)
    assert all(a["evidence"].strip() for a in found), \
        "measured with no evidence would mean the status moved without a run"
    assert all(a["cost_tier"] == "cheap" for a in found), \
        "the arm is still the cheap one; running it did not change its cost"


@pytest.mark.skipif(not ARMS.is_file(), reason="arms ledger absent")
def test_next_md_records_both_that_holding_it_was_right_and_that_it_then_RAN():
    """So a later reader can tell a discharged blocker from a tidied-away arm.

    Dropping the R6 half would leave the record saying only "it ran", losing why it was not
    killed when its sibling measured `NO_GAIN` -- which is the reusable half.
    """
    text = NEXT.read_text(encoding="utf-8")
    assert "R6" in text, "the reason it was not demoted must survive the update"
    assert "0.4404" in text, "the measured result must be in the record"
    assert "measured" in text.lower()
