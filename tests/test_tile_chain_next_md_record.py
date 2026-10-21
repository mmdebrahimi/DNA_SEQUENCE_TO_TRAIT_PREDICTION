"""`NEXT.md` records the tile chain's real blocker, and the arms ledger was NOT edited.

The arms-ledger assertion is the load-bearing one. `soraya_runs/arms-ledger.json` keeps
`orthogonal-modality-stack-promcalc-plus-N-channels` as `planned`, and that is CORRECT: per R6 a
declared arm cannot evaporate, and this arm's blocker is substrate, not merit. Moving its status to
reflect a substrate problem would destroy the information that it is still owed -- so the right record
of a substrate blocker is a memo plus an expiry constraint, never a ledger edit.
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
def test_the_arm_is_still_planned_and_the_ledger_was_not_edited():
    """R6: a declared arm cannot evaporate, and a substrate blocker is not a reason to demote one.

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
    assert all(a["status"] == "planned" for a in found)
    assert all(a["cost_tier"] == "cheap" for a in found), \
        "the arm is still the cheap one; a substrate blocker does not change its cost"


@pytest.mark.skipif(not ARMS.is_file(), reason="arms ledger absent")
def test_next_md_states_that_planned_is_the_CORRECT_status():
    """So a later reader does not 'tidy' the ledger by demoting a substrate-blocked arm."""
    text = NEXT.read_text(encoding="utf-8")
    assert "stays `planned`" in text
    assert "R6" in text
