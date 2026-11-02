"""The tile chain's blocker is registered with an auto-checkable expiry, and the store still parses.

This store is owned by a skill OUTSIDE this repo (`~/.claude/skills/soraya/scripts/constraint_expiry.py`)
and its schema is not pinned by this repo's tests, so a skill-side change could invalidate the row.
These tests bound that to DETECTION, not prevention -- which is the honest claim, and is why they
assert on shape rather than re-implementing the schema.

The expiry test itself has a limitation that is recorded rather than hidden, and
`test_the_rationale_states_the_expiry_limitation` is what keeps it recorded: a `file-exists` check on a
`D:` path fires only if that drive letter returns, so data arriving on another volume would not trip
it. The same honest limit was already recorded for R9, where the compute that eventually arrived was
cloud and a local-VRAM test could never have fired.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
STORE = REPO / "soraya_runs" / "constraint-expiry.json"
CID = "glm-tile-chain-promcalc-substrate-absent"

#: The four kinds the owning module defines. Restated here because this repo cannot import it.
KINDS = {"gpu-vram-gte", "file-exists", "test-exit-0", "manual"}


def _store() -> dict:
    return json.loads(STORE.read_text(encoding="utf-8"))


def test_the_store_still_parses():
    """A store that stopped parsing would take every constraint record with it."""
    d = _store()
    assert d["schema_version"] == 1
    assert isinstance(d["constraints"], dict)


def test_the_tile_chain_constraint_is_registered():
    assert CID in _store()["constraints"]


def test_the_expiry_kind_is_one_the_owning_module_defines():
    assert _store()["constraints"][CID]["expiry_kind"] in KINDS


def test_the_expiry_arg_is_a_path_not_prose():
    """A `file-exists` arg that is a sentence can never be evaluated by the checker."""
    arg = _store()["constraints"][CID]["expiry_arg"]
    assert "/" in arg
    assert " " not in arg
    assert arg.endswith(".gz")


def test_the_blocked_field_names_what_is_given_up_specifically_enough_to_act_on():
    blocked = _store()["constraints"][CID]["blocked"]
    assert "acceptance_bar" in blocked or "acceptance bar" in blocked.lower()
    assert "arms ledger" in blocked.lower() or "arms_ledger" in blocked


def test_the_blocked_field_also_records_what_is_NOT_blocked():
    """Half the chain was code-closable and is closed; a blocker that implied otherwise would mislead."""
    blocked = _store()["constraints"][CID]["blocked"]
    assert "--mpra-cache" in blocked
    assert "44106" in blocked


def test_the_rationale_states_the_expiry_limitation():
    """The limitation is the honest part. Deleting it would leave a check that looks stronger than it is."""
    r = _store()["constraints"][CID]["rationale"]
    assert "LIMITATION" in r
    assert "another volume" in r


def test_the_rationale_carries_the_drive_independent_complete_test():
    """`test-exit-0` is not auto-run by check_all, so the complete test lives in the text."""
    r = _store()["constraints"][CID]["rationale"]
    assert "promcalc_toolchain_probe" in r
    assert "ENV_PRESENT" in r


def test_the_rationale_names_the_anchor_risk():
    """Regeneration is reachable; reproducing the frozen anchors is a separate, stricter question."""
    assert "REFUSED_ANCHOR_MISMATCH" in _store()["constraints"][CID]["rationale"]


#: Every constraint id this repo has seen registered. NAMED INDIVIDUALLY on purpose -- the original
#: docstring is right that "a count alone would pass if one were swapped", and the names are the
#: actual guard. `glm-tile-chain-promcalc-substrate-absent` is `CID` and is also asserted by its own
#: test above; the redundancy is harmless and keeps this list a complete inventory.
#:
#: **The trailing total used to be `== 12` and had been RED against a registry of 14.** That is the
#: same defect class the arms guards carried: an exact count standing where an invariant belongs. The
#: count was 11 named + CID when written, and two rows have been appended since
#: (`gd-sweep-weights-absent`, `mvp-bar-parser-blind-to-table-format`). Both are now NAMED rather
#: than absorbed into a looser bound, so the file keeps its named-individually property, and the
#: total is a FLOOR TIED TO THIS LIST so it cannot go stale on the next append.
#:
#: The test's name was also made count-free. It read "eleven" while naming eleven of fourteen rows,
#: and a guard whose name misstates what it guards is the drift this project keeps finding in
#: orientation strings.
NAMED_CONSTRAINTS = (
    "north-star-step4-animals-deferred",
    "ecoff-values-still-unsourced",
    "f5-etd-route-retired-with-remainder",
    "atlas-pfizer-vivli-request-gated",
    "align-foundation-benchmark-is-a-proposal",
    "cabbage-genomes-tree-coverage-uncensused",
    "cabbage-clonality-g8-blocked-no-mash",
    "local-gpu-vram-blocks-dna-lm-finetune",
    "gd-finetune-gated-on-sampling-sweep-not-compute",
    "env-python-pinned-313-no-315-wheels",
    "glm-fragment-rung-kmer6-full-size-ram-bound",
    "gd-sweep-weights-absent",
    "mvp-bar-parser-blind-to-table-format",
    CID,
)


def assert_named_constraints_present(rows: dict) -> None:
    """The guard, taking its store as an argument so a synthetic one can prove it non-vacuous."""
    for cid in NAMED_CONSTRAINTS:
        assert cid in rows, f"pre-existing constraint {cid} vanished"
    assert len(rows) >= len(NAMED_CONSTRAINTS), (
        f"the registry holds {len(rows)} rows but {len(NAMED_CONSTRAINTS)} are named; a row was "
        f"removed without this list being updated"
    )


def test_the_pre_existing_constraints_are_still_present():
    assert_named_constraints_present(_store()["constraints"])


def test_the_presence_guard_is_NOT_VACUOUS_a_removed_constraint_is_detected():
    """What the floor must not cost. A bound alone would pass a swap; the names are the guard."""
    rows = dict(_store()["constraints"])
    victim = NAMED_CONSTRAINTS[0]
    rows.pop(victim)
    with pytest.raises(AssertionError, match=victim):
        assert_named_constraints_present(rows)


def test_the_floor_does_not_drift_on_a_further_append():
    """The whole reason the total is a floor rather than an equality."""
    rows = dict(_store()["constraints"])
    rows["a-synthetic-future-constraint"] = {"expiry_kind": "manual"}
    assert_named_constraints_present(rows)          # must not raise
