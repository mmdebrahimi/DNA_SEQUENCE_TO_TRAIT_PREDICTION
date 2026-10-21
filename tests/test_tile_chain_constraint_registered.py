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


def test_the_eleven_pre_existing_constraints_are_still_present():
    """This was an APPEND. Named individually, because a count alone would pass if one were swapped."""
    rows = _store()["constraints"]
    for cid in (
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
    ):
        assert cid in rows, f"pre-existing constraint {cid} vanished"
    assert len(rows) == 12
