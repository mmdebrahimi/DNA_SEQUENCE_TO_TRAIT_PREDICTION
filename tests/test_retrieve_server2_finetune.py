"""Tests for bringing the fine-tune result back from server2.

ONLY THE VERIFICATION IS TESTED, and that is the deliberate scope. The copy itself is `scp` against
a host that may not be reachable from wherever this runs, so wrapping it in a mock would test the
mock. What matters and IS testable is the half that decides whether what arrived can be trusted:
scp exits 0 on this host while transferring nothing, so `verify_local` is the only thing standing
between a silent non-transfer and a scored arm.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import scripts.retrieve_server2_finetune as R
from scripts.retrieve_server2_finetune import RUN_ARTIFACT, RetrieveRefused, verify_local

SENTINEL_REL = R.STATE_FILES[0]


def _stage(tmp_path, *, run=None, sentinel="EXIT=0\r\n", write_run=True, write_sentinel=True):
    art = tmp_path / RUN_ARTIFACT
    art.parent.mkdir(parents=True, exist_ok=True)
    if write_run:
        payload = {"schema": "finetune-arm-run-v1", "status": "COMPLETE",
                   "horizon_epochs": 21, "per_seed": {str(i): {} for i in range(10)}}
        if run is not None:
            payload.update(run)
        art.write_text(json.dumps(payload), encoding="utf-8")
    s = tmp_path / SENTINEL_REL
    s.parent.mkdir(parents=True, exist_ok=True)
    if write_sentinel:
        s.write_text(sentinel, encoding="utf-8")
    return tmp_path


def test_a_complete_run_verifies(tmp_path):
    info = verify_local(repo=_stage(tmp_path))
    assert info["schema"] == "finetune-arm-run-v1"
    assert info["status"] == "COMPLETE"
    assert info["n_seeds"] == 10
    assert info["horizon_epochs"] == 21
    assert info["sentinel"] == "EXIT=0"


def test_an_EMPTY_sentinel_is_REFUSED(tmp_path):
    """The documented `echo EXIT=%RC%>>file` trap: cmd reads the digit before `>>` as a file
    handle and leaves the sentinel empty, which is indistinguishable from a clean exit."""
    with pytest.raises(RetrieveRefused, match="EMPTY"):
        verify_local(repo=_stage(tmp_path, sentinel="   \r\n"))


def test_an_ABSENT_sentinel_is_REFUSED(tmp_path):
    """Without it the launcher state cannot read COMPLETE, so the collector would refuse the arm
    while the number sat on disk -- which looks like the run never happened."""
    with pytest.raises(RetrieveRefused, match="sentinel is absent"):
        verify_local(repo=_stage(tmp_path, write_sentinel=False))


def test_an_absent_run_artifact_is_REFUSED(tmp_path):
    with pytest.raises(RetrieveRefused, match="absent"):
        verify_local(repo=_stage(tmp_path, write_run=False))


def test_an_unreadable_run_artifact_is_REFUSED(tmp_path):
    p = tmp_path / RUN_ARTIFACT
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text('{"schema": "finetune-arm-run-v1"', encoding="utf-8")
    (tmp_path / SENTINEL_REL).parent.mkdir(parents=True, exist_ok=True)
    (tmp_path / SENTINEL_REL).write_text("EXIT=0", encoding="utf-8")
    with pytest.raises(RetrieveRefused, match="not readable JSON"):
        verify_local(repo=tmp_path)


def test_a_WRONG_SCHEMA_is_REFUSED(tmp_path):
    with pytest.raises(RetrieveRefused, match="schema"):
        verify_local(repo=_stage(tmp_path, run={"schema": "something-else-v9"}))


def test_a_NONZERO_sentinel_still_verifies_and_is_REPORTED(tmp_path):
    """A failed run is a fact to carry forward, not an error to hide: the scorer's own gates decide
    what a non-zero exit means, and they can only do that if the exit code arrives."""
    info = verify_local(repo=_stage(tmp_path, sentinel="EXIT=1\r\n"))
    assert info["sentinel"] == "EXIT=1"


def test_a_BOUNDED_run_verifies(tmp_path):
    """A bound in the ineligible tail leaves a readable verdict, so retrieval must not refuse it."""
    info = verify_local(repo=_stage(tmp_path, run={"status": "BOUNDED_BY_WALL_CLOCK"}))
    assert info["status"] == "BOUNDED_BY_WALL_CLOCK"


def test_the_state_files_include_the_sentinel_first():
    """`verify_local` indexes STATE_FILES[0] as the sentinel; a reorder would silently check the
    wrong file for emptiness."""
    assert "sentinel" in R.STATE_FILES[0]


def test_every_state_file_is_under_data_unattended():
    for rel in R.STATE_FILES:
        assert rel.startswith("data/unattended/"), rel
