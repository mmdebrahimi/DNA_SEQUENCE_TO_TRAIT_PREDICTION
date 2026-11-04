"""Tests for the atomic sidecar write and the post-crash reconstruction.

THE FIRST GROUP SIMULATES THE CRASH THE FIX EXISTS TO SURVIVE. A plain `write_text` lost the whole
resume anchor to a power cut on 2026-11-03 (sidecar returned as 1,193 NUL bytes while the
atomically-written heartbeat in the same directory survived). So the test that matters is not "the
write works" but "a write interrupted between content and rename leaves the PREVIOUS sidecar
intact".

THE SECOND GROUP IS ALL REFUSALS. Reconstructing a sidecar means asserting how many rows of a block
are trustworthy, and asserting too many would score uninitialised zeros as measurements -- which
passes every downstream gate and corrupts the arm's number with no symptom. Every disagreement
between the three sources is therefore a refusal rather than a reconciliation.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import scripts.expression_heads_score as S
import scripts.reconstruct_heads_sidecar as R
from scripts.reconstruct_heads_sidecar import ReconstructRefused

HEADS = ["promoteratlas-trspred-hossain2020", "promoteratlas-trspred-lafleur2022",
         "promoteratlas-trspred-urtecho2018", "promoteratlas-trspred-yu2021",
         "promoteratlas-tslpred-kosuri2013"]


# --------------------------------------------------------------------- the atomic write


def test_write_sidecar_round_trips(tmp_path):
    S.write_sidecar(tmp_path, {"schema": "x", "status": "PARTIAL", "n_done": 64, "n_rows": 128})
    assert S.read_sidecar(tmp_path)["n_done"] == 64


def test_write_sidecar_leaves_no_tmp_behind(tmp_path):
    S.write_sidecar(tmp_path, {"status": "PARTIAL", "n_done": 1})
    assert list(tmp_path.glob("*.tmp")) == []


def test_an_INTERRUPTED_write_leaves_the_PREVIOUS_sidecar_INTACT(tmp_path, monkeypatch):
    """The property the power cut destroyed, and the reason tmp+replace is used.

    A crash between writing the new content and renaming it must leave the old sidecar readable.
    With a plain write_text the target is truncated first, so the same crash leaves nothing.
    """
    S.write_sidecar(tmp_path, {"status": "PARTIAL", "n_done": 64, "n_rows": 128})

    real_replace = os.replace

    def crash(src, dst):
        raise OSError("simulated power loss between fsync and rename")

    monkeypatch.setattr(os, "replace", crash)
    with pytest.raises(OSError):
        S.write_sidecar(tmp_path, {"status": "PARTIAL", "n_done": 128, "n_rows": 128})
    monkeypatch.setattr(os, "replace", real_replace)

    survived = S.read_sidecar(tmp_path)
    assert survived is not None, "the previous sidecar was destroyed by the failed write"
    assert survived["n_done"] == 64, "the previous sidecar's content changed"


def test_the_atomic_write_never_produces_a_NUL_file(tmp_path):
    """The crash artifact was 1,193 NUL bytes. The writer must never emit one."""
    S.write_sidecar(tmp_path, {"status": "COMPLETE", "n_done": 44106, "n_rows": 44106})
    raw = (tmp_path / S.SIDECAR_NAME).read_bytes()
    assert b"\x00" not in raw
    assert json.loads(raw.decode("utf-8"))["n_done"] == 44106


def test_a_NUL_filled_sidecar_still_reads_as_ABSENT(tmp_path):
    """The crash artifact must not parse. If it ever did, extract would resume from a number it
    invented rather than restarting."""
    (tmp_path / S.SIDECAR_NAME).write_bytes(b"\x00" * 1193)
    assert S.read_sidecar(tmp_path) is None


# --------------------------------------------------------------------- reading the block


def _block(tmp_path, n_rows=256, n_done=192, n_cols=10, *, gap=False, nonfinite=False):
    p = tmp_path / S.BLOCK_NAME
    arr = np.lib.format.open_memmap(p, mode="w+", dtype="float64", shape=(n_rows, n_cols))
    rng = np.random.default_rng(0)
    arr[:n_done] = rng.normal(size=(n_done, n_cols)) + 5.0      # +5 so no row is accidentally 0
    if gap:
        arr[n_done // 2] = 0.0                                   # a hole inside the filled region
    if nonfinite:
        arr[0, 0] = np.nan
    arr.flush()
    del arr
    return p


def test_the_block_reader_finds_the_written_boundary(tmp_path):
    p = _block(tmp_path, n_rows=256, n_done=192)
    got = R.filled_rows_from_block(p)
    assert got["first_all_zero_row"] == 192
    assert got["n_rows"] == 256
    assert got["tail_is_contiguously_zero"] is True
    assert got["filled_region_is_finite"] is True


def test_a_fully_written_block_reports_no_zero_row(tmp_path):
    p = _block(tmp_path, n_rows=128, n_done=128)
    got = R.filled_rows_from_block(p)
    assert got["first_all_zero_row"] == 128
    assert got["tail_is_contiguously_zero"] is True


def test_the_block_reader_detects_a_HOLE_rather_than_an_end(tmp_path):
    p = _block(tmp_path, n_rows=256, n_done=192, gap=True)
    got = R.filled_rows_from_block(p)
    assert got["tail_is_contiguously_zero"] is False


# --------------------------------------------------------------------- corroboration refusals


def _blk(n_rows=256, first_zero=192, *, tail=True, finite=True, cols=10):
    return {"n_rows": n_rows, "n_cols": cols, "first_all_zero_row": first_zero,
            "n_nonzero_rows": first_zero, "tail_is_contiguously_zero": tail,
            "filled_region_is_finite": finite}


def test_corroborate_agrees_when_all_three_sources_agree():
    got = R.corroborate(_blk(first_zero=192), 192, batch=64)
    assert got["n_done"] == 192 and got["status"] == "PARTIAL"
    assert len(got["agreed_by"]) == 3


def test_corroborate_marks_a_full_block_COMPLETE():
    got = R.corroborate(_blk(n_rows=192, first_zero=192), 192, batch=64)
    assert got["status"] == "COMPLETE"


def test_corroborate_REFUSES_a_hole_in_the_filled_region():
    with pytest.raises(ReconstructRefused, match="gap rather than the end"):
        R.corroborate(_blk(tail=False), 192, batch=64)


def test_corroborate_REFUSES_a_non_finite_filled_region():
    with pytest.raises(ReconstructRefused, match="non-finite"):
        R.corroborate(_blk(finite=False), 192, batch=64)


def test_corroborate_REFUSES_an_empty_block():
    with pytest.raises(ReconstructRefused, match="no written rows"):
        R.corroborate(_blk(first_zero=0), 0, batch=64)


def test_corroborate_REFUSES_a_count_that_is_not_a_batch_multiple():
    """The sidecar is written once per completed batch, so a non-multiple means a different loop
    produced this block and the tool's premise does not hold."""
    with pytest.raises(ReconstructRefused, match="not a multiple of the batch size"):
        R.corroborate(_blk(first_zero=190), 190, batch=64)


def test_a_non_multiple_is_allowed_only_when_it_is_the_FULL_row_count():
    """The last batch is short, so a complete block legitimately ends off-multiple."""
    got = R.corroborate(_blk(n_rows=190, first_zero=190), 190, batch=64)
    assert got["status"] == "COMPLETE" and got["n_done"] == 190


def test_corroborate_REFUSES_with_no_heartbeat_to_check_against():
    """One source is not corroboration. Writing a sidecar off the block alone would assert how
    many rows are trustworthy on unconfirmed evidence."""
    with pytest.raises(ReconstructRefused, match="ONE source"):
        R.corroborate(_blk(), None, batch=64)


def test_corroborate_REFUSES_a_DISAGREEING_heartbeat():
    """And explicitly does not reconcile by taking the minimum, which would be a guess."""
    with pytest.raises(ReconstructRefused, match="cannot tell which") as e:
        R.corroborate(_blk(first_zero=192), 128, batch=64)
    assert "guess" in str(e.value)


# --------------------------------------------------------------------- payload


def test_build_payload_derives_the_layout_from_the_head_names():
    payload = R.build_payload(R.corroborate(_blk(), 192, batch=64), HEADS, block=_blk())
    assert payload["layout"] == S.column_layout(HEADS)
    assert payload["layout"]["_total"] == 10
    assert payload["reconstructed"] is True


def test_build_payload_REFUSES_a_layout_that_disagrees_with_the_block():
    """A sidecar whose layout has the wrong column count maps every head to the wrong column."""
    with pytest.raises(ReconstructRefused, match="columns but the block has"):
        R.build_payload(R.corroborate(_blk(cols=6), 192, batch=64), HEADS, block=_blk(cols=6))


def test_the_payload_records_that_it_was_reconstructed_and_why():
    payload = R.build_payload(R.corroborate(_blk(), 192, batch=64), HEADS, block=_blk())
    note = payload["reconstruction_note"]
    assert "REBUILT" in note and "Kernel-Power 41" in note
    assert "block's own data is untouched" in note


def test_the_payload_is_readable_by_the_scorers_own_reader(tmp_path):
    """End to end: the reconstructed payload must satisfy `read_sidecar`, which is what `extract`
    actually consults."""
    payload = R.build_payload(R.corroborate(_blk(), 192, batch=64), HEADS, block=_blk())
    S.write_sidecar(tmp_path, payload)
    back = S.read_sidecar(tmp_path)
    assert back["n_done"] == 192 and back["status"] == "PARTIAL"


# --------------------------------------------------------------------- the real damaged block


REAL_BLOCK = S.BLOCK_DIR / S.BLOCK_NAME


def test_the_recovered_count_43328_is_a_clean_batch_MULTIPLE():
    """One of the three sources that agreed during the recovery, and the only one that is pure
    arithmetic rather than a property of a file: 43,328 is exactly 677 x 64. Kept as a standing
    fact because the other two sources no longer exist -- see the test below."""
    assert 43328 % S.BATCH == 0
    assert 43328 // S.BATCH == 677


@pytest.mark.skipif(not REAL_BLOCK.is_file(), reason="the real block is absent on this checkout")
def test_the_REAL_block_is_now_COMPLETE_and_no_longer_the_damaged_one():
    """CONVERTED 2026-11-04 from `test_the_REAL_damaged_block_corroborates_at_43328`.

    The original pinned a TRANSIENT state: the block stopped at 43,328 after a power cut, which
    three independent sources agreed on (the first all-zero row, the surviving heartbeat, and the
    batch arithmetic). The heads run has since finished -- `expression_heads` reads COMPLETE at
    44,106/44,106 -- so the block is fully written, `first_all_zero_row` returns the row count
    meaning "no zero tail", and the old assertion failed as `44106 == 43328`.

    That is the recovery having SUCCEEDED, not a regression, and the damaged block cannot be
    recreated. So this now pins the repaired state, and the recovery LOGIC stays pinned where it
    belongs: on the synthetic blocks above, which cover the boundary finder, holes, non-finite
    regions, non-multiple counts and a disagreeing heartbeat. A pin on a file that was damaged for
    one afternoon could never have survived, and keeping it red would have taught a reader to
    discount a red suite.
    """
    block = R.filled_rows_from_block(REAL_BLOCK)
    assert block["n_rows"] == 44106
    assert block["first_all_zero_row"] == 44106, "the block has a zero tail again"
    assert block["filled_region_is_finite"] is True
    agreed = R.corroborate(block, 44106, batch=S.BATCH)
    assert agreed["n_done"] == 44106 and agreed["status"] == "COMPLETE"
