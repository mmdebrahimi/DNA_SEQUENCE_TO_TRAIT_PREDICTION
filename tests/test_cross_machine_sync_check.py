"""Tests for scripts/cross_machine_sync_check.py pure-logic helpers."""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pytest


def test_render_markdown_in_sync():
    from scripts.cross_machine_sync_check import DriftReport, render_markdown
    r = DriftReport(run_date="2026-05-24", repo_path="/tmp/foo", in_sync=True)
    r.pytest_collected = 700
    md = render_markdown(r)
    assert "IN SYNC" in md
    assert "Clean." in md
    assert "**700** tests" in md


def test_render_markdown_drift_with_modifications():
    from scripts.cross_machine_sync_check import DriftReport, render_markdown
    r = DriftReport(
        run_date="2026-05-24",
        repo_path="/tmp/foo",
        in_sync=False,
        commit_gap_ahead=2,
        commit_gap_behind=5,
        working_tree_dirty=True,
        modified_files=["scripts/pipeline.py", "wiki/decoder_v0_ux_and_success_criterion.md"],
        untracked_files=["plans/new_plan.md"],
        downloads_recent=["cipro_bounded_falsifier_results_2026-05-24.json"],
        spec_divergences=[
            {"file": "scripts/pipeline.py", "marker": "cv_strategy", "status": "MISSING_LIKELY_STALE"},
        ],
    )
    md = render_markdown(r)
    assert "DRIFT DETECTED" in md
    assert "AHEAD of origin: **2**" in md
    assert "BEHIND origin: **5**" in md
    assert "Modified: 2" in md
    assert "scripts/pipeline.py" in md
    assert "Untracked: 1" in md
    assert "cipro_bounded_falsifier_results_2026-05-24.json" in md
    assert "MISSING_LIKELY_STALE" in md


def test_render_markdown_handles_large_modifications():
    """Truncates after 20 modified files + adds 'and N more' line."""
    from scripts.cross_machine_sync_check import DriftReport, render_markdown
    r = DriftReport(
        run_date="2026-05-24",
        repo_path="/tmp/foo",
        in_sync=False,
        working_tree_dirty=True,
        modified_files=[f"file_{i}.py" for i in range(30)],
    )
    md = render_markdown(r)
    assert "Modified: 30" in md
    assert "and 10 more" in md
    # First 20 listed
    assert "file_0.py" in md
    assert "file_19.py" in md


def test_check_downloads_recent_filters_by_name_pattern(tmp_path: Path):
    """Only files matching the project-artifact pattern get listed."""
    from scripts.cross_machine_sync_check import check_downloads_recent
    # Create test files
    (tmp_path / "cipro_results_2026-05-24.md").write_text("x")
    (tmp_path / "vacation_photo.jpg").write_text("x")
    (tmp_path / "dna_decoder_v0_handoff.md").write_text("x")
    (tmp_path / "unrelated.txt").write_text("x")  # filtered by name pattern
    (tmp_path / "falsifier_log.json").write_text("x")
    hits = check_downloads_recent(tmp_path)
    names = set(hits)
    assert "cipro_results_2026-05-24.md" in names
    assert "dna_decoder_v0_handoff.md" in names
    assert "falsifier_log.json" in names
    assert "vacation_photo.jpg" not in names
    assert "unrelated.txt" not in names


def test_check_downloads_recent_missing_dir_returns_empty(tmp_path: Path):
    from scripts.cross_machine_sync_check import check_downloads_recent
    missing = tmp_path / "does_not_exist"
    assert check_downloads_recent(missing) == []


def test_known_divergence_targets_pinned():
    """Pin the list so accidental edits show up in reviews."""
    from scripts.cross_machine_sync_check import KNOWN_DIVERGENCE_TARGETS
    expected = {
        # 2026-05-23 v0 closeout RELOCK markers
        ("wiki/decoder_v0_ux_and_success_criterion.md", "RELOCKED"),
        ("scripts/pipeline.py", "cv_strategy"),
        ("scripts/pipeline.py", "leave_one_accession_out"),
        # 2026-05-25 v0.1 slice 1 (genome-input cipro) markers
        ("scripts/pipeline.py", "--genome-fasta"),
        ("scripts/pipeline.py", "--allow-missing-audit"),
    }
    assert set(KNOWN_DIVERGENCE_TARGETS) == expected


def test_drift_report_json_roundtrip():
    """DriftReport -> dict -> json -> dict yields equivalent payload."""
    from scripts.cross_machine_sync_check import DriftReport
    r = DriftReport(
        run_date="2026-05-24",
        repo_path="/tmp/foo",
        commit_gap_ahead=1,
        modified_files=["a", "b"],
    )
    d = asdict(r)
    s = json.dumps(d)
    d2 = json.loads(s)
    assert d2["commit_gap_ahead"] == 1
    assert d2["modified_files"] == ["a", "b"]
