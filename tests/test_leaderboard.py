"""Tests for Step 17 — Shell-loop leaderboard."""
from __future__ import annotations

import pickle
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from scripts.leaderboard import (
    LeaderboardReport,
    LeaderboardRow,
    _fmt,
    run_classical_baseline_for_drug,
    run_foundation_model_for_drug,
)


def test_fmt_handles_none_and_signed():
    assert _fmt(None) == "—"
    assert _fmt(0.85) == "0.850"
    assert _fmt(0.12, signed=True) == "+0.120"
    assert _fmt(-0.05, signed=True) == "-0.050"


def test_leaderboard_report_markdown_groups_by_drug():
    report = LeaderboardReport()
    report.rows = [
        LeaderboardRow(drug="cipro", model_name="evo", is_classical=False, auroc=0.85),
        LeaderboardRow(drug="cipro", model_name="dnabert2", is_classical=False, auroc=0.78),
        LeaderboardRow(drug="cipro", model_name="kmer", is_classical=True, auroc=0.70),
        LeaderboardRow(drug="ceftriaxone", model_name="evo", is_classical=False, auroc=0.65),
    ]
    md = report.as_markdown()
    # Drug headers present
    assert "## cipro" in md
    assert "## ceftriaxone" in md
    # AUROC sorted descending within drug
    cipro_block = md.split("## cipro")[1].split("##")[0]
    evo_pos = cipro_block.find("| evo |")
    dnabert_pos = cipro_block.find("| dnabert2 |")
    assert evo_pos < dnabert_pos < cipro_block.find("| kmer |")
    # Type column populated
    assert "foundation" in md
    assert "classical" in md


def test_run_foundation_model_handles_train_failure(tmp_path: Path):
    """Failed subprocess → row with failure note + AUROC=None."""
    fail_result = MagicMock()
    fail_result.returncode = 1
    fail_result.stderr = "missing dependency"
    with patch("subprocess.run", return_value=fail_result):
        row = run_foundation_model_for_drug("evo", "cipro", tmp_path)
    assert row.auroc is None
    assert "train failed" in row.notes


def test_run_foundation_model_reads_bundle_metrics(tmp_path: Path):
    """Successful subprocess + bundle on disk → row with AUROC populated."""
    bundle_dir = tmp_path / "data" / "processed" / "models"
    bundle_dir.mkdir(parents=True)
    bundle_path = bundle_dir / "cipro_evo.pkl"
    with open(bundle_path, "wb") as f:
        pickle.dump({"auroc_loso": 0.91, "drug": "cipro", "model_name": "evo"}, f)

    ok_result = MagicMock()
    ok_result.returncode = 0
    ok_result.stderr = ""
    with patch("subprocess.run", return_value=ok_result):
        row = run_foundation_model_for_drug("evo", "cipro", tmp_path)
    assert row.auroc == pytest.approx(0.91)


def test_run_classical_baseline_missing_bundle(tmp_path: Path):
    """Missing bundle → row marked 'pending'."""
    row = run_classical_baseline_for_drug("classical_kmer", "cipro", tmp_path)
    assert row.is_classical is True
    assert "pending" in row.notes
    assert row.auroc is None


def test_run_classical_baseline_with_bundle(tmp_path: Path):
    """Bundle present → row with AUROC populated."""
    bundle_dir = tmp_path / "data" / "processed" / "models"
    bundle_dir.mkdir(parents=True)
    bundle_path = bundle_dir / "cipro_classical_kmer.pkl"
    with open(bundle_path, "wb") as f:
        pickle.dump({"auroc_loso": 0.72}, f)

    row = run_classical_baseline_for_drug("classical_kmer", "cipro", tmp_path)
    assert row.auroc == pytest.approx(0.72)
    assert row.is_classical is True
