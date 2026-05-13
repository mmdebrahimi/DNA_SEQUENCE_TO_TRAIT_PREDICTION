"""Tests for Step 15 — End-to-end smoke pipeline."""
from __future__ import annotations

from pathlib import Path

import pytest

xgboost = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
biopython = pytest.importorskip("Bio")


def test_smoke_pipeline_runs_to_completion(tmp_path: Path, project_root: Path):
    """End-to-end smoke pipeline completes in <60s on CPU + asserts AUROC + top-1."""
    from scripts.smoke_pipeline import run_smoke, SIGNAL_GENE

    results = run_smoke(
        fixtures_dir=project_root / "tests" / "fixtures" / "ecoli_mini",
        output_dir=tmp_path / "out",
    )

    # The seeded resistance signal in g1 should produce > random AUROC on
    # this synthetic 12-strain setup; we use a conservative 0.6 threshold
    # (LOSO on small N is noisy)
    assert results["auroc"] is not None
    assert results["auroc"] >= 0.6, f"AUROC {results['auroc']} too low on seeded signal"

    # Top-1 attribution should hit the seeded gene
    assert results["top_1_gene"] == SIGNAL_GENE, (
        f"top-1 attribution is {results['top_1_gene']!r}, expected {SIGNAL_GENE!r}"
    )

    # Report file written
    assert Path(results["report_path"]).exists()


def test_smoke_cli_exits_zero_on_success(tmp_path: Path, project_root: Path, monkeypatch):
    """CLI entry point exits 0 on a successful smoke run."""
    monkeypatch.chdir(project_root)

    from scripts.smoke_pipeline import main

    # Use the canonical fixture path + tmp output
    exit_code = main(
        [
            "--fixtures-dir",
            str(project_root / "tests" / "fixtures" / "ecoli_mini"),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    # 0 = PASS; 2 = wrong top-1 (synthetic signal didn't propagate); 3 = AUROC too low
    # We accept 0 only — anything else means the smoke regression breaks
    assert exit_code == 0, f"smoke CLI exited {exit_code}, expected 0"
