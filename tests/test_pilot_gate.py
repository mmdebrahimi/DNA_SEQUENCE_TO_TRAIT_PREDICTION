"""Tests for Step 0.5 — Real-data pilot gate."""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.data import pilot as pilot_module
from dna_decode.data.pilot import (
    CohortSelectionCriteria,
    PilotGateError,
    PilotReport,
    apply_assembly_filter,
    apply_method_filter,
    estimate_three_drug_intersection,
    run_pilot_gate,
    write_pilot_report,
)


# ---- filter-math unit tests ----


def test_apply_method_filter_default_factor():
    raw = {"cipro": 1000, "ceftriaxone": 800, "tet": 500}
    out = apply_method_filter(raw)
    # Default 0.6 survival
    assert out["cipro"] == 600
    assert out["ceftriaxone"] == 480
    assert out["tet"] == 300


def test_apply_assembly_filter_default_factor():
    method_filtered = {"cipro": 600, "ceftriaxone": 480}
    out = apply_assembly_filter(method_filtered)
    # Default 0.85 survival
    assert out["cipro"] == 510
    assert out["ceftriaxone"] == 408


def test_estimate_three_drug_intersection_min_strain_basis():
    """Intersection bounded by min-strain drug × overlap factor."""
    counts = {"cipro": 510, "ceftriaxone": 408, "tet": 280}
    # min=280, factor=0.45 -> 126
    assert estimate_three_drug_intersection(counts) == 126


def test_estimate_three_drug_intersection_empty():
    assert estimate_three_drug_intersection({}) == 0


# ---- run_pilot_gate integration tests (mocked HTTP) ----


def test_run_pilot_gate_go_verdict(monkeypatch, project_root: Path):
    """All drugs above target → GO."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg: {d: 1000 for d in drugs},
    )

    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    report = run_pilot_gate(
        drugs=drugs,
        criteria=CohortSelectionCriteria(target_per_drug=150, three_drug_intersection_target=75),
        config_path=project_root / "config" / "datasources.yaml",
    )

    assert report.go_no_go == "GO"
    assert report.failure_reasons == []
    # 1000 * 0.6 * 0.85 = 510 after both filters
    for stage in report.per_drug:
        assert stage.after_assembly_filter == 510
    # 3-drug intersection: min=510, factor=0.45 -> 229
    assert report.three_drug_intersection == 229


def test_run_pilot_gate_no_go_per_drug_short(monkeypatch, project_root: Path):
    """At least one drug below target → NO-GO with diagnostic."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg: {"ciprofloxacin": 1000, "ceftriaxone": 200, "tetracycline": 800},
    )

    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    report = run_pilot_gate(
        drugs=drugs,
        criteria=CohortSelectionCriteria(target_per_drug=150),
        config_path=project_root / "config" / "datasources.yaml",
    )

    # ceftriaxone: 200 * 0.6 * 0.85 = 102, below 150 target
    assert report.go_no_go == "NO-GO"
    assert any("ceftriaxone" in r for r in report.failure_reasons)


def test_run_pilot_gate_no_go_intersection_short(monkeypatch, project_root: Path):
    """All drugs pass per-drug but intersection fails."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg: {d: 400 for d in drugs},
    )

    # 400 * 0.6 * 0.85 = 204 per drug (passes 150)
    # intersection: 204 * 0.45 = 91. With a high target it fails:
    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    report = run_pilot_gate(
        drugs=drugs,
        criteria=CohortSelectionCriteria(
            target_per_drug=150, three_drug_intersection_target=200
        ),
        config_path=project_root / "config" / "datasources.yaml",
    )

    assert report.go_no_go == "NO-GO"
    assert any("intersection" in r for r in report.failure_reasons)


def test_run_pilot_gate_config_missing_raises(tmp_path: Path):
    """Bad config path raises PilotGateError."""
    bogus = tmp_path / "missing.yaml"
    with pytest.raises(PilotGateError, match="config not found"):
        run_pilot_gate(drugs=("cipro",), config_path=bogus)


# ---- report-writing tests ----


def test_write_pilot_report_go(tmp_path: Path):
    """GO report has expected sections."""
    report = PilotReport(
        drugs=("cipro",),
        criteria=CohortSelectionCriteria(),
        per_drug=[],
        three_drug_intersection=200,
        go_no_go="GO",
    )
    out = write_pilot_report(report, tmp_path / "report.md")
    content = out.read_text(encoding="utf-8")
    assert "Pilot Gate Report" in content
    assert "GO" in content
    assert "Failure reasons" not in content


def test_write_pilot_report_no_go_includes_failures(tmp_path: Path):
    """NO-GO report includes failure-reasons section + remediation hint."""
    report = PilotReport(
        drugs=("ceftriaxone",),
        criteria=CohortSelectionCriteria(),
        per_drug=[],
        three_drug_intersection=50,
        go_no_go="NO-GO",
        failure_reasons=["ceftriaxone: 100 strains after filters < target 150"],
    )
    out = write_pilot_report(report, tmp_path / "report.md")
    content = out.read_text(encoding="utf-8")
    assert "NO-GO" in content
    assert "Failure reasons" in content
    assert "ceftriaxone" in content
    assert "remediation" in content.lower()


# ---- CLI exit-code semantics ----


def test_cli_go_exits_zero(monkeypatch, tmp_path: Path, project_root: Path):
    """Successful pilot gate exits 0."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg: {d: 1000 for d in drugs},
    )
    monkeypatch.chdir(project_root)

    from scripts.pilot_gate import main

    exit_code = main(
        [
            "--drugs",
            "ciprofloxacin,ceftriaxone,tetracycline",
            "--output",
            str(tmp_path / "report.md"),
        ]
    )
    assert exit_code == 0


def test_cli_no_go_exits_nonzero(monkeypatch, tmp_path: Path, project_root: Path):
    """Failed pilot gate exits 1 (HARD-gate semantic)."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg: {d: 100 for d in drugs},  # too few
    )
    monkeypatch.chdir(project_root)

    from scripts.pilot_gate import main

    exit_code = main(
        [
            "--drugs",
            "ciprofloxacin,ceftriaxone,tetracycline",
            "--output",
            str(tmp_path / "report.md"),
        ]
    )
    assert exit_code == 1


def test_cli_scaffold_exits_three(tmp_path: Path, project_root: Path, monkeypatch):
    """Without monkeypatched fetcher, scaffold raises NotImplementedError → exit 3."""
    monkeypatch.chdir(project_root)

    from scripts.pilot_gate import main

    exit_code = main(
        [
            "--drugs",
            "ciprofloxacin,ceftriaxone,tetracycline",
            "--output",
            str(tmp_path / "report.md"),
        ]
    )
    assert exit_code == 3
