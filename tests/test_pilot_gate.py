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


def test_apply_method_filter_pass_through_by_default():
    """Wave 1.5: local-TSV path applies method filter in load_bvbrc_ast; pilot
    apply_method_filter is pass-through (factor=1.0) by default."""
    raw = {"cipro": 1000, "ceftriaxone": 800, "tet": 500}
    out = apply_method_filter(raw)
    assert out["cipro"] == 1000
    assert out["ceftriaxone"] == 800
    assert out["tet"] == 500


def test_apply_method_filter_explicit_factor():
    """When an explicit factor is supplied (e.g., for future live-API path)."""
    raw = {"cipro": 1000}
    out = apply_method_filter(raw, method_filter_factor=0.6)
    assert out["cipro"] == 600


def test_apply_assembly_filter_default_factor():
    method_filtered = {"cipro": 1000, "ceftriaxone": 800}
    out = apply_assembly_filter(method_filtered)
    # Default 0.85 survival
    assert out["cipro"] == 850
    assert out["ceftriaxone"] == 680


def test_estimate_three_drug_intersection_min_strain_basis():
    """Intersection bounded by min-strain drug × overlap factor."""
    counts = {"cipro": 850, "ceftriaxone": 680, "tet": 425}
    # min=425, factor=0.45 -> 191
    assert estimate_three_drug_intersection(counts) == 191


def test_estimate_three_drug_intersection_empty():
    assert estimate_three_drug_intersection({}) == 0


# ---- run_pilot_gate integration tests (mocked HTTP) ----


def test_run_pilot_gate_go_verdict(monkeypatch, project_root: Path):
    """All drugs above target → GO."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg, ast_tsv_path=None: {d: 1000 for d in drugs},
    )

    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    report = run_pilot_gate(
        drugs=drugs,
        criteria=CohortSelectionCriteria(target_per_drug=150, three_drug_intersection_target=75),
        config_path=project_root / "config" / "datasources.yaml",
    )

    assert report.go_no_go == "GO"
    assert report.failure_reasons == []
    # Wave 1.5: local-TSV path already filters; method-filter is pass-through
    # 1000 * 1.0 * 0.85 = 850 after both filters
    for stage in report.per_drug:
        assert stage.after_assembly_filter == 850
    # 3-drug intersection: min=850, factor=0.45 -> 382
    assert report.three_drug_intersection == 382


def test_run_pilot_gate_no_go_per_drug_short(monkeypatch, project_root: Path):
    """At least one drug below target → NO-GO with diagnostic."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg, ast_tsv_path=None: {
            "ciprofloxacin": 1000,
            "ceftriaxone": 100,
            "tetracycline": 800,
        },
    )

    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    report = run_pilot_gate(
        drugs=drugs,
        criteria=CohortSelectionCriteria(target_per_drug=150),
        config_path=project_root / "config" / "datasources.yaml",
    )

    # ceftriaxone: 100 * 1.0 * 0.85 = 85, below 150 target
    assert report.go_no_go == "NO-GO"
    assert any("ceftriaxone" in r for r in report.failure_reasons)


def test_run_pilot_gate_no_go_intersection_short(monkeypatch, project_root: Path):
    """All drugs pass per-drug but intersection fails."""
    monkeypatch.setattr(
        pilot_module,
        "fetch_bvbrc_drug_counts",
        lambda drugs, cfg, ast_tsv_path=None: {d: 250 for d in drugs},
    )

    # 250 * 1.0 * 0.85 = 212 per drug (passes 150)
    # intersection: 212 * 0.45 = 95. With a high target it fails:
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


# ---- Local-TSV fetch (Wave 1.5 hardening) ----


_MOCK_BVBRC_TSV = """\
genome_id\tgenome_name\tantibiotic\tresistant_phenotype\tmeasurement\tmeasurement_unit\tlaboratory_typing_method\ttesting_standard
strain001\tEscherichia coli\tciprofloxacin\tResistant\t8\tmg/L\tBroth microdilution\tCLSI
strain002\tEscherichia coli\tciprofloxacin\tSusceptible\t0.06\tmg/L\tBroth microdilution\tEUCAST
strain003\tEscherichia coli\tceftriaxone\tResistant\t32\tmg/L\tBroth microdilution\tCLSI
strain004\tEscherichia coli\ttetracycline\tResistant\t16\tmg/L\tBroth microdilution\tCLSI
strain005\tEscherichia coli\tciprofloxacin\tResistant\t8\tmg/L\tDisk diffusion\tCLSI
"""


def test_fetch_bvbrc_drug_counts_from_local_tsv(tmp_path: Path, project_root: Path):
    """Local-TSV path counts per-drug + applies broth-microdilution filter."""
    tsv = tmp_path / "ast.tsv"
    tsv.write_text(_MOCK_BVBRC_TSV, encoding="utf-8")

    import yaml
    with open(project_root / "config" / "datasources.yaml") as f:
        cfg = yaml.safe_load(f)

    counts = pilot_module.fetch_bvbrc_drug_counts(
        ("ciprofloxacin", "ceftriaxone", "tetracycline"), cfg, ast_tsv_path=tsv
    )
    # strain005 is disk-diffusion → filtered. Cipro: 2 (strain001, strain002)
    # Ceftriaxone: 1 (strain003). Tetracycline: 1 (strain004).
    assert counts["ciprofloxacin"] == 2
    assert counts["ceftriaxone"] == 1
    assert counts["tetracycline"] == 1


def test_fetch_bvbrc_drug_counts_via_env_var(tmp_path: Path, project_root: Path, monkeypatch):
    """BVBRC_AST_TSV env var triggers local path."""
    tsv = tmp_path / "ast.tsv"
    tsv.write_text(_MOCK_BVBRC_TSV, encoding="utf-8")
    monkeypatch.setenv("BVBRC_AST_TSV", str(tsv))

    import yaml
    with open(project_root / "config" / "datasources.yaml") as f:
        cfg = yaml.safe_load(f)

    counts = pilot_module.fetch_bvbrc_drug_counts(("ciprofloxacin",), cfg)
    assert counts["ciprofloxacin"] == 2


def test_fetch_bvbrc_drug_counts_no_source_raises():
    """No TSV path AND no env var AND no config entry → NotImplementedError."""
    cfg = {"bvbrc_ast": {}}
    with pytest.raises(NotImplementedError, match="Live BV-BRC API"):
        pilot_module.fetch_bvbrc_drug_counts(("cipro",), cfg)


def test_fetch_bvbrc_drug_counts_missing_tsv_raises(tmp_path: Path):
    """Path provided but file doesn't exist → PilotGateError."""
    cfg = {"bvbrc_ast": {}}
    with pytest.raises(PilotGateError, match="not found"):
        pilot_module.fetch_bvbrc_drug_counts(
            ("cipro",), cfg, ast_tsv_path=tmp_path / "missing.tsv"
        )


def test_run_pilot_gate_end_to_end_local_tsv(tmp_path: Path, project_root: Path):
    """End-to-end pilot gate run from a local TSV (no monkeypatching)."""
    # Build a TSV with enough resistant strains for cipro
    rows = ["\t".join(["genome_id", "genome_name", "antibiotic", "resistant_phenotype",
                       "measurement", "measurement_unit", "laboratory_typing_method",
                       "testing_standard"])]
    for i in range(200):
        rows.append(
            f"s{i:04d}\tEscherichia coli\tciprofloxacin\tResistant\t8\tmg/L\tBroth microdilution\tCLSI"
        )
    for i in range(200):
        rows.append(
            f"s{i+200:04d}\tEscherichia coli\tceftriaxone\tResistant\t32\tmg/L\tBroth microdilution\tCLSI"
        )
    for i in range(200):
        rows.append(
            f"s{i+400:04d}\tEscherichia coli\ttetracycline\tResistant\t16\tmg/L\tBroth microdilution\tCLSI"
        )
    tsv = tmp_path / "ast.tsv"
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")

    report = run_pilot_gate(
        drugs=("ciprofloxacin", "ceftriaxone", "tetracycline"),
        criteria=CohortSelectionCriteria(
            target_per_drug=150, three_drug_intersection_target=75
        ),
        config_path=project_root / "config" / "datasources.yaml",
        ast_tsv_path=tsv,
    )
    # 200 per drug; 200 * 0.85 = 170 after assembly filter → above 150 target
    assert report.go_no_go == "GO"
    for stage in report.per_drug:
        assert stage.raw == 200
        assert stage.after_assembly_filter == 170


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
        lambda drugs, cfg, ast_tsv_path=None: {d: 1000 for d in drugs},
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
        lambda drugs, cfg, ast_tsv_path=None: {d: 100 for d in drugs},  # too few
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


def test_cli_no_ast_source_exits_three(tmp_path: Path, project_root: Path, monkeypatch):
    """No --ast-tsv + no env var + no config entry → NotImplementedError → exit 3."""
    monkeypatch.chdir(project_root)
    # Make sure no env var leaks in
    monkeypatch.delenv("BVBRC_AST_TSV", raising=False)

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


def test_cli_ast_tsv_flag_drives_real_path(tmp_path: Path, project_root: Path, monkeypatch):
    """--ast-tsv flag → real local-TSV path; no monkeypatching of fetch_bvbrc_drug_counts."""
    tsv = tmp_path / "ast.tsv"
    # Build a TSV that produces a GO verdict (200 per drug)
    rows = ["\t".join(["genome_id", "genome_name", "antibiotic", "resistant_phenotype",
                       "measurement", "measurement_unit", "laboratory_typing_method",
                       "testing_standard"])]
    for drug in ("ciprofloxacin", "ceftriaxone", "tetracycline"):
        for i in range(200):
            rows.append(
                f"s_{drug}_{i:04d}\tEscherichia coli\t{drug}\tResistant\t8\tmg/L\tBroth microdilution\tCLSI"
            )
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    monkeypatch.chdir(project_root)

    from scripts.pilot_gate import main

    exit_code = main(
        [
            "--drugs",
            "ciprofloxacin,ceftriaxone,tetracycline",
            "--target-per-drug",
            "150",
            "--ast-tsv",
            str(tsv),
            "--output",
            str(tmp_path / "report.md"),
        ]
    )
    assert exit_code == 0
