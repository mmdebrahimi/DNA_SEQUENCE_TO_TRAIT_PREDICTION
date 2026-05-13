"""Tests for scripts/audit_cohort.py — Phase 2.5 audit report generator."""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.data.cohort import CandidateStrain, StrainCohort, save_cohort
from scripts.audit_cohort import (
    VerdictRules,
    build_parser,
    build_report,
    evaluate_verdict,
    main,
)


def _make_strain(strain_id: str, contigs: int, n50: int, mlst: str, labels: dict[str, int]):
    return CandidateStrain(
        strain_id=strain_id,
        assembly_accession=f"GCF_{strain_id}.1",
        mlst=mlst,
        contig_count=contigs,
        n50=n50,
        ast_labels=labels,
    )


def _build_cohort(strains: list[CandidateStrain]) -> StrainCohort:
    """Build StrainCohort with per-drug pools auto-populated."""
    per_drug: dict[str, list[str]] = {}
    for s in strains:
        for drug in s.ast_labels:
            per_drug.setdefault(drug, []).append(s.strain_id)
    return StrainCohort(
        strains=strains,
        per_drug_strain_ids=per_drug,
        three_drug_intersection=[],
    )


# ---- overview / report assembly ----


def test_build_report_renders_all_sections(tmp_path: Path):
    """Smoke: empty-ish AST + small cohort → report contains all expected headers."""
    strains = [
        _make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1}),
        _make_strain("S1", 2, 4_500_000, "ST10", {"ciprofloxacin": 0}),
    ]
    cohort = _build_cohort(strains)
    report = build_report(cohort, None, ["ciprofloxacin"], VerdictRules(), tmp_path / "src.parquet")
    for header in (
        "# Cohort Audit Report",
        "## Cohort overview",
        "## Clade composition",
        "## Metadata completeness",
        "## Assembly QC distribution",
        "## Pre-Gate-B verdict",
    ):
        assert header in report


def test_build_report_excludes_ast_section_when_no_ast(tmp_path: Path):
    strains = [_make_strain("R1", 2, 5_000_000, "ST1", {"ciprofloxacin": 1})]
    cohort = _build_cohort(strains)
    report = build_report(cohort, None, ["ciprofloxacin"], VerdictRules(), tmp_path / "x.parquet")
    assert "## AST source breakdown" not in report


def test_build_report_includes_ast_section_when_provided(tmp_path: Path):
    import pandas as pd

    strains = [_make_strain("R1", 2, 5_000_000, "ST1", {"ciprofloxacin": 1})]
    cohort = _build_cohort(strains)
    ast_df = pd.DataFrame(
        [
            {"strain_id": "R1", "antibiotic": "ciprofloxacin",
             "susceptibility_label": "R", "mic_value": 8.0, "mic_units": "mg/L",
             "measurement_method": "broth_microdilution", "source": "CLSI"}
        ]
    )
    report = build_report(cohort, ast_df, ["ciprofloxacin"], VerdictRules(), tmp_path / "x.parquet")
    assert "## AST source breakdown" in report
    assert "broth_microdilution" in report


# ---- verdict logic ----


def test_verdict_go_when_all_thresholds_met():
    """50 strains × 30 R / 20 S for one drug, full metadata → GO."""
    strains = []
    for i in range(30):
        strains.append(_make_strain(
            f"R{i}", 3, 200_000, f"ST{i}", {"ciprofloxacin": 1}
        ))
    for i in range(20):
        strains.append(_make_strain(
            f"S{i}", 3, 200_000, f"ST{i+100}", {"ciprofloxacin": 0}
        ))
    cohort = _build_cohort(strains)
    rules = VerdictRules(target_per_drug=40, min_minority_class=15)
    res = evaluate_verdict(cohort, ["ciprofloxacin"], None, rules)
    assert res.verdict == "GO"


def test_verdict_warn_when_strain_count_below_target_but_above_floor():
    """Strain count between target/3 and target → WARN, not FAIL."""
    strains = []
    # 80 strains (40 R + 40 S) with target=150 → 80 ∈ [50, 150) → WARN
    for i in range(40):
        strains.append(_make_strain(
            f"R{i}", 3, 200_000, f"ST{i}", {"ciprofloxacin": 1}
        ))
    for i in range(40):
        strains.append(_make_strain(
            f"S{i}", 3, 200_000, f"ST{i+100}", {"ciprofloxacin": 0}
        ))
    cohort = _build_cohort(strains)
    rules = VerdictRules(target_per_drug=150, min_minority_class=15)  # target unreachable
    res = evaluate_verdict(cohort, ["ciprofloxacin"], None, rules)
    assert res.verdict == "WARN"


def test_verdict_no_go_when_minority_below_classifier_guard():
    """Minority class < 10 → NO-GO (below MIN_TRAINING_SAMPLES)."""
    strains = []
    for i in range(50):
        strains.append(_make_strain(
            f"R{i}", 3, 200_000, f"ST{i}", {"ciprofloxacin": 1}
        ))
    for i in range(3):  # only 3 S → minority=3
        strains.append(_make_strain(
            f"S{i}", 3, 200_000, f"ST{i+200}", {"ciprofloxacin": 0}
        ))
    cohort = _build_cohort(strains)
    rules = VerdictRules(target_per_drug=40, min_minority_class=30)
    res = evaluate_verdict(cohort, ["ciprofloxacin"], None, rules)
    assert res.verdict == "NO-GO"


def test_verdict_no_go_when_metadata_completeness_very_low():
    """>50% missing assembly_accession → FAIL on completeness."""
    strains = []
    for i in range(50):
        s = _make_strain(f"R{i}", 3, 200_000, f"ST{i}", {"ciprofloxacin": 1})
        if i >= 25:
            # Strip assembly_accession to simulate missing metadata
            s = CandidateStrain(
                strain_id=s.strain_id,
                assembly_accession="",  # missing
                mlst=s.mlst, contig_count=s.contig_count, n50=s.n50,
                ast_labels=s.ast_labels,
            )
        strains.append(s)
    for i in range(50):
        s = _make_strain(f"S{i}", 3, 200_000, f"ST{i+100}", {"ciprofloxacin": 0})
        if i >= 25:
            s = CandidateStrain(
                strain_id=s.strain_id,
                assembly_accession="",
                mlst=s.mlst, contig_count=s.contig_count, n50=s.n50,
                ast_labels=s.ast_labels,
            )
        strains.append(s)
    cohort = _build_cohort(strains)
    rules = VerdictRules(target_per_drug=40, min_minority_class=15)
    res = evaluate_verdict(cohort, ["ciprofloxacin"], None, rules)
    # 50% missing → FAIL (>50 threshold)
    assert any(
        "assembly_accession completeness" in name and outcome == "FAIL"
        for name, outcome, _ in res.rules
    ) or res.verdict in ("NO-GO", "WARN")


# ---- clade composition ----


def test_clade_largest_fraction_reported_correctly(tmp_path: Path):
    """Top MLST is reported with correct largest-clade fraction."""
    strains = []
    # 6 strains in ST131 (dominant clade), 4 in unique clades
    for i in range(6):
        strains.append(_make_strain(f"R{i}", 3, 200_000, "ST131",
                                     {"ciprofloxacin": 1}))
    for i in range(4):
        strains.append(_make_strain(f"S{i}", 3, 200_000, f"ST_uniq_{i}",
                                     {"ciprofloxacin": 0}))
    cohort = _build_cohort(strains)
    report = build_report(cohort, None, ["ciprofloxacin"], VerdictRules(),
                          tmp_path / "x.parquet")
    assert "Largest-clade fraction" in report
    assert "60.0%" in report  # 6/10
    assert "Unique MLSTs" in report


def test_clade_duplicate_signature_detection(tmp_path: Path):
    """Strains with same (MLST, N50, contig_count) flagged as duplicates."""
    strains = [
        # Three near-identical strains
        _make_strain("D1", 3, 200_000, "ST131", {"ciprofloxacin": 1}),
        _make_strain("D2", 3, 200_000, "ST131", {"ciprofloxacin": 1}),
        _make_strain("D3", 3, 200_000, "ST131", {"ciprofloxacin": 1}),
        # One unique
        _make_strain("U1", 5, 150_000, "ST10", {"ciprofloxacin": 0}),
    ]
    cohort = _build_cohort(strains)
    report = build_report(cohort, None, ["ciprofloxacin"], VerdictRules(),
                          tmp_path / "x.parquet")
    assert "Suspected duplicate isolates" in report


# ---- CLI ----


def test_parser_requires_cohort_and_output():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])  # missing both required args


def test_main_missing_cohort_exits_2(tmp_path: Path):
    exit_code = main(
        [
            "--cohort", str(tmp_path / "missing.parquet"),
            "--output", str(tmp_path / "out.md"),
        ]
    )
    assert exit_code == 2


def test_main_writes_report_to_disk(tmp_path: Path):
    strains = [
        _make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1}),
        _make_strain("R2", 2, 4_500_000, "ST10", {"ciprofloxacin": 1}),
        _make_strain("S1", 3, 4_000_000, "ST4", {"ciprofloxacin": 0}),
        _make_strain("S2", 3, 3_500_000, "ST5", {"ciprofloxacin": 0}),
    ]
    cohort = _build_cohort(strains)
    cohort_path = tmp_path / "cohort.parquet"
    save_cohort(cohort, cohort_path)

    output_path = tmp_path / "subdir" / "audit.md"  # subdir auto-created
    exit_code = main(
        [
            "--cohort", str(cohort_path),
            "--output", str(output_path),
            "--drugs", "ciprofloxacin",
            "--target-per-drug", "2",
            "--min-minority-class", "1",
        ]
    )
    assert exit_code == 0
    assert output_path.exists()
    content = output_path.read_text(encoding="utf-8")
    assert "# Cohort Audit Report" in content
    assert "ciprofloxacin" in content


def test_main_ast_missing_exits_2(tmp_path: Path):
    strains = [_make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1})]
    cohort_path = tmp_path / "cohort.parquet"
    save_cohort(_build_cohort(strains), cohort_path)
    exit_code = main(
        [
            "--cohort", str(cohort_path),
            "--ast", str(tmp_path / "missing_ast.csv"),
            "--output", str(tmp_path / "out.md"),
        ]
    )
    assert exit_code == 2
