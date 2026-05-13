"""Tests for scripts/audit_cohort.py — Phase 2.5 audit report generator."""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.data.cohort import CandidateStrain, StrainCohort, save_cohort
from scripts.audit_cohort import (
    PHASE1_DEFAULT_RULES,
    VerdictRules,
    build_parser,
    build_report,
    evaluate_verdict,
    main,
    non_default_threshold_fields,
    relaxed_flags_warning,
    stdout_warning_lines,
    thresholds_block,
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


# ---- Phase 2.5 calibration discipline (post-/brainstorm patches D6+D7+D4) ----


def test_default_thresholds_warn_on_undersized_cohort():
    """Regression lock: a 50R/50S cohort under Phase 1 defaults must yield WARN.

    Origin: commit b646fc9 shipped a misleading 'GO' on the 67-strain Gate B cohort
    because no test exercised default semantics — only custom-threshold cases.
    This test pins the documented Phase 1 thresholds so silent drift fires."""
    strains = []
    for i in range(50):
        strains.append(_make_strain(
            f"R{i}", 3, 200_000, f"ST{i}", {"ciprofloxacin": 1}
        ))
    for i in range(50):
        strains.append(_make_strain(
            f"S{i}", 3, 200_000, f"ST{i+100}", {"ciprofloxacin": 0}
        ))
    cohort = _build_cohort(strains)
    res = evaluate_verdict(cohort, ["ciprofloxacin"], None, VerdictRules())
    assert res.verdict == "WARN"
    # Specifically: strain count warns because 50 < 150
    assert any(
        "strain count" in name and outcome == "WARN" and "150" in detail
        for name, outcome, detail in res.rules
    )


def test_thresholds_block_renders_all_six_fields(tmp_path: Path):
    block = thresholds_block(VerdictRules())
    rendered = "\n".join(block)
    assert "**Thresholds applied:**" in rendered
    for field_name in (
        "target_per_drug: 150",
        "min_minority_class: 30",
        "max_pct_missing_metadata: 20.0",
        "min_pct_broth_microdilution: 80.0",
        "n50_min: 50000",
        "contig_count_max: 500",
    ):
        assert field_name in rendered


def test_thresholds_block_reflects_overrides(tmp_path: Path):
    custom = VerdictRules(target_per_drug=12, min_minority_class=6)
    rendered = "\n".join(thresholds_block(custom))
    assert "target_per_drug: 12" in rendered
    assert "min_minority_class: 6" in rendered


def test_non_default_threshold_fields_empty_for_defaults():
    assert non_default_threshold_fields(VerdictRules()) == []


def test_non_default_threshold_fields_detects_target_per_drug_relax():
    relaxed = VerdictRules(target_per_drug=50)
    fields = non_default_threshold_fields(relaxed)
    names = {f[0] for f in fields}
    assert "target_per_drug" in names
    # Direction = "lower" (50 < 150)
    assert any(f[3] == "lower" for f in fields if f[0] == "target_per_drug")


def test_non_default_threshold_fields_detects_all_six_directions():
    """Cover all six knobs in their permissive direction."""
    permissive = VerdictRules(
        target_per_drug=10,           # lower
        min_minority_class=2,         # lower
        max_pct_missing_metadata=50,  # higher
        min_pct_broth_microdilution=10,  # lower
        n50_min=1000,                 # lower
        contig_count_max=10000,       # higher
    )
    fields = non_default_threshold_fields(permissive)
    names = {f[0] for f in fields}
    assert names == {
        "target_per_drug",
        "min_minority_class",
        "max_pct_missing_metadata",
        "min_pct_broth_microdilution",
        "n50_min",
        "contig_count_max",
    }


def test_non_default_threshold_fields_ignores_tightened_thresholds():
    """Only permissive deviations trigger the warning. Tighter is fine."""
    tightened = VerdictRules(target_per_drug=200, min_minority_class=50)
    assert non_default_threshold_fields(tightened) == []


def test_relaxed_flags_warning_empty_on_defaults():
    assert relaxed_flags_warning(VerdictRules()) == []


def test_relaxed_flags_warning_lists_deviated_fields():
    relaxed = VerdictRules(target_per_drug=50, max_pct_missing_metadata=40.0)
    banner = "\n".join(relaxed_flags_warning(relaxed))
    assert "WARNING" in banner
    assert "NOT COMPARABLE TO PHASE 1 GATE" in banner
    assert "target_per_drug" in banner and "50" in banner
    assert "max_pct_missing_metadata" in banner and "40" in banner


def test_stdout_warning_lines_match_banner_coverage():
    """stdout warning fires on the same conditions as the report banner."""
    relaxed = VerdictRules(target_per_drug=50)
    lines = stdout_warning_lines(relaxed)
    assert any("WARNING: relaxed thresholds" in line for line in lines)
    assert any("target_per_drug=50" in line for line in lines)
    # Defaults → no stdout warning
    assert stdout_warning_lines(VerdictRules()) == []


def test_build_report_includes_thresholds_block(tmp_path: Path):
    strains = [
        _make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1}),
        _make_strain("S1", 2, 4_500_000, "ST10", {"ciprofloxacin": 0}),
    ]
    cohort = _build_cohort(strains)
    report = build_report(cohort, None, ["ciprofloxacin"], VerdictRules(), tmp_path / "x.parquet")
    assert "**Thresholds applied:**" in report
    assert "target_per_drug: 150" in report
    assert "min_minority_class: 30" in report


def test_build_report_omits_warning_banner_on_defaults(tmp_path: Path):
    strains = [_make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1})]
    cohort = _build_cohort(strains)
    report = build_report(cohort, None, ["ciprofloxacin"], VerdictRules(), tmp_path / "x.parquet")
    assert "WARNING: NON-DEFAULT THRESHOLDS" not in report


def test_build_report_includes_warning_banner_on_relaxed_thresholds(tmp_path: Path):
    strains = [_make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1})]
    cohort = _build_cohort(strains)
    relaxed = VerdictRules(target_per_drug=50, min_minority_class=10)
    report = build_report(cohort, None, ["ciprofloxacin"], relaxed, tmp_path / "x.parquet")
    assert "WARNING: NON-DEFAULT THRESHOLDS" in report
    assert "target_per_drug" in report and "50" in report
    assert "min_minority_class" in report and "10" in report


def test_main_stdout_includes_warning_on_relaxed_flags(tmp_path: Path, capsys):
    strains = [
        _make_strain("R1", 2, 5_000_000, "ST131", {"ciprofloxacin": 1}),
        _make_strain("S1", 2, 4_500_000, "ST10", {"ciprofloxacin": 0}),
    ]
    cohort = _build_cohort(strains)
    cohort_path = tmp_path / "cohort.parquet"
    save_cohort(cohort, cohort_path)

    exit_code = main(
        [
            "--cohort", str(cohort_path),
            "--output", str(tmp_path / "audit.md"),
            "--drugs", "ciprofloxacin",
            "--target-per-drug", "1",
            "--min-minority-class", "1",
        ]
    )
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "WARNING: relaxed thresholds" in out
    # Verdict echo also present
    assert "[audit_cohort] verdict:" in out


def test_main_stdout_no_warning_on_defaults(tmp_path: Path, capsys):
    """Default flags → no WARNING stdout line."""
    # Build a 50R/50S cohort so verdict is WARN under defaults but the cohort itself
    # has enough strains that the CLI runs without --target-per-drug overrides
    strains = []
    for i in range(50):
        strains.append(_make_strain(f"R{i}", 3, 200_000, f"ST{i}", {"ciprofloxacin": 1}))
    for i in range(50):
        strains.append(_make_strain(f"S{i}", 3, 200_000, f"ST{i+100}", {"ciprofloxacin": 0}))
    cohort = _build_cohort(strains)
    cohort_path = tmp_path / "cohort.parquet"
    save_cohort(cohort, cohort_path)

    exit_code = main(
        [
            "--cohort", str(cohort_path),
            "--output", str(tmp_path / "audit.md"),
            "--drugs", "ciprofloxacin",
            # No --target-per-drug etc. → defaults apply
        ]
    )
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "WARNING: relaxed thresholds" not in out
    # Verdict still WARN (50 < 150), but no relaxed-threshold warning
    assert "verdict: WARN" in out


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
