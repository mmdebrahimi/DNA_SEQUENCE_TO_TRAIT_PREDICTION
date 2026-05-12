"""Tests for Step 11.5 — Quantization-fidelity check.

Covers the three pure functions that don't touch subprocess / external models:
- compare_top_k_attribution: top-K intersection + Spearman on rank vectors
- build_fidelity_report: aggregate GO/NO-GO verdict from per-strain results
- write_report: markdown serialization

The script's CLI main() path (manifest loading + pd.read_csv + report writing)
is exercised indirectly via build_fidelity_report + write_report; we do NOT
mock subprocess invocations or full classifier runs (Step 11.5 was committed
as a scaffold — real runs require rented A100 + two trained classifiers, per
the module docstring).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts.quantize_fidelity_check import (
    FidelityReport,
    FidelityResult,
    SPEARMAN_GO_THRESHOLD,
    TOP_K_INTERSECTION_GO_THRESHOLD,
    build_fidelity_report,
    compare_top_k_attribution,
    write_report,
)


# ---- compare_top_k_attribution ----


def _make_gene_table(gene_ids: list[str]) -> pd.DataFrame:
    """GeneEffectTable shape: gene_id column + decreasing prediction_delta."""
    return pd.DataFrame(
        {
            "gene_id": gene_ids,
            "prediction_delta": [1.0 - 0.01 * i for i in range(len(gene_ids))],
        }
    )


def test_compare_top_k_identical_tables_full_concordance():
    """Identical top-K → intersection = K and Spearman = 1.0."""
    genes = [f"g{i}" for i in range(20)]
    full = _make_gene_table(genes)
    quant = _make_gene_table(genes)
    intersection, spearman = compare_top_k_attribution(full, quant, top_k=20)
    assert intersection == 20
    assert spearman == pytest.approx(1.0)


def test_compare_top_k_partial_overlap():
    """5 of 10 genes shared → intersection = 5; Spearman drops below 1."""
    full = _make_gene_table([f"g{i}" for i in range(10)])
    quant = _make_gene_table([f"g{i}" for i in range(5, 15)])
    intersection, spearman = compare_top_k_attribution(full, quant, top_k=10)
    assert intersection == 5
    assert spearman is not None
    assert spearman < 1.0


def test_compare_top_k_disjoint_tables():
    """Zero overlap → intersection = 0; Spearman defined but expected low/negative."""
    full = _make_gene_table([f"g{i}" for i in range(10)])
    quant = _make_gene_table([f"h{i}" for i in range(10)])
    intersection, spearman = compare_top_k_attribution(full, quant, top_k=10)
    assert intersection == 0
    # Disjoint top-K with all-absent ranks → low correlation
    assert spearman is None or spearman <= 0.5


def test_compare_top_k_empty_table_returns_zero_none():
    """Empty input short-circuits to (0, None) — guards against degenerate runs."""
    empty = pd.DataFrame({"gene_id": [], "prediction_delta": []})
    nonempty = _make_gene_table(["g1", "g2"])
    assert compare_top_k_attribution(empty, nonempty, top_k=5) == (0, None)
    assert compare_top_k_attribution(nonempty, empty, top_k=5) == (0, None)


def test_compare_top_k_respects_top_k_argument():
    """top_k=3 only inspects first 3 rows of each table."""
    full = _make_gene_table(["a", "b", "c", "d", "e"])
    quant = _make_gene_table(["a", "b", "c", "x", "y"])
    intersection, _ = compare_top_k_attribution(full, quant, top_k=3)
    # First 3 are identical → intersection = 3 regardless of trailing rows
    assert intersection == 3


# ---- build_fidelity_report ----


def _result(strain_id: str, intersection: int, spearman: float | None, top_k: int = 20) -> FidelityResult:
    frac = intersection / top_k
    sp = spearman if spearman is not None else 0.0
    return FidelityResult(
        strain_id=strain_id,
        drug="cipro",
        top_k=top_k,
        intersection_size=intersection,
        spearman=spearman,
        fidelity_score=(frac + max(sp, 0.0)) / 2.0,
    )


def test_build_report_empty_list_is_no_go():
    """No per-strain results → NO-GO + explicit failure reason."""
    report = build_fidelity_report([], top_k=20)
    assert report.go_no_go == "NO-GO"
    assert any("no per-strain" in r.lower() for r in report.failure_reasons)


def test_build_report_strong_concordance_is_go():
    """High intersection + Spearman across strains → GO."""
    results = [_result(f"s{i}", intersection=18, spearman=0.9) for i in range(3)]
    report = build_fidelity_report(results, top_k=20)
    assert report.go_no_go == "GO"
    assert report.failure_reasons == []
    assert report.overall_intersection_mean == pytest.approx(0.9)
    assert report.overall_spearman_mean == pytest.approx(0.9)


def test_build_report_low_intersection_is_no_go_with_reason():
    """Mean intersection below threshold → NO-GO + threshold-specific reason."""
    # intersection_fraction = 4/20 = 0.2, well below TOP_K_INTERSECTION_GO_THRESHOLD (0.6)
    results = [_result(f"s{i}", intersection=4, spearman=0.9) for i in range(3)]
    report = build_fidelity_report(results, top_k=20)
    assert report.go_no_go == "NO-GO"
    assert any("intersection-fraction" in r for r in report.failure_reasons)


def test_build_report_low_spearman_is_no_go_with_reason():
    """Mean Spearman below threshold → NO-GO + threshold-specific reason."""
    # Intersection passes (0.9 ≥ 0.6) but Spearman 0.3 < 0.7
    results = [_result(f"s{i}", intersection=18, spearman=0.3) for i in range(3)]
    report = build_fidelity_report(results, top_k=20)
    assert report.go_no_go == "NO-GO"
    assert any("spearman" in r.lower() for r in report.failure_reasons)


def test_build_report_handles_all_none_spearman():
    """Strains with None Spearman → mean Spearman = 0.0, NO-GO on Spearman threshold."""
    results = [_result(f"s{i}", intersection=18, spearman=None) for i in range(3)]
    report = build_fidelity_report(results, top_k=20)
    assert report.overall_spearman_mean == 0.0
    assert report.go_no_go == "NO-GO"


def test_build_report_thresholds_match_module_constants():
    """Guard against silent threshold drift between code and docstring."""
    assert SPEARMAN_GO_THRESHOLD == 0.7
    assert TOP_K_INTERSECTION_GO_THRESHOLD == 0.6


# ---- write_report ----


def test_write_report_creates_markdown_with_verdict(tmp_path: Path):
    """Markdown contains drug, verdict, mean stats, per-strain table."""
    results = [_result("s1", intersection=18, spearman=0.9)]
    report = build_fidelity_report(results, top_k=20)
    out_path = write_report(report, tmp_path / "subdir" / "report.md", drug="cipro")
    assert out_path.exists()
    text = out_path.read_text(encoding="utf-8")
    assert "# Quantization-Fidelity Check" in text
    assert "**Drug:** cipro" in text
    assert "GO" in text
    assert "| s1 |" in text


def test_write_report_emits_failure_section_on_no_go(tmp_path: Path):
    """NO-GO verdict → Failure reasons + Remediation sections present."""
    results = [_result(f"s{i}", intersection=2, spearman=0.1) for i in range(2)]
    report = build_fidelity_report(results, top_k=20)
    out_path = write_report(report, tmp_path / "report.md", drug="ceftriaxone")
    text = out_path.read_text(encoding="utf-8")
    assert "NO-GO" in text
    assert "## Failure reasons" in text
    assert "Remediation" in text


def test_write_report_none_spearman_renders_em_dash(tmp_path: Path):
    """Per-strain rows with Spearman=None render '—' rather than crash."""
    results = [_result("s1", intersection=15, spearman=None)]
    report = build_fidelity_report(results, top_k=20)
    out_path = write_report(report, tmp_path / "report.md", drug="cipro")
    text = out_path.read_text(encoding="utf-8")
    assert "| s1 | 15/20 | — |" in text
