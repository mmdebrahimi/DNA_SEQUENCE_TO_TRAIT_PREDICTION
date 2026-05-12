"""Tests for Step 13 — matplotlib + TSV viz."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")

from dna_decode.viz.browser import export_attribution_tsv, render_attribution_plot  # noqa: E402


# ---- TSV export ----


def test_export_attribution_tsv_writes_tab_separated(tmp_path: Path):
    df = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "locus_tag": ["TAG_001", "TAG_002"],
            "prediction_delta": [0.5, -0.3],
            "baseline_probability": [0.7, 0.7],
            "knockout_probability": [0.2, 1.0],
        }
    )
    out, _ = export_attribution_tsv(df, tmp_path / "ge.tsv")
    assert out.exists()
    content = out.read_text(encoding="utf-8")
    # Tab-separated; header line has 5 tab-delimited columns
    first_line = content.splitlines()[0]
    assert first_line.count("\t") == 4


def test_export_attribution_tsv_optional_position_panel(tmp_path: Path):
    gene = pd.DataFrame({"gene_id": ["g1"], "prediction_delta": [0.5]})
    pos = pd.DataFrame(
        {
            "gene_id": ["g1", "g1"],
            "position": [10, 20],
            "ref_base": ["A", "T"],
            "alt_base": ["C", "G"],
            "prediction_delta": [0.1, 0.2],
        }
    )
    gene_path, pos_path = export_attribution_tsv(
        gene,
        tmp_path / "gene.tsv",
        position_effects=pos,
        position_output_path=tmp_path / "pos.tsv",
    )
    assert gene_path.exists() and pos_path.exists()


def test_export_attribution_tsv_position_without_path_raises(tmp_path: Path):
    gene = pd.DataFrame({"gene_id": ["g1"], "prediction_delta": [0.5]})
    pos = pd.DataFrame(
        {"gene_id": ["g1"], "position": [10], "ref_base": ["A"], "alt_base": ["C"], "prediction_delta": [0.1]}
    )
    with pytest.raises(ValueError, match="position_output_path"):
        export_attribution_tsv(
            gene, tmp_path / "gene.tsv", position_effects=pos, position_output_path=None
        )


# ---- matplotlib rendering ----


def test_render_attribution_plot_gene_only(tmp_path: Path):
    df = pd.DataFrame(
        {
            "gene_id": [f"g{i}" for i in range(5)],
            "locus_tag": [f"TAG_{i:03d}" for i in range(5)],
            "prediction_delta": [0.5, -0.4, 0.3, -0.2, 0.1],
            "baseline_probability": [0.7] * 5,
            "knockout_probability": [0.2, 1.0, 0.4, 0.9, 0.6],
        }
    )
    out = render_attribution_plot(df, tmp_path / "viz.png", top_k=5, drug_name="cipro")
    assert out.exists()
    # PNG file starts with the standard magic bytes
    head = out.read_bytes()[:8]
    assert head.startswith(b"\x89PNG\r\n\x1a\n")


def test_render_attribution_plot_with_position_panel(tmp_path: Path):
    gene = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "locus_tag": ["TAG_001", "TAG_002"],
            "prediction_delta": [0.5, 0.2],
            "baseline_probability": [0.7, 0.7],
            "knockout_probability": [0.2, 0.5],
        }
    )
    pos = pd.DataFrame(
        {
            "gene_id": ["g1"] * 6,
            "position": [10, 10, 50, 50, 100, 100],
            "ref_base": ["A"] * 6,
            "alt_base": ["C", "G", "A", "T", "A", "C"],
            "prediction_delta": [0.5, 0.6, 0.01, 0.02, 0.4, 0.3],
        }
    )
    out = render_attribution_plot(
        gene, tmp_path / "viz_pos.png", top_k=2, drug_name="cipro", position_effects=pos
    )
    assert out.exists()


def test_render_attribution_plot_empty_position_panel_falls_back(tmp_path: Path):
    """Empty position_effects → top panel still renders; no crash."""
    gene = pd.DataFrame(
        {"gene_id": ["g1"], "locus_tag": [""], "prediction_delta": [0.5],
         "baseline_probability": [0.7], "knockout_probability": [0.2]}
    )
    empty_pos = pd.DataFrame(
        columns=["gene_id", "position", "ref_base", "alt_base", "prediction_delta"]
    )
    out = render_attribution_plot(
        gene, tmp_path / "viz_empty.png", position_effects=empty_pos
    )
    assert out.exists()
