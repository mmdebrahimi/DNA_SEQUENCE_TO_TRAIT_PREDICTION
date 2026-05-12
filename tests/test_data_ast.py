"""Tests for Step 5 — BV-BRC AST phenotype data loader."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from dna_decode.data.ast_data import (
    AST_COLUMNS,
    binarize_susceptibility,
    get_binary_labels,
    get_drug_list,
    load_bvbrc_ast,
)


MOCK_BVBRC_TSV = """\
genome_id\tgenome_name\tantibiotic\tresistant_phenotype\tmeasurement\tmeasurement_unit\tlaboratory_typing_method\ttesting_standard
strain001\tEscherichia coli K-12\tciprofloxacin\tResistant\t8\tmg/L\tBroth microdilution\tCLSI
strain002\tEscherichia coli ST131\tciprofloxacin\tSusceptible\t0.06\tmg/L\tBroth microdilution\tEUCAST
strain003\tEscherichia coli\tceftriaxone\tResistant\t>32\tmg/L\tBroth microdilution\tCLSI
strain004\tEscherichia coli\ttetracycline\tIntermediate\t8\tmg/L\tBroth microdilution\tCLSI
strain005\tEscherichia coli\tciprofloxacin\tResistant\t16\tmg/L\tDisk diffusion\tCLSI
strain006\tSalmonella enterica\tciprofloxacin\tResistant\t8\tmg/L\tBroth microdilution\tCLSI
strain007\tEscherichia coli\tciprofloxacin\tNOT_VALID_LABEL\t2\tmg/L\tBroth microdilution\tCLSI
strain008\tEscherichia coli\tceftriaxone\tSusceptible\t0.125\tmg/L\tBroth microdilution\tCLSI
"""


@pytest.fixture
def ast_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "ast.tsv"
    p.write_text(MOCK_BVBRC_TSV, encoding="utf-8")
    return p


# ---- binarize_susceptibility ----


def test_binarize_resistant():
    assert binarize_susceptibility("Resistant") == 1
    assert binarize_susceptibility("R") == 1


def test_binarize_susceptible():
    assert binarize_susceptibility("Susceptible") == 0
    assert binarize_susceptibility("S") == 0


def test_binarize_intermediate_treated_as_susceptible():
    """v1 binary task treats I as S; MIC regression in Phase 2 uses raw MIC."""
    assert binarize_susceptibility("Intermediate") == 0
    assert binarize_susceptibility("I") == 0


def test_binarize_unknown_raises():
    with pytest.raises(ValueError):
        binarize_susceptibility("UNKNOWN_LABEL")


def test_binarize_handles_whitespace_and_case():
    assert binarize_susceptibility("  r  ") == 1
    assert binarize_susceptibility("\tS\n") == 0


# ---- load_bvbrc_ast ----


def test_load_bvbrc_ast_returns_stable_columns(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    assert list(df.columns) == list(AST_COLUMNS)


def test_load_bvbrc_ast_filters_to_ecoli(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv, organism="Escherichia coli")
    assert all("strain00" in s for s in df["strain_id"])
    # strain006 is Salmonella → excluded
    assert "strain006" not in df["strain_id"].values


def test_load_bvbrc_ast_filters_to_broth_microdilution(ast_tsv: Path):
    """Disk-diffusion rows filtered out per failure-mode #4."""
    df = load_bvbrc_ast(ast_tsv)
    # strain005 used disk diffusion → excluded
    assert "strain005" not in df["strain_id"].values
    # All retained rows have broth_microdilution method
    assert set(df["measurement_method"]) == {"broth_microdilution"}


def test_load_bvbrc_ast_excludes_invalid_labels(ast_tsv: Path):
    """strain007 has 'NOT_VALID_LABEL' → silently dropped."""
    df = load_bvbrc_ast(ast_tsv)
    assert "strain007" not in df["strain_id"].values


def test_load_bvbrc_ast_parses_mic_values(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    # strain001 has '8' MIC → 8.0
    # strain003 has '>32' MIC → 32.0 (lstrip < > =)
    row1 = df[df["strain_id"] == "strain001"].iloc[0]
    row3 = df[df["strain_id"] == "strain003"].iloc[0]
    assert row1["mic_value"] == 8.0
    assert row3["mic_value"] == 32.0


def test_load_bvbrc_ast_lowercase_antibiotic(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    assert "ciprofloxacin" in df["antibiotic"].values
    # Verify no leftover capital-case
    assert not any(a != a.lower() for a in df["antibiotic"])


def test_load_bvbrc_ast_missing_file_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_bvbrc_ast(tmp_path / "missing.tsv")


def test_load_bvbrc_ast_empty_file_returns_empty_table(tmp_path: Path):
    p = tmp_path / "empty.tsv"
    p.write_text("genome_id\tantibiotic\n", encoding="utf-8")  # header only
    df = load_bvbrc_ast(p)
    assert len(df) == 0
    assert list(df.columns) == list(AST_COLUMNS)


# ---- get_drug_list ----


def test_get_drug_list_thresholds_correctly(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    # With min_strains=1 → all drugs present
    drugs = get_drug_list(df, min_strains=1)
    assert "ciprofloxacin" in drugs
    assert "ceftriaxone" in drugs
    assert "tetracycline" in drugs


def test_get_drug_list_high_threshold_excludes(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    # No drug has ≥100 strains in mock data
    assert get_drug_list(df, min_strains=100) == []


def test_get_drug_list_empty_returns_empty():
    import pandas as pd
    empty = pd.DataFrame(columns=list(AST_COLUMNS))
    assert get_drug_list(empty) == []


# ---- get_binary_labels ----


def test_get_binary_labels_returns_strain_to_int(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    cipro_labels = get_binary_labels(df, "ciprofloxacin")
    # strain001 R=1, strain002 S=0
    assert cipro_labels["strain001"] == 1
    assert cipro_labels["strain002"] == 0


def test_get_binary_labels_is_case_insensitive(ast_tsv: Path):
    df = load_bvbrc_ast(ast_tsv)
    labels_upper = get_binary_labels(df, "CIPROFLOXACIN")
    labels_lower = get_binary_labels(df, "ciprofloxacin")
    assert labels_upper == labels_lower
