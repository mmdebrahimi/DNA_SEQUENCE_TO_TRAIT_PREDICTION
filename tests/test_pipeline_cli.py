"""Tests for Step 14 — Single scripts/pipeline.py CLI dispatcher."""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.data.cohort import CandidateStrain, StrainCohort, save_cohort

xgboost = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")


def test_build_parser_has_all_four_subcommands():
    """Parser exposes ingest / train / predict / attribute subcommands."""
    from scripts.pipeline import build_parser

    parser = build_parser()
    # Try parsing minimal valid invocations for each subcommand
    parser.parse_args(["ingest", "--drugs", "cipro", "--ast-tsv", "fake.tsv"])
    parser.parse_args(["train", "--drug", "cipro"])
    parser.parse_args(["predict", "--model-path", "fake.pkl"])
    parser.parse_args(["attribute", "--model-path", "fake.pkl", "--strain-id", "s001"])


def test_pipeline_missing_config_exits_2(tmp_path: Path, monkeypatch):
    """Bad --config path → exit 2."""
    monkeypatch.chdir(tmp_path)
    from scripts.pipeline import main

    exit_code = main(
        [
            "--config",
            str(tmp_path / "missing.yaml"),
            "ingest",
            "--drugs",
            "cipro",
            "--ast-tsv",
            "fake.tsv",
        ]
    )
    assert exit_code == 2


def test_ingest_without_ast_tsv_arg_argparse_fails(monkeypatch, project_root: Path):
    """Argparse enforces --ast-tsv required for ingest."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    # Argparse exits 2 (its own convention) on missing required args
    with pytest.raises(SystemExit) as exc_info:
        main(["ingest", "--drugs", "cipro"])
    assert exc_info.value.code == 2


def test_train_missing_cohort_exits_2(tmp_path: Path, project_root: Path, monkeypatch):
    """train fails clean when cohort.parquet doesn't exist."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    exit_code = main(
        [
            "train",
            "--drug",
            "cipro",
            "--cohort",
            str(tmp_path / "missing.parquet"),
        ]
    )
    assert exit_code == 2


def test_predict_missing_model_exits_2(tmp_path: Path, project_root: Path, monkeypatch):
    """predict fails clean when --model-path doesn't exist."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    exit_code = main(
        [
            "predict",
            "--model-path",
            str(tmp_path / "missing.pkl"),
            "--strain-id",
            "s001",
        ]
    )
    assert exit_code == 2


def test_predict_missing_audit_merge_exits_2(tmp_path: Path, project_root: Path, monkeypatch):
    """Canonical predict should refuse audit-free invocation before output generation."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    model_path = tmp_path / "model.pkl"
    model_path.write_bytes(b"not-a-real-model")

    exit_code = main(
        [
            "predict",
            "--model-path",
            str(model_path),
            "--strain-id",
            "s001",
        ]
    )
    assert exit_code == 2


def test_predict_rejects_both_input_modes(tmp_path: Path, project_root: Path, monkeypatch):
    """predict requires exactly one of --strain-id or --genome-fasta."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    model_path = tmp_path / "model.pkl"
    model_path.write_bytes(b"not-a-real-model")
    fasta = tmp_path / "sample.fna"
    fasta.write_text(">x\nATGC\n", encoding="utf-8")

    exit_code = main(
        [
            "predict",
            "--model-path",
            str(model_path),
            "--strain-id",
            "s001",
            "--genome-fasta",
            str(fasta),
        ]
    )
    assert exit_code == 2


def test_predict_genome_input_requires_annotations(
    tmp_path: Path, project_root: Path, monkeypatch
):
    """Genome-input mode fails closed without --annotations."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    model_path = tmp_path / "model.pkl"
    model_path.write_bytes(b"not-a-real-model")
    fasta = tmp_path / "sample.fna"
    fasta.write_text(">x\nATGC\n", encoding="utf-8")

    exit_code = main(
        [
            "predict",
            "--model-path",
            str(model_path),
            "--genome-fasta",
            str(fasta),
        ]
    )
    assert exit_code == 2


def test_attribute_missing_model_exits_2(tmp_path: Path, project_root: Path, monkeypatch):
    """attribute fails clean when --model-path doesn't exist."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    exit_code = main(
        [
            "attribute",
            "--model-path",
            str(tmp_path / "missing.pkl"),
            "--strain-id",
            "s001",
        ]
    )
    assert exit_code == 2


# ---- Option B (Phase 2): --assembly-metadata-csv ----


def _write_tiny_ast(tmp_path: Path) -> Path:
    """Write a 2-strain AST CSV that exercises the load_bvbrc_ast loader."""
    ast = (
        "Genome ID,Genome Name,Antibiotic,Resistant Phenotype,Measurement,"
        "Measurement Unit,Laboratory Typing Method,Testing Standard\n"
        "562.1,Escherichia coli A,ciprofloxacin,Resistant,8,mg/L,Broth dilution,CLSI\n"
        "562.2,Escherichia coli B,ciprofloxacin,Susceptible,0.06,mg/L,Broth dilution,EUCAST\n"
    )
    p = tmp_path / "ast.csv"
    p.write_text(ast, encoding="utf-8")
    return p


def _write_tiny_genome_csv(tmp_path: Path) -> Path:
    """Write a 1-strain genome CSV that covers half the AST cohort."""
    g = (
        "Genome ID,Genome Name,Species,MLST,Assembly Accession,Contigs,"
        "Contig N50,Collection Year,Isolation Country\n"
        "562.1,Escherichia coli A,Escherichia coli,ST1,GCF_A.1,10,200000,2020,USA\n"
    )
    p = tmp_path / "genome.csv"
    p.write_text(g, encoding="utf-8")
    return p


def test_ingest_accepts_assembly_metadata_csv_flag(
    tmp_path: Path, project_root: Path, monkeypatch, capsys
):
    """Happy path: --assembly-metadata-csv loads + emits coverage-log line."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    ast_path = _write_tiny_ast(tmp_path)
    genome_csv = _write_tiny_genome_csv(tmp_path)

    exit_code = main(
        [
            "ingest",
            "--drugs", "ciprofloxacin",
            "--ast-tsv", str(ast_path),
            "--assembly-metadata-csv", str(genome_csv),
            "--target-per-drug", "1",
            "--intersection-target", "0",
            "--cohort-out", str(tmp_path / "cohort.parquet"),
        ]
    )
    # Cohort construction may fail (not enough strains for Phase 1 defaults) but
    # the assembly-metadata wire should run cleanly to coverage-log emission.
    out = capsys.readouterr().out
    assert "loaded assembly metadata for 1 strains" in out
    assert "from CSV" in out
    assert "assembly_meta covers 1 / 2 AST strain_ids (50.0%)" in out
    # exit code is whatever build_cohort returns; the wire validation is in stdout


def test_ingest_rejects_both_metadata_flags_simultaneously(
    tmp_path: Path, project_root: Path, monkeypatch
):
    """argparse mutex group: passing both flags exits 2."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    ast_path = _write_tiny_ast(tmp_path)
    yaml_path = tmp_path / "meta.yaml"
    yaml_path.write_text("{}", encoding="utf-8")
    csv_path = _write_tiny_genome_csv(tmp_path)

    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "ingest",
                "--drugs", "ciprofloxacin",
                "--ast-tsv", str(ast_path),
                "--assembly-metadata", str(yaml_path),
                "--assembly-metadata-csv", str(csv_path),
            ]
        )
    assert exc_info.value.code == 2


def test_ingest_assembly_metadata_csv_missing_exits_2(
    tmp_path: Path, project_root: Path, monkeypatch
):
    """File-not-found on --assembly-metadata-csv → clean exit 2."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    ast_path = _write_tiny_ast(tmp_path)

    exit_code = main(
        [
            "ingest",
            "--drugs", "ciprofloxacin",
            "--ast-tsv", str(ast_path),
            "--assembly-metadata-csv", str(tmp_path / "missing.csv"),
        ]
    )
    assert exit_code == 2


def test_train_duplicate_accessions_exit_2(tmp_path: Path, project_root: Path, monkeypatch, capsys):
    """Explicit strain_id CV refuses silent leakage when the same assembly appears twice."""
    monkeypatch.chdir(project_root)
    from scripts.pipeline import main

    cohort = StrainCohort(
        strains=[
            CandidateStrain(
                strain_id="s1",
                assembly_accession="GCA_DUP.1",
                ast_labels={"ciprofloxacin": 1},
            ),
            CandidateStrain(
                strain_id="s2",
                assembly_accession="GCA_DUP.1",
                ast_labels={"ciprofloxacin": 0},
            ),
        ],
        per_drug_strain_ids={"ciprofloxacin": ["s1", "s2"]},
        three_drug_intersection=[],
    )
    cohort_path = save_cohort(cohort, tmp_path / "cohort.parquet")
    cache_path = tmp_path / "cache.h5"
    cache_path.write_bytes(b"")

    exit_code = main(
        [
            "train",
            "--drug",
            "ciprofloxacin",
            "--cohort",
            str(cohort_path),
            "--cache",
            str(cache_path),
            "--model",
            "nucleotide_transformer",
            "--cv-grouping",
            "strain_id",
        ]
    )
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "duplicate assembly_accession values detected" in err
    assert "GCA_DUP.1" in err
