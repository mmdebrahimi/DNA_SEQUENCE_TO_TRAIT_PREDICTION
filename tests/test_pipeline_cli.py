"""Tests for Step 14 — Single scripts/pipeline.py CLI dispatcher."""
from __future__ import annotations

from pathlib import Path

import pytest

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
