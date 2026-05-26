"""End-to-end integration test for `pipeline.py predict` (v0 schema).

Validates the v0 success criterion #1 ("Functional, HARD: predict runs
end-to-end without crashing and emits JSON + markdown outputs") using
synthetic fixtures — no real Databricks cache or trained model needed.

Fixtures synthesized:
- A trained XGBoost pickle with v0 provenance fields (training_cohort,
  trained_on, auroc_loso, model_name, drug)
- An HDF5 embedding cache with 3 genes for 1 strain (mean-pooled features)
- A GFF3 annotation file with 3 gene rows
- A merge-gate JSON with the strain's noise class + cohort gate verdict

Skipped only if the synthetic foundation-model config can't be wired
(catalog mismatch); else runs the full predict path.
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def synthetic_fixtures(tmp_path: Path, project_root: Path):
    """Build all the synthetic inputs needed for cmd_predict to run.

    Returns a dict with paths + the strain_id under test.
    """
    from dna_decode.models.cache import EmbeddingCache
    from dna_decode.models.classifiers import train_xgboost_classifier
    from datetime import date as _date

    strain_id = "TEST_STRAIN_001"
    embedding_dim = 8
    gene_ids = ["gene-gyrA", "gene-parC", "gene-marR"]

    # --- 1. Synthetic embedding cache ---
    cache_path = tmp_path / "synthetic_cache.h5"
    cache = EmbeddingCache(
        cache_path,
        model_name="mock_nt",
        model_version="mock_v1",
        embedding_dim=embedding_dim,
    )
    rng = np.random.default_rng(seed=42)
    for gene_id in gene_ids:
        # Per-gene embedding stored as 1D (the foundation-model + cache layer
        # mean-pools windows before cache.put — see EmbeddingCache.populate).
        embedding = rng.standard_normal(embedding_dim).astype(np.float32)
        cache.put(strain_id, gene_id, embedding)

    # --- 2. Synthetic trained classifier pickle (v0 provenance) ---
    # Mean-pooled feature per strain = embedding_dim-length vector
    n_train = 50
    X_train = rng.standard_normal((n_train, embedding_dim))
    y_train = (X_train[:, 0] + 0.3 * rng.standard_normal(n_train) > 0).astype(int)
    clf = train_xgboost_classifier(X_train, y_train, drug_name="ciprofloxacin", calibrate=True)

    model_path = tmp_path / "test_model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump({
            "classifier": clf,
            "drug": "ciprofloxacin",
            "model_name": "mock_nt",
            "feature_dim": embedding_dim,
            "auroc_loso": 0.78,
            "auroc_lomo_clade_out": None,
            "strain_id_order": [f"TRAIN_{i}" for i in range(n_train)],
            "n_strains": n_train,
            "training_cohort": "synthetic_test_cohort",
            "trained_on": _date.today().isoformat(),
        }, f)

    # --- 3. Synthetic GFF3 annotation ---
    gff_path = tmp_path / "test_annotations.gff3"
    gff_path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene-gyrA;Name=gyrA;gene=gyrA\n"
        "chr1\t.\tCDS\t1\t100\t.\t+\t0\tID=cds-gyrA;Parent=gene-gyrA;gene=gyrA;locus_tag=b2231\n"
        "chr1\t.\tgene\t200\t300\t.\t+\t.\tID=gene-parC;Name=parC;gene=parC\n"
        "chr1\t.\tCDS\t200\t300\t.\t+\t0\tID=cds-parC;Parent=gene-parC;gene=parC;locus_tag=b3019\n"
        "chr1\t.\tgene\t400\t500\t.\t+\t.\tID=gene-marR;Name=marR;gene=marR\n"
        "chr1\t.\tCDS\t400\t500\t.\t+\t0\tID=cds-marR;Parent=gene-marR;gene=marR;locus_tag=b1530\n",
        encoding="utf-8",
    )

    # --- 4. Synthetic merge-gate JSON with the strain's audit verdict ---
    merge_json_path = tmp_path / "test_merge.json"
    merge_json_path.write_text(json.dumps({
        "gate_verdict": "SUSPEND_CONDITION_4",
        "recommended_next_step": "SUSPEND_CONDITION_4",
        "per_strain": [{
            "strain_id": strain_id,
            "noise_class": "OPAQUE_R_no_mechanism",
            "mic_tier": "HIGH_R",
            "mechanism_opacity_flag": True,
            "primary_mechanisms": [],
            "co_resistance_modifiers": ["efflux"],
        }],
    }), encoding="utf-8")

    # --- 5. Outputs ---
    out_json_path = tmp_path / "prediction.json"

    return {
        "strain_id": strain_id,
        "model_path": model_path,
        "cache_path": cache_path,
        "gff_path": gff_path,
        "merge_json_path": merge_json_path,
        "out_json_path": out_json_path,
        "config_path": project_root / "config" / "datasources.yaml",
    }


def _synthetic_config_with_mock_nt(real_config_path: Path, tmp_path: Path) -> Path:
    """Append a mock_nt foundation-model entry to a config copy."""
    import yaml
    cfg = yaml.safe_load(real_config_path.read_text(encoding="utf-8"))
    cfg.setdefault("foundation_models", {})["mock_nt"] = {
        "huggingface_id": "mock_v1",
        "embedding_dim": 8,
        "max_context": 1024,
    }
    copy_path = tmp_path / "test_config.yaml"
    copy_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return copy_path


def test_predict_e2e_emits_v0_json_schema(synthetic_fixtures, tmp_path: Path):
    """v0 criterion #1 (HARD): predict runs end-to-end + JSON has all v0 fields."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    exit_code = main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--annotations", str(synthetic_fixtures["gff_path"]),
        "--audit-merge-json", str(synthetic_fixtures["merge_json_path"]),
        "--output", str(synthetic_fixtures["out_json_path"]),
        "--top-k", "5",
    ])
    assert exit_code == 0, "predict should complete cleanly on synthetic fixtures"

    # JSON output exists + has v0 schema
    assert synthetic_fixtures["out_json_path"].exists()
    result = json.loads(synthetic_fixtures["out_json_path"].read_text(encoding="utf-8"))
    for required_field in (
        "strain_id", "drug", "prediction", "calibrated_probability",
        "confidence_tier", "top_k_attribution", "audit_verdict", "provenance",
    ):
        assert required_field in result, f"missing v0 schema field: {required_field}"

    # Strain + drug round-trip
    assert result["strain_id"] == synthetic_fixtures["strain_id"]
    assert result["drug"] == "ciprofloxacin"
    assert result["prediction"] in {"R", "S"}
    assert 0.0 <= result["calibrated_probability"] <= 1.0
    assert result["confidence_tier"] in {"HIGH", "MEDIUM", "LOW"}


def test_predict_e2e_emits_markdown_sidecar(synthetic_fixtures, tmp_path: Path):
    """v0 criterion #1 (HARD): markdown sidecar generated alongside JSON."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    exit_code = main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--annotations", str(synthetic_fixtures["gff_path"]),
        "--audit-merge-json", str(synthetic_fixtures["merge_json_path"]),
        "--output", str(synthetic_fixtures["out_json_path"]),
    ])
    assert exit_code == 0

    # Markdown sidecar = <output>.md
    md_path = synthetic_fixtures["out_json_path"].with_suffix(".md")
    assert md_path.exists()
    md = md_path.read_text(encoding="utf-8")
    assert "strain `" + synthetic_fixtures["strain_id"] + "`" in md
    assert "**Prediction:**" in md
    assert "Top-K gene attribution" in md
    assert "Provenance" in md
    assert "Not a clinical decision support tool" in md  # honest-output disclaimer


def test_predict_e2e_propagates_suspend_verdict(synthetic_fixtures, tmp_path: Path):
    """v0 criterion #4 (HARD): audit-verdict SUSPEND status propagates to output."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--annotations", str(synthetic_fixtures["gff_path"]),
        "--audit-merge-json", str(synthetic_fixtures["merge_json_path"]),
        "--output", str(synthetic_fixtures["out_json_path"]),
    ])

    result = json.loads(synthetic_fixtures["out_json_path"].read_text(encoding="utf-8"))
    av = result["audit_verdict"]
    assert av is not None
    assert av["cohort_gate_verdict"] == "SUSPEND_CONDITION_4"
    assert av["suspend_gate_fired"] is True
    assert "informational only" in av["verdict_explanation"]

    # Markdown should also surface SUSPEND
    md = synthetic_fixtures["out_json_path"].with_suffix(".md").read_text(encoding="utf-8")
    assert "SUSPEND gate fired" in md


def test_predict_e2e_requires_audit_merge_json_by_default(synthetic_fixtures, tmp_path: Path):
    """Canonical predict should fail closed when audit framing is omitted."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    exit_code = main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--annotations", str(synthetic_fixtures["gff_path"]),
        "--output", str(synthetic_fixtures["out_json_path"]),
    ])
    assert exit_code == 2
    assert not synthetic_fixtures["out_json_path"].exists()


def test_predict_e2e_allows_missing_audit_in_debug_mode(synthetic_fixtures, tmp_path: Path):
    """Non-canonical internal/debug mode can still emit audit-free output explicitly."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    exit_code = main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--annotations", str(synthetic_fixtures["gff_path"]),
        "--allow-missing-audit",
        "--output", str(synthetic_fixtures["out_json_path"]),
    ])
    assert exit_code == 0

    result = json.loads(synthetic_fixtures["out_json_path"].read_text(encoding="utf-8"))
    assert result["audit_verdict"] is None
    assert result["provenance"]["reporting_mode"] == "non_canonical_missing_audit"
    md = synthetic_fixtures["out_json_path"].with_suffix(".md").read_text(encoding="utf-8")
    assert "Non-canonical internal/debug run." in md


def test_predict_e2e_no_attribution_skips_top_k(synthetic_fixtures, tmp_path: Path):
    """--no-attribution flag should leave top_k_attribution empty."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--audit-merge-json", str(synthetic_fixtures["merge_json_path"]),
        "--no-attribution",
        "--output", str(synthetic_fixtures["out_json_path"]),
    ])

    result = json.loads(synthetic_fixtures["out_json_path"].read_text(encoding="utf-8"))
    assert result["top_k_attribution"] == []


def test_predict_e2e_provenance_block_populated(synthetic_fixtures, tmp_path: Path):
    """Provenance fields from the train pickle should propagate to predict output."""
    from scripts.pipeline import main

    cfg_path = _synthetic_config_with_mock_nt(
        synthetic_fixtures["config_path"], tmp_path
    )

    main([
        "--config", str(cfg_path),
        "predict",
        "--model-path", str(synthetic_fixtures["model_path"]),
        "--strain-id", synthetic_fixtures["strain_id"],
        "--cache", str(synthetic_fixtures["cache_path"]),
        "--audit-merge-json", str(synthetic_fixtures["merge_json_path"]),
        "--output", str(synthetic_fixtures["out_json_path"]),
        "--no-attribution",
    ])

    prov = json.loads(synthetic_fixtures["out_json_path"].read_text(encoding="utf-8"))["provenance"]
    assert prov["training_cohort"] == "synthetic_test_cohort"
    assert prov["reporting_mode"] == "canonical_audit_aware"
    assert prov["loso_auroc"] == 0.78
    assert prov["trained_on"]  # date string from pickle
    assert "XGBoost" in prov["model"]
