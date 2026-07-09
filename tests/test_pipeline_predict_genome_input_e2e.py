"""End-to-end integration tests for genome-input `pipeline.py predict`."""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pytest

# xgboost lives in the [ml] extra (not in a default `uv sync`); the fixtures train
# a real classifier, so skip the whole module when it is absent — matching the
# other classifier-dependent test modules.
xgboost = pytest.importorskip("xgboost")


@pytest.fixture
def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def genome_input_fixtures(tmp_path: Path, project_root: Path):
    from datetime import date as _date

    from dna_decode.models.classifiers import train_xgboost_classifier

    fasta_path = tmp_path / "external_sample.fna"
    fasta_path.write_text(
        ">contig1\n"
        "ATGCGTAAACCCTTTGGGCATCATGAATTCGCGCGAGGGCCCAAATTTCCCGGGAGGCGCTAGAGCAGTAGGCATATTGCTCAGCATCGCATGATAAAGTC"
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        encoding="utf-8",
    )

    gff_path = tmp_path / "external_sample.gff3"
    gff_path.write_text(
        "##gff-version 3\n"
        "contig1\t.\tgene\t1\t99\t.\t+\t.\tID=gene-gyrA;Name=gyrA;gene=gyrA\n"
        "contig1\t.\tCDS\t1\t99\t.\t+\t0\tID=gene-gyrA;Parent=gene-gyrA;gene=gyrA;locus_tag=b2231\n"
        "contig1\t.\tgene\t200\t299\t.\t+\t.\tID=gene-parC;Name=parC;gene=parC\n"
        "contig1\t.\tCDS\t200\t299\t.\t+\t0\tID=gene-parC;Parent=gene-parC;gene=parC;locus_tag=b3019\n"
        "contig1\t.\tgene\t350\t449\t.\t+\t.\tID=gene-marR;Name=marR;gene=marR\n"
        "contig1\t.\tCDS\t350\t449\t.\t+\t0\tID=gene-marR;Parent=gene-marR;gene=marR;locus_tag=b1530\n",
        encoding="utf-8",
    )

    rng = np.random.default_rng(seed=42)
    embedding_dim = 128  # mock foundation-model default
    X_train = rng.standard_normal((50, embedding_dim))
    y_train = (X_train[:, 0] + 0.3 * rng.standard_normal(50) > 0).astype(int)
    clf = train_xgboost_classifier(X_train, y_train, drug_name="ciprofloxacin", calibrate=True)

    model_path = tmp_path / "test_model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump(
            {
                "classifier": clf,
                "drug": "ciprofloxacin",
                "model_name": "mock",
                "feature_dim": embedding_dim,
                "auroc_loso": 0.78,
                "auroc_lomo_clade_out": None,
                "strain_id_order": [f"TRAIN_{i}" for i in range(50)],
                "n_strains": 50,
                "training_cohort": "synthetic_test_cohort",
                "trained_on": _date.today().isoformat(),
            },
            f,
        )

    merge_json_path = tmp_path / "test_merge.json"
    merge_json_path.write_text(
        json.dumps(
            {
                "gate_verdict": "SUSPEND_CONDITION_4",
                "per_strain": [
                    {
                        "strain_id": "TRAIN_ONLY_SAMPLE",
                        "noise_class": "OPAQUE_R_no_mechanism",
                        "mic_tier": "HIGH_R",
                        "mechanism_opacity_flag": True,
                        "primary_mechanisms": [],
                        "co_resistance_modifiers": ["efflux"],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    return {
        "config_path": project_root / "config" / "datasources.yaml",
        "fasta_path": fasta_path,
        "gff_path": gff_path,
        "model_path": model_path,
        "merge_json_path": merge_json_path,
        "out_json_path": tmp_path / "prediction.json",
    }


def test_predict_genome_input_e2e_emits_schema_and_cohort_audit_fallback(
    genome_input_fixtures,
):
    from scripts.pipeline import main

    exit_code = main(
        [
            "--config",
            str(genome_input_fixtures["config_path"]),
            "predict",
            "--model-path",
            str(genome_input_fixtures["model_path"]),
            "--genome-fasta",
            str(genome_input_fixtures["fasta_path"]),
            "--annotations",
            str(genome_input_fixtures["gff_path"]),
            "--audit-merge-json",
            str(genome_input_fixtures["merge_json_path"]),
            "--output",
            str(genome_input_fixtures["out_json_path"]),
            "--top-k",
            "5",
        ]
    )
    assert exit_code == 0

    result = json.loads(genome_input_fixtures["out_json_path"].read_text(encoding="utf-8"))
    assert result["strain_id"] == "external_sample"
    assert result["drug"] == "ciprofloxacin"
    assert result["prediction"] in {"R", "S"}
    assert 0.0 <= result["calibrated_probability"] <= 1.0
    assert result["provenance"]["input_mode"] == "genome_input"
    assert result["provenance"]["reporting_mode"] == "canonical_audit_aware"
    assert result["audit_verdict"] is not None
    assert result["audit_verdict"]["cohort_gate_verdict"] == "SUSPEND_CONDITION_4"
    assert result["audit_verdict"]["noise_class"] is None
    assert result["audit_verdict"]["suspend_gate_fired"] is True

    md = genome_input_fixtures["out_json_path"].with_suffix(".md").read_text(encoding="utf-8")
    assert "strain `external_sample`" in md
    assert "SUSPEND gate fired" in md
    assert "external genome is not present in the audit cohort" in md


def test_predict_genome_input_e2e_respects_sample_id_override(genome_input_fixtures):
    from scripts.pipeline import main

    exit_code = main(
        [
            "--config",
            str(genome_input_fixtures["config_path"]),
            "predict",
            "--model-path",
            str(genome_input_fixtures["model_path"]),
            "--genome-fasta",
            str(genome_input_fixtures["fasta_path"]),
            "--annotations",
            str(genome_input_fixtures["gff_path"]),
            "--audit-merge-json",
            str(genome_input_fixtures["merge_json_path"]),
            "--sample-id",
            "EXT_SAMPLE_001",
            "--no-attribution",
            "--output",
            str(genome_input_fixtures["out_json_path"]),
        ]
    )
    assert exit_code == 0

    result = json.loads(genome_input_fixtures["out_json_path"].read_text(encoding="utf-8"))
    assert result["strain_id"] == "EXT_SAMPLE_001"
    assert result["top_k_attribution"] == []


def test_predict_genome_input_e2e_allows_debug_mode_without_audit(genome_input_fixtures):
    from scripts.pipeline import main

    exit_code = main(
        [
            "--config",
            str(genome_input_fixtures["config_path"]),
            "predict",
            "--model-path",
            str(genome_input_fixtures["model_path"]),
            "--genome-fasta",
            str(genome_input_fixtures["fasta_path"]),
            "--annotations",
            str(genome_input_fixtures["gff_path"]),
            "--allow-missing-audit",
            "--no-attribution",
            "--output",
            str(genome_input_fixtures["out_json_path"]),
        ]
    )
    assert exit_code == 0

    result = json.loads(genome_input_fixtures["out_json_path"].read_text(encoding="utf-8"))
    assert result["audit_verdict"] is None
    assert result["provenance"]["reporting_mode"] == "non_canonical_missing_audit"
