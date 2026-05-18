"""Tests for scripts/cipro_curated_baseline.py — module-level contracts.

The full LOSO + LR + XGB run is ML orchestration (skipped). Tests pin the
module-level constants + ABLATION_FEATURE_SETS contract + the
_load_amrfinder_features JSON-parsing logic that drives the 2-layer verdict.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.cipro_curated_baseline import (
    ABLATION_FEATURE_SETS,
    AMENDED_MECHANISM_ONLY_GATE_AUROC,
    AMENDED_NO_POINT_GATE_AUROC,
    BASELINE_ABSOLUTE_GATE_AUROC,
    BASELINE_COMPARATIVE_GATE_PP,
    STAGE1B_KMER_XGB_AUROC,
    STAGE1B_NT_LR_AUROC,
    STAGE1B_NT_XGB_AUROC,
    _load_amrfinder_features,
)


# ---- Constants contract (frozen reference points) ---------------------------


def test_baseline_absolute_gate_locked_at_0_80():
    # Stage 2 burst threshold; do not soften
    assert BASELINE_ABSOLUTE_GATE_AUROC == 0.80


def test_baseline_comparative_gate_locked_at_10_pp():
    assert BASELINE_COMPARATIVE_GATE_PP == 10.0


def test_stage1b_baselines_pinned():
    # These reference values frame the amended condition 4 gate; freezing them
    # so any future Stage 1b re-run that changes the AUROC must explicitly
    # update this constant + re-derive AMENDED_NO_POINT_GATE_AUROC
    assert STAGE1B_NT_LR_AUROC == 0.673
    assert STAGE1B_NT_XGB_AUROC == 0.615
    assert STAGE1B_KMER_XGB_AUROC == 0.648


def test_amended_no_point_gate_consistent_with_comparative_framing():
    # AMENDED_NO_POINT_GATE = max(0.75, STAGE1B_NT_LR_AUROC + 0.10)
    # = max(0.75, 0.773) = 0.773
    expected = max(0.75, STAGE1B_NT_LR_AUROC + 0.10)
    assert AMENDED_NO_POINT_GATE_AUROC == expected
    assert AMENDED_NO_POINT_GATE_AUROC >= 0.773 - 1e-9


def test_amended_mechanism_only_gate_equals_absolute_baseline():
    # mechanism-only path uses the absolute gate (does AMRFinder panel alone
    # match the Stage 2 burst threshold?)
    assert AMENDED_MECHANISM_ONLY_GATE_AUROC == BASELINE_ABSOLUTE_GATE_AUROC


# ---- ABLATION_FEATURE_SETS contract -----------------------------------------


def test_ablation_feature_sets_includes_required_keys():
    # The 2-layer verdict gating reads "all", "no_POINT", "mechanism_only"
    required = {"all", "no_POINT", "mechanism_only"}
    assert required <= set(ABLATION_FEATURE_SETS.keys())


def test_ablation_no_POINT_excludes_point_block():
    # Amended verdict's no_POINT gate is load-bearing because POINT mutations
    # are essentially labels-in-genome-form (gyrA_S83L IS the cipro label).
    # MUST exclude "point" from this set or the test is circular.
    blocks = ABLATION_FEATURE_SETS["no_POINT"]
    assert "point" not in blocks
    assert "acquired" in blocks
    assert "kmer" in blocks
    assert "mlst" in blocks


def test_ablation_mechanism_only_includes_point_and_acquired():
    blocks = ABLATION_FEATURE_SETS["mechanism_only"]
    assert set(blocks) == {"point", "acquired"}


def test_ablation_all_includes_all_four_blocks():
    blocks = ABLATION_FEATURE_SETS["all"]
    assert set(blocks) == {"point", "acquired", "kmer", "mlst"}


def test_ablation_single_block_sets_have_one_entry_each():
    assert ABLATION_FEATURE_SETS["POINT_only"] == ("point",)
    assert ABLATION_FEATURE_SETS["kmer_only"] == ("kmer",)
    assert ABLATION_FEATURE_SETS["MLST_only"] == ("mlst",)


# ---- _load_amrfinder_features (JSON parsing) --------------------------------


def _write_audit_json(path: Path, per_strain: list[dict]) -> None:
    path.write_text(json.dumps({"per_strain": per_strain}), encoding="utf-8")


def test_load_amrfinder_features_extracts_quinolone_points(tmp_path: Path):
    json_path = tmp_path / "mech.json"
    _write_audit_json(json_path, [
        {"strain_id": "S1", "status": "OK", "hits": [
            {"kind": "mutation", "class": "QUINOLONE", "symbol": "gyrA_S83L"},
        ]},
    ])
    points, acquired = _load_amrfinder_features(json_path)
    assert "gyrA_S83L" in points["S1"]
    assert acquired["S1"] == set()


def test_load_amrfinder_features_extracts_acquired_genes(tmp_path: Path):
    json_path = tmp_path / "mech.json"
    _write_audit_json(json_path, [
        {"strain_id": "S1", "status": "OK", "hits": [
            {"kind": "acquired", "class": "BETA-LACTAM", "symbol": "blaCMY-2"},
        ]},
    ])
    points, acquired = _load_amrfinder_features(json_path)
    assert "blaCMY-2" in acquired["S1"]
    assert points["S1"] == set()


def test_load_amrfinder_features_skips_non_ok_status(tmp_path: Path):
    # Strains with status MISSING_FASTA / DOCKER_FAIL are dropped
    json_path = tmp_path / "mech.json"
    _write_audit_json(json_path, [
        {"strain_id": "S1", "status": "MISSING_FASTA", "hits": []},
        {"strain_id": "S2", "status": "OK", "hits": [
            {"kind": "mutation", "class": "QUINOLONE", "symbol": "gyrA_S83L"},
        ]},
    ])
    points, _ = _load_amrfinder_features(json_path)
    assert "S1" not in points
    assert "S2" in points


def test_load_amrfinder_features_filters_non_quinolone_mutations(tmp_path: Path):
    # mutation rows with non-QUINOLONE class are NOT included in points
    json_path = tmp_path / "mech.json"
    _write_audit_json(json_path, [
        {"strain_id": "S1", "status": "OK", "hits": [
            {"kind": "mutation", "class": "MULTIDRUG", "symbol": "marR_V84WfsTer"},
            {"kind": "mutation", "class": "QUINOLONE", "symbol": "gyrA_S83L"},
        ]},
    ])
    points, _ = _load_amrfinder_features(json_path)
    # MULTIDRUG mutations are NOT loaded as POINT features here (different
    # filter than the mechanism audit; curated baseline focuses on cipro-direct)
    assert "gyrA_S83L" in points["S1"]
    assert "marR_V84WfsTer" not in points["S1"]


def test_load_amrfinder_features_fluoroquinolone_class_also_kept(tmp_path: Path):
    json_path = tmp_path / "mech.json"
    _write_audit_json(json_path, [
        {"strain_id": "S1", "status": "OK", "hits": [
            {"kind": "mutation", "class": "FLUOROQUINOLONE", "symbol": "gyrA_D87N"},
        ]},
    ])
    points, _ = _load_amrfinder_features(json_path)
    assert "gyrA_D87N" in points["S1"]
