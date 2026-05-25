"""Tests for the --build mode of scripts/bvbrc_strict_mic_4drug_census.py.

EP-2.5 (per plans/Post_V0_EP_Ladder_Plan.md): the census script now does
feasibility report + optional cohort build in ONE invocation. Tests the
build_cohort_from_census helper against synthetic census results + metadata.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest


def _synthetic_metadata(n_r: int = 5, n_s: int = 5) -> dict:
    """Build a metadata dict with N_R R-pool strains + N_S S-pool strains.
    Each strain has a unique MLST + accession + adequate assembly QC.
    """
    meta = {}
    for i in range(n_r):
        sid = f"strain_R_{i}"
        meta[sid] = {
            "assembly_accession": f"GCA_{i:09d}.1",
            "mlst": f"ST{100 + i}",
            "country": "USA",
            "year": "2024",
            "contig_count": 100,
            "n50": 150_000,
        }
    for i in range(n_s):
        sid = f"strain_S_{i}"
        meta[sid] = {
            "assembly_accession": f"GCA_{1000 + i:09d}.1",
            "mlst": f"ST{200 + i}",
            "country": "USA",
            "year": "2024",
            "contig_count": 100,
            "n50": 150_000,
        }
    return meta


def _synthetic_census_result(drug: str = "ciprofloxacin", n_strict_r: int = 5, n_strict_s: int = 5,
                              n_decisive_r: int = 2, n_decisive_s: int = 2) -> dict:
    """Build a census_drug-style result with explicit feasible_strain_ids."""
    return {
        "drug": drug,
        "feasible_strain_ids": {
            "strict_r": [f"strain_R_{i}" for i in range(n_strict_r)],
            "strict_s": [f"strain_S_{i}" for i in range(n_strict_s)],
            "decisive_r": [f"strain_R_decisive_{i}" for i in range(n_decisive_r)],
            "decisive_s": [f"strain_S_decisive_{i}" for i in range(n_decisive_s)],
        },
    }


def test_build_strict_mode_emits_cohort_parquet(tmp_path: Path):
    """Strict mode picks HIGH_R + HIGH_S only; ignores DECISIVE."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    metadata = _synthetic_metadata(n_r=5, n_s=5)
    census = _synthetic_census_result(n_strict_r=5, n_strict_s=5, n_decisive_r=2, n_decisive_s=2)
    out = tmp_path / "cipro_strict_cohort.parquet"
    result = build_cohort_from_census(
        census, metadata,
        label_quality="strict",
        target_per_class=5,
        balance_slack=1,
        output_path=out,
    )
    assert result["r_count"] == 5
    assert result["s_count"] == 5
    assert result["pool_sizes"]["r_pool"] == 5  # strict only; decisive ignored
    assert result["pool_sizes"]["s_pool"] == 5
    assert out.exists()


def test_build_relaxed_mode_includes_decisive(tmp_path: Path):
    """Relaxed mode pulls from HIGH + DECISIVE in both directions."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    # Add decisive metadata
    metadata = _synthetic_metadata(n_r=5, n_s=5)
    metadata["strain_R_decisive_0"] = {
        "assembly_accession": "GCA_999999900.1", "mlst": "ST_DR_0",
        "country": "USA", "year": "2024", "contig_count": 100, "n50": 150_000,
    }
    metadata["strain_R_decisive_1"] = {
        "assembly_accession": "GCA_999999901.1", "mlst": "ST_DR_1",
        "country": "USA", "year": "2024", "contig_count": 100, "n50": 150_000,
    }
    metadata["strain_S_decisive_0"] = {
        "assembly_accession": "GCA_999999902.1", "mlst": "ST_DS_0",
        "country": "USA", "year": "2024", "contig_count": 100, "n50": 150_000,
    }
    metadata["strain_S_decisive_1"] = {
        "assembly_accession": "GCA_999999903.1", "mlst": "ST_DS_1",
        "country": "USA", "year": "2024", "contig_count": 100, "n50": 150_000,
    }
    census = _synthetic_census_result(n_strict_r=5, n_strict_s=5, n_decisive_r=2, n_decisive_s=2)
    out = tmp_path / "cipro_relaxed_cohort.parquet"
    result = build_cohort_from_census(
        census, metadata,
        label_quality="relaxed",
        target_per_class=7,
        balance_slack=0,
        output_path=out,
    )
    assert result["pool_sizes"]["r_pool"] == 7  # 5 strict + 2 decisive
    assert result["pool_sizes"]["s_pool"] == 7
    assert result["r_count"] == 7
    assert result["s_count"] == 7
    assert out.exists()


def test_build_fails_when_pool_too_small(tmp_path: Path):
    """Target > pool - slack -> ValueError, not silent under-target cohort."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    metadata = _synthetic_metadata(n_r=3, n_s=5)
    census = _synthetic_census_result(n_strict_r=3, n_strict_s=5, n_decisive_r=0, n_decisive_s=0)
    out = tmp_path / "cipro_strict_cohort.parquet"
    with pytest.raises(ValueError, match="R pool size"):
        build_cohort_from_census(
            census, metadata,
            label_quality="strict",
            target_per_class=5,
            balance_slack=0,
            output_path=out,
        )


def test_build_drops_strains_missing_mlst(tmp_path: Path):
    """Strain in feasible list but no MLST in metadata -> dropped from candidate pool."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    metadata = _synthetic_metadata(n_r=5, n_s=5)
    # Strip MLST from one strain
    metadata["strain_R_0"]["mlst"] = ""
    census = _synthetic_census_result(n_strict_r=5, n_strict_s=5, n_decisive_r=0, n_decisive_s=0)
    out = tmp_path / "cipro_strict_cohort.parquet"
    result = build_cohort_from_census(
        census, metadata,
        label_quality="strict",
        target_per_class=4,
        balance_slack=0,
        output_path=out,
    )
    assert result["pool_sizes"]["r_pool"] == 4  # 5 strict_r - 1 missing MLST
    assert result["r_count"] == 4


def test_build_drops_strains_missing_accession(tmp_path: Path):
    """Strain with empty assembly_accession -> dropped."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    metadata = _synthetic_metadata(n_r=5, n_s=5)
    metadata["strain_R_0"]["assembly_accession"] = ""
    census = _synthetic_census_result(n_strict_r=5, n_strict_s=5, n_decisive_r=0, n_decisive_s=0)
    out = tmp_path / "cipro_strict_cohort.parquet"
    result = build_cohort_from_census(
        census, metadata,
        label_quality="strict",
        target_per_class=4,
        balance_slack=0,
        output_path=out,
    )
    assert result["pool_sizes"]["r_pool"] == 4
    assert result["r_count"] == 4


def test_build_rejects_unknown_label_quality(tmp_path: Path):
    """label_quality outside {strict, relaxed} -> ValueError."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    metadata = _synthetic_metadata()
    census = _synthetic_census_result()
    with pytest.raises(ValueError, match="unknown label_quality"):
        build_cohort_from_census(
            census, metadata,
            label_quality="permissive",  # not a real tier
            target_per_class=5,
            balance_slack=0,
            output_path=tmp_path / "x.parquet",
        )


def test_build_cohort_output_dir_created(tmp_path: Path):
    """Parent dirs of cohort_output are created if missing."""
    from scripts.bvbrc_strict_mic_4drug_census import build_cohort_from_census
    metadata = _synthetic_metadata(n_r=5, n_s=5)
    census = _synthetic_census_result(n_strict_r=5, n_strict_s=5, n_decisive_r=0, n_decisive_s=0)
    nested = tmp_path / "deep" / "nest" / "cipro.parquet"
    assert not nested.parent.exists()
    build_cohort_from_census(
        census, metadata,
        label_quality="strict",
        target_per_class=5,
        balance_slack=0,
        output_path=nested,
    )
    assert nested.exists()
