"""Tests for scripts/build_mini_cohort.py — Gate B mini-cohort selector."""
from __future__ import annotations

from pathlib import Path

import pytest


def _make_strain(strain_id: str, contigs: int, n50: int, mlst: str, label: int):
    """Construct a CandidateStrain with cipro AST label."""
    from dna_decode.data.cohort import CandidateStrain

    return CandidateStrain(
        strain_id=strain_id,
        assembly_accession=f"GCF_{strain_id}.1",
        mlst=mlst,
        contig_count=contigs,
        n50=n50,
        ast_labels={"ciprofloxacin": label},
    )


def _write_cohort(tmp_path: Path, strains, drug="ciprofloxacin") -> Path:
    """Build a tiny StrainCohort parquet for selection tests."""
    from dna_decode.data.cohort import StrainCohort, save_cohort

    drug_ids = [s.strain_id for s in strains if drug in s.ast_labels]
    cohort = StrainCohort(
        strains=strains,
        per_drug_strain_ids={drug: drug_ids},
        three_drug_intersection=drug_ids,
    )
    p = tmp_path / "source_cohort.parquet"
    save_cohort(cohort, p)
    return p


# ---- happy path ----


def test_select_mini_cohort_balanced(tmp_path: Path):
    """Pick 3 R + 3 S highest-quality strains for cipro."""
    from dna_decode.data.cohort import load_cohort
    from scripts.build_mini_cohort import select_mini_cohort

    strains = [
        _make_strain("R1", 2, 5_000_000, "ST131", 1),
        _make_strain("R2", 10, 200_000, "ST10", 1),
        _make_strain("R3", 4, 4_000_000, "ST167", 1),
        _make_strain("S1", 3, 4_500_000, "ST4", 0),
        _make_strain("S2", 8, 100_000, "ST155", 0),
        _make_strain("S3", 5, 3_800_000, "ST69", 0),
    ]
    src = _write_cohort(tmp_path, strains)
    mini = select_mini_cohort(load_cohort(src), "ciprofloxacin", per_class=3)
    assert len(mini) == 6
    drug_ids = mini.per_drug_strain_ids["ciprofloxacin"]
    assert len(drug_ids) == 6
    r_count = sum(1 for s in mini.strains if s.ast_labels.get("ciprofloxacin") == 1)
    s_count = sum(1 for s in mini.strains if s.ast_labels.get("ciprofloxacin") == 0)
    assert r_count == 3
    assert s_count == 3


def test_quality_sort_picks_low_contig_high_n50_first(tmp_path: Path):
    """Quality sort: contig_count asc, n50 desc."""
    from dna_decode.data.cohort import load_cohort
    from scripts.build_mini_cohort import select_mini_cohort

    strains = [
        _make_strain("R1", 50, 100_000, "STa", 1),    # low quality
        _make_strain("R2", 2, 5_000_000, "STb", 1),   # high quality (lowest contig, highest n50)
        _make_strain("R3", 5, 2_000_000, "STc", 1),   # mid
        _make_strain("S1", 50, 80_000, "STd", 0),
        _make_strain("S2", 2, 4_900_000, "STe", 0),   # high quality
        _make_strain("S3", 8, 1_500_000, "STf", 0),
    ]
    src = _write_cohort(tmp_path, strains)
    mini = select_mini_cohort(load_cohort(src), "ciprofloxacin", per_class=1)
    assert {s.strain_id for s in mini.strains} == {"R2", "S2"}


def test_per_class_2(tmp_path: Path):
    """Per-class > 1: still balanced, still quality-sorted."""
    from dna_decode.data.cohort import load_cohort
    from scripts.build_mini_cohort import select_mini_cohort

    strains = [
        _make_strain("R1", 2, 5_000_000, "STa", 1),
        _make_strain("R2", 3, 4_000_000, "STb", 1),
        _make_strain("R3", 5, 1_000_000, "STc", 1),  # excluded (3rd quality rank)
        _make_strain("S1", 4, 4_500_000, "STd", 0),
        _make_strain("S2", 6, 2_000_000, "STe", 0),
        _make_strain("S3", 8, 1_000_000, "STf", 0),  # excluded
    ]
    src = _write_cohort(tmp_path, strains)
    mini = select_mini_cohort(load_cohort(src), "ciprofloxacin", per_class=2)
    assert {s.strain_id for s in mini.strains} == {"R1", "R2", "S1", "S2"}


# ---- error paths ----


def test_unknown_drug_raises(tmp_path: Path):
    from dna_decode.data.cohort import load_cohort
    from scripts.build_mini_cohort import select_mini_cohort

    strains = [_make_strain("R1", 2, 5_000_000, "STa", 1)]
    src = _write_cohort(tmp_path, strains)
    with pytest.raises(ValueError, match="not in source"):
        select_mini_cohort(load_cohort(src), "novobiocin", per_class=1)


def test_insufficient_resistant_strains_raises(tmp_path: Path):
    from dna_decode.data.cohort import load_cohort
    from scripts.build_mini_cohort import select_mini_cohort

    strains = [
        _make_strain("R1", 2, 5_000_000, "STa", 1),
        _make_strain("S1", 3, 4_500_000, "STb", 0),
        _make_strain("S2", 4, 4_000_000, "STc", 0),
    ]
    src = _write_cohort(tmp_path, strains)
    with pytest.raises(ValueError, match="resistant"):
        select_mini_cohort(load_cohort(src), "ciprofloxacin", per_class=2)


def test_insufficient_susceptible_strains_raises(tmp_path: Path):
    from dna_decode.data.cohort import load_cohort
    from scripts.build_mini_cohort import select_mini_cohort

    strains = [
        _make_strain("R1", 2, 5_000_000, "STa", 1),
        _make_strain("R2", 3, 4_500_000, "STb", 1),
        _make_strain("S1", 4, 4_000_000, "STc", 0),
    ]
    src = _write_cohort(tmp_path, strains)
    with pytest.raises(ValueError, match="susceptible"):
        select_mini_cohort(load_cohort(src), "ciprofloxacin", per_class=2)


# ---- CLI surface ----


def test_main_missing_source_exits_2(tmp_path: Path):
    from scripts.build_mini_cohort import main

    exit_code = main(
        [
            "--source", str(tmp_path / "missing.parquet"),
            "--output", str(tmp_path / "out.parquet"),
            "--per-class", "1",
        ]
    )
    assert exit_code == 2


def test_main_happy_path_writes_output(tmp_path: Path):
    from scripts.build_mini_cohort import main

    strains = [
        _make_strain("R1", 2, 5_000_000, "STa", 1),
        _make_strain("R2", 3, 4_500_000, "STb", 1),
        _make_strain("S1", 4, 4_000_000, "STc", 0),
        _make_strain("S2", 5, 3_500_000, "STd", 0),
    ]
    src = _write_cohort(tmp_path, strains)
    out = tmp_path / "mini.parquet"
    exit_code = main(
        [
            "--source", str(src),
            "--output", str(out),
            "--per-class", "2",
        ]
    )
    assert exit_code == 0
    assert out.exists()
