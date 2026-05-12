"""Tests for Step 6 — Strain / AST cohort catalog (drug-first)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from dna_decode.data.cohort import (
    CandidateStrain,
    CohortConstructionError,
    CohortSelectionCriteria,
    StrainCohort,
    _filter_by_assembly_quality,
    _mlst_balanced_selection,
    build_cohort,
    candidates_from_bvbrc_ast,
    load_cohort,
    save_cohort,
)


# ---- helpers ----


def _mk_candidate(
    strain_id: str,
    mlst: str = "ST10",
    contig_count: int = 50,
    n50: int = 100_000,
    drugs: tuple[str, ...] = ("ciprofloxacin", "ceftriaxone", "tetracycline"),
    labels: tuple[int, ...] | None = None,
) -> CandidateStrain:
    """Build a CandidateStrain with optional per-drug labels."""
    if labels is None:
        labels = (1,) * len(drugs)
    ast = dict(zip(drugs, labels))
    return CandidateStrain(
        strain_id=strain_id,
        mlst=mlst,
        contig_count=contig_count,
        n50=n50,
        ast_labels=ast,
    )


# ---- assembly-quality filter ----


def test_filter_excludes_high_contig_count():
    bad = _mk_candidate("s_bad", contig_count=1000)
    good = _mk_candidate("s_good", contig_count=50)
    out = _filter_by_assembly_quality(
        [bad, good], CohortSelectionCriteria(assembly_contig_count_max=500)
    )
    assert [s.strain_id for s in out] == ["s_good"]


def test_filter_excludes_low_n50():
    bad = _mk_candidate("s_low_n50", n50=10_000)
    good = _mk_candidate("s_good", n50=200_000)
    out = _filter_by_assembly_quality(
        [bad, good], CohortSelectionCriteria(assembly_n50_min=50_000)
    )
    assert [s.strain_id for s in out] == ["s_good"]


def test_filter_empty_returns_empty():
    assert _filter_by_assembly_quality([], CohortSelectionCriteria()) == []


# ---- MLST-balanced selection ----


def test_mlst_balanced_round_robin_across_buckets():
    """Picks should rotate across MLST buckets, not pack one bucket."""
    pool = [
        _mk_candidate("a1", mlst="ST131"),
        _mk_candidate("a2", mlst="ST131"),
        _mk_candidate("a3", mlst="ST131"),
        _mk_candidate("b1", mlst="ST10"),
        _mk_candidate("b2", mlst="ST10"),
        _mk_candidate("c1", mlst="ST73"),
    ]
    selected = _mlst_balanced_selection(pool, target_n=3)
    mlst_seen = {s.mlst for s in selected}
    # All 3 MLSTs represented in the first 3 picks
    assert mlst_seen == {"ST10", "ST131", "ST73"}


def test_mlst_balanced_target_exceeds_pool():
    pool = [_mk_candidate("a1", mlst="ST131"), _mk_candidate("a2", mlst="ST10")]
    selected = _mlst_balanced_selection(pool, target_n=5)
    assert len(selected) == 2  # exhausted pool


def test_mlst_balanced_unknown_mlst_handled():
    pool = [_mk_candidate("a1", mlst=""), _mk_candidate("b1", mlst="ST10")]
    selected = _mlst_balanced_selection(pool, target_n=2)
    assert {s.strain_id for s in selected} == {"a1", "b1"}


def test_mlst_balanced_empty_pool():
    assert _mlst_balanced_selection([], target_n=5) == []


def test_mlst_balanced_zero_target():
    pool = [_mk_candidate("a1", mlst="ST10")]
    assert _mlst_balanced_selection(pool, target_n=0) == []


# ---- build_cohort happy path ----


def test_build_cohort_basic_success():
    """All 3 drugs have 5+ qualifying strains with MLST diversity → GO."""
    candidates = []
    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    # 10 strains across 3 MLSTs, all 3 drugs labeled
    for i in range(10):
        mlst = ["ST131", "ST10", "ST73"][i % 3]
        candidates.append(_mk_candidate(f"s{i:03d}", mlst=mlst, drugs=drugs))

    criteria = CohortSelectionCriteria(
        target_per_drug=5, three_drug_intersection_target=3
    )
    cohort = build_cohort(candidates, drugs, criteria)

    for drug in drugs:
        assert len(cohort.per_drug_strain_ids[drug]) == 5
    assert len(cohort.three_drug_intersection) >= 3


def test_build_cohort_per_drug_diversity():
    """The per-drug pool should pick diverse MLSTs first."""
    candidates = []
    # 6 ST131 strains (favored by raw order) + 1 ST10 + 1 ST73, all drug-labeled
    for i in range(6):
        candidates.append(_mk_candidate(f"st131_{i}", mlst="ST131"))
    candidates.append(_mk_candidate("st10_0", mlst="ST10"))
    candidates.append(_mk_candidate("st73_0", mlst="ST73"))

    # target=3 → MLST-balanced pick covers all 3 MLSTs
    cohort = build_cohort(
        candidates,
        ("ciprofloxacin",),
        CohortSelectionCriteria(target_per_drug=3, three_drug_intersection_target=1),
    )
    picked_ids = set(cohort.per_drug_strain_ids["ciprofloxacin"])
    picked_mlsts = {c.mlst for c in candidates if c.strain_id in picked_ids}
    assert picked_mlsts == {"ST131", "ST10", "ST73"}


# ---- build_cohort failure paths ----


def test_build_cohort_no_assembly_passes_raises():
    """All candidates fail assembly filter → CohortConstructionError."""
    candidates = [_mk_candidate("s_bad", contig_count=10_000)]
    with pytest.raises(CohortConstructionError, match="assembly-quality"):
        build_cohort(
            candidates,
            ("ciprofloxacin",),
            CohortSelectionCriteria(assembly_contig_count_max=500, target_per_drug=1),
        )


def test_build_cohort_per_drug_short_raises():
    """One drug below target → CohortConstructionError naming the drug."""
    candidates = []
    for i in range(10):
        # All 10 have cipro labels; only 1 has ceftriaxone
        labels = (1, 0, 1) if i == 0 else (1, None, 1)
        ast = {"ciprofloxacin": 1, "tetracycline": 1}
        if i == 0:
            ast["ceftriaxone"] = 0
        candidates.append(CandidateStrain(
            strain_id=f"s{i:03d}", mlst="ST10", contig_count=50, n50=100_000, ast_labels=ast
        ))

    drugs = ("ciprofloxacin", "ceftriaxone", "tetracycline")
    with pytest.raises(CohortConstructionError, match="ceftriaxone"):
        build_cohort(
            candidates,
            drugs,
            CohortSelectionCriteria(target_per_drug=5, three_drug_intersection_target=1),
        )


def test_build_cohort_intersection_short_raises():
    """All per-drug targets pass but intersection fails → error."""
    # 10 strains, cipro-only. 10 different strains, ceftriaxone-only.
    candidates = []
    for i in range(10):
        candidates.append(CandidateStrain(
            strain_id=f"c{i:03d}", mlst="ST10", contig_count=50, n50=100_000,
            ast_labels={"ciprofloxacin": 1},
        ))
    for i in range(10):
        candidates.append(CandidateStrain(
            strain_id=f"t{i:03d}", mlst="ST131", contig_count=50, n50=100_000,
            ast_labels={"ceftriaxone": 1},
        ))

    with pytest.raises(CohortConstructionError, match="intersection"):
        build_cohort(
            candidates,
            ("ciprofloxacin", "ceftriaxone"),
            CohortSelectionCriteria(target_per_drug=5, three_drug_intersection_target=5),
        )


# ---- save_cohort ----


def test_save_cohort_writes_parquet_with_per_drug_columns(tmp_path: Path):
    candidates = [
        CandidateStrain(
            strain_id="s001",
            assembly_accession="GCF_X",
            mlst="ST10",
            country="USA",
            year=2024,
            contig_count=50,
            n50=200_000,
            ast_labels={"ciprofloxacin": 1, "ceftriaxone": 0},
            plasmid_resistance_genes=("blaCTX-M-15",),
        )
    ]
    cohort = StrainCohort(
        strains=candidates,
        per_drug_strain_ids={"ciprofloxacin": ["s001"]},
        three_drug_intersection=[],
    )
    out = save_cohort(cohort, tmp_path / "cohort.parquet")
    assert out.exists()
    df = pd.read_parquet(out)
    assert list(df["strain_id"]) == ["s001"]
    assert df["assembly_accession"].iloc[0] == "GCF_X"
    assert "ast_ciprofloxacin" in df.columns
    assert "ast_ceftriaxone" in df.columns
    assert "in_pool_ciprofloxacin" in df.columns
    assert "in_three_drug_intersection" in df.columns
    assert df["plasmid_resistance_genes"].iloc[0] == "blaCTX-M-15"


# ---- Wave 2.5 hardening: load_cohort round-trip ----


def test_load_cohort_round_trips_save(tmp_path: Path):
    """save_cohort → load_cohort returns equivalent StrainCohort."""
    candidates = [
        CandidateStrain(
            strain_id="s001",
            assembly_accession="GCF_X",
            mlst="ST10",
            contig_count=50,
            n50=100_000,
            ast_labels={"ciprofloxacin": 1, "ceftriaxone": 0},
            plasmid_resistance_genes=("blaCTX-M-15", "tetA"),
        ),
        CandidateStrain(
            strain_id="s002",
            assembly_accession="GCF_Y",
            mlst="ST131",
            contig_count=100,
            n50=80_000,
            ast_labels={"ciprofloxacin": 0},
        ),
    ]
    orig = StrainCohort(
        strains=candidates,
        per_drug_strain_ids={"ciprofloxacin": ["s001", "s002"], "ceftriaxone": ["s001"]},
        three_drug_intersection=["s001"],
    )
    save_cohort(orig, tmp_path / "cohort.parquet")
    loaded = load_cohort(tmp_path / "cohort.parquet")

    assert len(loaded.strains) == 2
    assert {s.strain_id for s in loaded.strains} == {"s001", "s002"}
    assert loaded.per_drug_strain_ids["ciprofloxacin"] == ["s001", "s002"]
    assert loaded.per_drug_strain_ids["ceftriaxone"] == ["s001"]
    assert loaded.three_drug_intersection == ["s001"]
    # Per-strain fields preserved
    s001 = loaded.strain_by_id("s001")
    assert s001.assembly_accession == "GCF_X"
    assert s001.ast_labels == {"ciprofloxacin": 1, "ceftriaxone": 0}
    assert "blaCTX-M-15" in s001.plasmid_resistance_genes


def test_load_cohort_missing_file_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_cohort(tmp_path / "missing.parquet")


# ---- Wave 2.5 hardening: candidates_from_bvbrc_ast adapter ----


def test_candidates_from_bvbrc_ast_groups_by_strain():
    """Multiple AST rows per strain → one CandidateStrain with per-drug labels."""
    ast = pd.DataFrame(
        {
            "strain_id": ["s001", "s001", "s002", "s002"],
            "antibiotic": ["ciprofloxacin", "ceftriaxone", "ciprofloxacin", "tetracycline"],
            "susceptibility_label": ["RESISTANT", "SUSCEPTIBLE", "SUSCEPTIBLE", "RESISTANT"],
        }
    )
    cands = candidates_from_bvbrc_ast(ast)
    by_id = {c.strain_id: c for c in cands}
    assert by_id["s001"].ast_labels == {"ciprofloxacin": 1, "ceftriaxone": 0}
    assert by_id["s002"].ast_labels == {"ciprofloxacin": 0, "tetracycline": 1}


def test_candidates_from_bvbrc_ast_uses_assembly_metadata():
    """assembly_metadata supplies accession + assembly QC + MLST per strain."""
    ast = pd.DataFrame(
        {
            "strain_id": ["s001"],
            "antibiotic": ["ciprofloxacin"],
            "susceptibility_label": ["R"],
        }
    )
    meta = {
        "s001": {
            "assembly_accession": "GCF_000005845.2",
            "contig_count": 50,
            "n50": 200_000,
            "mlst": "ST10",
            "country": "USA",
            "year": 2023,
            "plasmid_resistance_genes": ["blaCTX-M-15"],
        }
    }
    cands = candidates_from_bvbrc_ast(ast, assembly_metadata=meta)
    assert len(cands) == 1
    assert cands[0].assembly_accession == "GCF_000005845.2"
    assert cands[0].contig_count == 50
    assert cands[0].mlst == "ST10"
    assert cands[0].plasmid_resistance_genes == ("blaCTX-M-15",)


def test_candidates_from_bvbrc_ast_majority_vote_on_duplicates():
    """Multiple AST rows for same (strain, drug) → majority vote."""
    ast = pd.DataFrame(
        {
            "strain_id": ["s001", "s001", "s001"],
            "antibiotic": ["ciprofloxacin"] * 3,
            "susceptibility_label": ["R", "R", "S"],  # 2 R + 1 S → majority R
        }
    )
    cands = candidates_from_bvbrc_ast(ast)
    assert cands[0].ast_labels["ciprofloxacin"] == 1


def test_candidates_from_bvbrc_ast_no_metadata_falls_back_to_zero():
    """Missing assembly_metadata → contig_count=0 etc. (will fail assembly filter)."""
    ast = pd.DataFrame(
        {
            "strain_id": ["s001"],
            "antibiotic": ["ciprofloxacin"],
            "susceptibility_label": ["R"],
        }
    )
    cands = candidates_from_bvbrc_ast(ast)
    assert cands[0].assembly_accession == ""
    assert cands[0].contig_count == 0
    assert cands[0].n50 == 0


def test_candidates_from_bvbrc_ast_missing_column_raises():
    bad = pd.DataFrame({"foo": ["bar"]})
    with pytest.raises(ValueError, match="strain_id"):
        candidates_from_bvbrc_ast(bad)


# ---- Wave 2.5 hardening: download_cohort_genomes accession routing ----


def test_download_cohort_genomes_requires_assembly_accession():
    """Strains without assembly_accession raise CohortConstructionError."""
    from unittest.mock import patch

    strains = [CandidateStrain(strain_id="s001", assembly_accession="")]
    cohort = StrainCohort(
        strains=strains,
        per_drug_strain_ids={"ciprofloxacin": ["s001"]},
        three_drug_intersection=[],
    )

    from dna_decode.data.cohort import download_cohort_genomes

    with pytest.raises(CohortConstructionError, match="assembly_accession"):
        download_cohort_genomes(cohort, "/tmp/cache_dir")


def test_download_cohort_genomes_uses_assembly_accession_not_strain_id(tmp_path: Path):
    """refseq.download_genome is called with assembly_accession, NOT strain_id."""
    from unittest.mock import MagicMock, patch

    from dna_decode.data import cohort as cohort_mod

    strains = [
        CandidateStrain(strain_id="bv_511145.12", assembly_accession="GCF_000005845.2")
    ]
    cohort = StrainCohort(
        strains=strains,
        per_drug_strain_ids={"ciprofloxacin": ["bv_511145.12"]},
        three_drug_intersection=[],
    )

    fake_path = tmp_path / "GCF_000005845.2"
    fake_path.mkdir()
    mock_download = MagicMock(return_value=fake_path)
    with patch("dna_decode.data.refseq.download_genome", mock_download):
        out = cohort_mod.download_cohort_genomes(cohort, tmp_path)

    # Called with the accession, not the BV-BRC strain_id
    call_args = mock_download.call_args
    assert call_args[0][0] == "GCF_000005845.2"
    # But the output map keys by the BV-BRC strain_id (consistent ID space downstream)
    assert "bv_511145.12" in out


# ---- StrainCohort accessors ----


def test_strain_cohort_strain_by_id_round_trip():
    s = _mk_candidate("s001")
    cohort = StrainCohort(strains=[s], per_drug_strain_ids={}, three_drug_intersection=[])
    assert cohort.strain_by_id("s001") is s
    assert cohort.strain_by_id("missing") is None


def test_strain_cohort_len_matches_strain_count():
    cohort = StrainCohort(
        strains=[_mk_candidate(f"s{i}") for i in range(7)],
        per_drug_strain_ids={},
        three_drug_intersection=[],
    )
    assert len(cohort) == 7
