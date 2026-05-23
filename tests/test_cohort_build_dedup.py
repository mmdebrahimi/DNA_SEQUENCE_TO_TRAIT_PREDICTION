"""Regression tests for assembly_accession uniqueness assertion in build_cohort.

Pins the 2026-05-22 LESSON on duplicate-accession LOSO leakage by construction:
when two strain_ids share an assembly_accession, LOSO puts each strain in the
other's training set -- same-genome train/test leakage. The fix asserts at the
cohort-builder layer; this test file pins the assertion.

Override via build_cohort(..., allow_duplicate_accessions=True) for the rare
case where intentional same-genome registration is wanted.
"""
from __future__ import annotations

import pytest

from dna_decode.data.cohort import (
    CandidateStrain,
    CohortConstructionError,
    CohortSelectionCriteria,
    build_cohort,
)


def _mk_candidate(
    strain_id: str,
    accession: str = "",
    mlst: str = "ST10",
    contig_count: int = 50,
    n50: int = 100_000,
    drugs: tuple[str, ...] = ("ciprofloxacin",),
    labels: tuple[int, ...] | None = None,
) -> CandidateStrain:
    if labels is None:
        labels = (1,) * len(drugs)
    ast = dict(zip(drugs, labels))
    return CandidateStrain(
        strain_id=strain_id,
        assembly_accession=accession,
        mlst=mlst,
        contig_count=contig_count,
        n50=n50,
        ast_labels=ast,
    )


def test_build_cohort_raises_on_duplicate_accession():
    """The 2026-05-22 leakage-by-construction class -- pinned regression test."""
    candidates = [
        _mk_candidate("562.109860", accession="GCA_025200635.1", mlst="ST131"),
        _mk_candidate("562.111036", accession="GCA_025200635.1", mlst="ST131"),
    ]
    with pytest.raises(CohortConstructionError, match="Duplicate assembly_accession"):
        build_cohort(
            candidates,
            drugs=("ciprofloxacin",),
            criteria=CohortSelectionCriteria(target_per_drug=1, three_drug_intersection_target=0),
        )


def test_build_cohort_accepts_unique_accessions():
    """Unique accessions across all strain_ids -- builder proceeds normally."""
    candidates = [
        _mk_candidate(f"strain_{i}", accession=f"GCA_{i:09d}.1", mlst=f"ST{i}")
        for i in range(5)
    ]
    cohort = build_cohort(
        candidates,
        drugs=("ciprofloxacin",),
        criteria=CohortSelectionCriteria(target_per_drug=5, three_drug_intersection_target=0),
    )
    assert len(cohort.strains) == 5


def test_build_cohort_allows_override_when_explicit():
    """User can override the assertion via allow_duplicate_accessions=True."""
    candidates = [
        _mk_candidate("a", accession="GCA_025200635.1", mlst="ST131"),
        _mk_candidate("b", accession="GCA_025200635.1", mlst="ST131"),
    ]
    cohort = build_cohort(
        candidates,
        drugs=("ciprofloxacin",),
        criteria=CohortSelectionCriteria(target_per_drug=2, three_drug_intersection_target=0),
        allow_duplicate_accessions=True,
    )
    assert len(cohort.strains) == 2


def test_build_cohort_ignores_empty_accession_in_dedup_check():
    """Empty assembly_accession is excluded from the dedup check (legit case)."""
    candidates = [
        _mk_candidate("a", accession="", mlst="ST10"),
        _mk_candidate("b", accession="", mlst="ST20"),
        _mk_candidate("c", accession="GCA_000005845.2", mlst="ST30"),
    ]
    cohort = build_cohort(
        candidates,
        drugs=("ciprofloxacin",),
        criteria=CohortSelectionCriteria(target_per_drug=3, three_drug_intersection_target=0),
    )
    assert len(cohort.strains) == 3


def test_build_cohort_raises_on_three_way_duplicate():
    """Three strain_ids sharing one accession -- still flagged."""
    candidates = [
        _mk_candidate(f"strain_{i}", accession="GCA_025200635.1", mlst=f"ST{i}")
        for i in range(3)
    ]
    with pytest.raises(CohortConstructionError, match="Duplicate assembly_accession"):
        build_cohort(
            candidates,
            drugs=("ciprofloxacin",),
            criteria=CohortSelectionCriteria(target_per_drug=1, three_drug_intersection_target=0),
        )


def test_build_cohort_error_message_includes_sample_strain_ids():
    """Error message surfaces the duplicate accession + at least one strain_id pair."""
    candidates = [
        _mk_candidate("562.109860", accession="GCA_025200635.1"),
        _mk_candidate("562.111036", accession="GCA_025200635.1"),
        _mk_candidate("562.222222", accession="GCA_111111111.1"),  # unique, control
    ]
    with pytest.raises(CohortConstructionError) as exc_info:
        build_cohort(
            candidates,
            drugs=("ciprofloxacin",),
            criteria=CohortSelectionCriteria(target_per_drug=1, three_drug_intersection_target=0),
        )
    msg = str(exc_info.value)
    assert "GCA_025200635.1" in msg
    assert "562.109860" in msg or "562.111036" in msg
