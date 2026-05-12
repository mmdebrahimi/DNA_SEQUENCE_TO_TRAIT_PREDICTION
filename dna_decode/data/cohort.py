"""Step 6 — Strain / AST cohort catalog (drug-first).

Builds a Phase 1 cohort by selecting strains drug-first: for each target drug,
ensures >= target_per_drug strains with AST labels for that drug; maximizes
the 3-drug intersection (strains labeled across cipro + ceftriaxone + tet);
maximizes MLST diversity within each per-drug pool.

Assembly-quality filter replaces the rejected "complete-circle only" gate
(post-tech-plan brainstorm M1). AMRFinder plasmid/chromosome metadata is
informational — helps interpret attribution downstream, doesn't gate
selection.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


COHORT_COLUMNS = (
    "strain_id",
    "mlst",
    "country",
    "year",
    "contig_count",
    "n50",
    "plasmid_resistance_genes",
    "chromosome_resistance_genes",
)


class CohortConstructionError(Exception):
    """Cohort cannot be built that satisfies the criteria (insufficient strains)."""


@dataclass(frozen=True)
class CandidateStrain:
    """A strain candidate considered for cohort inclusion.

    `ast_labels`: drug name (lowercase) -> binary label (0=susceptible, 1=resistant).
    Strains without an AST label for a given drug are simply not in that dict —
    no NaN sentinel needed.
    """

    strain_id: str
    mlst: str = ""
    country: str = ""
    year: int = 0
    contig_count: int = 0
    n50: int = 0
    ast_labels: dict[str, int] = field(default_factory=dict)
    plasmid_resistance_genes: tuple[str, ...] = ()
    chromosome_resistance_genes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CohortSelectionCriteria:
    """Filter + selection knobs for drug-first cohort construction."""

    target_per_drug: int = 150
    three_drug_intersection_target: int = 75
    assembly_contig_count_max: int = 500
    assembly_n50_min: int = 50_000
    # Plasmid annotation is informational only in v1; does not gate selection
    plasmid_localization_required: bool = False


@dataclass
class StrainCohort:
    """Output of build_cohort.

    `per_drug_strain_ids`: drug name (lowercase) -> list of strain_ids included
    for that drug's training pool.

    `three_drug_intersection`: sorted list of strain_ids labeled for all 3
    Phase 1 drugs (cipro + ceftriaxone + tet). Used by Step 9's joint training
    fallback when single-drug pools are sparse.
    """

    strains: list[CandidateStrain]
    per_drug_strain_ids: dict[str, list[str]]
    three_drug_intersection: list[str]

    def __len__(self) -> int:
        return len(self.strains)

    def strain_by_id(self, strain_id: str) -> CandidateStrain | None:
        for s in self.strains:
            if s.strain_id == strain_id:
                return s
        return None


def _filter_by_assembly_quality(
    candidates: list[CandidateStrain],
    criteria: CohortSelectionCriteria,
) -> list[CandidateStrain]:
    """Apply contig-count + N50 thresholds."""
    return [
        c
        for c in candidates
        if c.contig_count <= criteria.assembly_contig_count_max
        and c.n50 >= criteria.assembly_n50_min
    ]


def _mlst_balanced_selection(
    pool: list[CandidateStrain],
    target_n: int,
) -> list[CandidateStrain]:
    """Select target_n strains while maximizing MLST diversity.

    Round-robin pick: cycle through MLSTs and take one strain at a time from
    each MLST bucket until target_n strains are picked OR pool exhausted.

    Strains with empty `mlst` field go into a single "unknown" bucket — they
    still get picked but don't get over-prioritized.
    """
    if target_n <= 0 or not pool:
        return []

    # Group by MLST
    by_mlst: dict[str, list[CandidateStrain]] = {}
    for strain in pool:
        key = strain.mlst or "unknown"
        by_mlst.setdefault(key, []).append(strain)

    # Round-robin pick until target_n satisfied or all buckets empty
    selected: list[CandidateStrain] = []
    while len(selected) < target_n and any(by_mlst.values()):
        for mlst_key in sorted(by_mlst.keys()):  # deterministic order
            bucket = by_mlst[mlst_key]
            if not bucket:
                continue
            selected.append(bucket.pop(0))
            if len(selected) >= target_n:
                break
    return selected


def build_cohort(
    candidates: list[CandidateStrain],
    drugs: tuple[str, ...],
    criteria: CohortSelectionCriteria | None = None,
) -> StrainCohort:
    """Drug-first cohort construction.

    Algorithm:
      1. Apply assembly-quality threshold (contig_count + N50).
      2. For each target drug, select up to target_per_drug strains with an
         AST label for that drug, balancing MLST diversity.
      3. Raise CohortConstructionError per drug if target unmet.
      4. Compute 3-drug intersection (strains in ALL drug pools). Raise if
         below intersection target.
      5. Materialize the union of all per-drug pools as the cohort's strain
         list.

    Returns:
        StrainCohort with per_drug pools + intersection + metadata.
    """
    criteria = criteria or CohortSelectionCriteria()

    qualified = _filter_by_assembly_quality(candidates, criteria)
    if not qualified:
        raise CohortConstructionError(
            f"No candidate strains passed assembly-quality threshold "
            f"(contig_count <= {criteria.assembly_contig_count_max}, "
            f"n50 >= {criteria.assembly_n50_min})"
        )

    per_drug_strain_ids: dict[str, list[str]] = {}
    for drug in drugs:
        drug_lower = drug.lower()
        labeled = [c for c in qualified if drug_lower in c.ast_labels]
        selected = _mlst_balanced_selection(labeled, criteria.target_per_drug)
        if len(selected) < criteria.target_per_drug:
            raise CohortConstructionError(
                f"Drug {drug!r}: only {len(selected)} qualifying strains "
                f"with AST labels (target {criteria.target_per_drug})"
            )
        per_drug_strain_ids[drug_lower] = [s.strain_id for s in selected]

    # 3-drug intersection: strains present in EVERY per-drug pool
    if len(drugs) >= 2:
        intersection: set[str] = set(per_drug_strain_ids[drugs[0].lower()])
        for drug in drugs[1:]:
            intersection &= set(per_drug_strain_ids[drug.lower()])
    else:
        intersection = set(per_drug_strain_ids[drugs[0].lower()])

    if len(intersection) < criteria.three_drug_intersection_target:
        raise CohortConstructionError(
            f"3-drug intersection: {len(intersection)} strains "
            f"(target {criteria.three_drug_intersection_target})"
        )

    all_selected_ids: set[str] = set()
    for ids in per_drug_strain_ids.values():
        all_selected_ids.update(ids)
    selected_strains = [c for c in qualified if c.strain_id in all_selected_ids]

    return StrainCohort(
        strains=selected_strains,
        per_drug_strain_ids=per_drug_strain_ids,
        three_drug_intersection=sorted(intersection),
    )


def save_cohort(cohort: StrainCohort, parquet_path: Path | str) -> Path:
    """Persist cohort metadata as a parquet file.

    Schema follows COHORT_COLUMNS plus per-drug binary label columns (e.g.,
    `ast_ciprofloxacin`, `ast_ceftriaxone`, `ast_tetracycline`). Missing labels
    are NaN.
    """
    out_path = Path(parquet_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    all_drugs = sorted({d for s in cohort.strains for d in s.ast_labels.keys()})
    for s in cohort.strains:
        row: dict[str, object] = {
            "strain_id": s.strain_id,
            "mlst": s.mlst,
            "country": s.country,
            "year": s.year,
            "contig_count": s.contig_count,
            "n50": s.n50,
            "plasmid_resistance_genes": ";".join(s.plasmid_resistance_genes),
            "chromosome_resistance_genes": ";".join(s.chromosome_resistance_genes),
        }
        for drug in all_drugs:
            row[f"ast_{drug}"] = s.ast_labels.get(drug, None)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_parquet(out_path, index=False)
    return out_path


def download_cohort_genomes(
    cohort: StrainCohort,
    cache_root: Path | str,
    max_workers: int = 4,
) -> dict[str, Path]:
    """Materialize genomes for every strain in the cohort.

    Returns mapping strain_id -> per-accession cache directory (the per-strain
    output of `refseq.download_genome`).

    Uses `concurrent.futures.ThreadPoolExecutor` with `max_workers=4` by default
    (NCBI Datasets API politeness limit). Strains that fail to download surface
    via the exception type from `refseq` — the caller decides whether to
    skip + continue or abort the cohort run.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from dna_decode.data import refseq

    out: dict[str, Path] = {}
    if not cohort.strains:
        return out

    def _fetch(strain_id: str) -> tuple[str, Path]:
        path = refseq.download_genome(strain_id, cache_root)
        return strain_id, path

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(_fetch, s.strain_id): s.strain_id for s in cohort.strains}
        for fut in as_completed(futures):
            sid, path = fut.result()
            out[sid] = path

    return out
