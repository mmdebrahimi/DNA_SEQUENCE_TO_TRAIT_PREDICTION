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
    "assembly_accession",
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

    `strain_id`: BV-BRC genome id (e.g., "511145.12") — the AST-data primary key.
    `assembly_accession`: NCBI Datasets accession (e.g., "GCF_000005845.2") used
        by refseq.download_genome. Wave 2.5 hardening C4: BV-BRC and NCBI use
        different ID spaces; download paths require the accession, not the BV-BRC
        ID. Empty string means "not yet resolved" — download will raise.

    `ast_labels`: drug name (lowercase) -> binary label (0=susceptible, 1=resistant).
    Strains without an AST label for a given drug are simply not in that dict —
    no NaN sentinel needed.
    """

    strain_id: str
    assembly_accession: str = ""
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
    allow_duplicate_accessions: bool = False,
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
      6. Assert assembly_accession uniqueness (LOSO same-genome leakage
         guard; 2026-05-22 LESSON). Override via `allow_duplicate_accessions`.

    Returns:
        StrainCohort with per_drug pools + intersection + metadata.
    """
    criteria = criteria or CohortSelectionCriteria()

    # Pre-filter duplicate-accession check on input candidates -- catches the
    # bug class before filter steps drop one of the pair.
    if not allow_duplicate_accessions:
        accession_to_strain_ids: dict[str, list[str]] = {}
        for c in candidates:
            if c.assembly_accession:
                accession_to_strain_ids.setdefault(c.assembly_accession, []).append(c.strain_id)
        dups = {acc: sids for acc, sids in accession_to_strain_ids.items() if len(sids) > 1}
        if dups:
            sample = next(iter(dups.items()))
            raise CohortConstructionError(
                f"Duplicate assembly_accession in candidate set: {len(dups)} accession(s) "
                f"shared by multiple strain_ids (sample: {sample[0]} -> {sample[1]}). "
                f"Same-genome LOSO leakage would result. Dedup at source OR pass "
                f"allow_duplicate_accessions=True if intentional."
            )

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


def find_duplicate_accessions(
    cohort: StrainCohort,
    restrict_to_strain_ids: list[str] | None = None,
) -> dict[str, list[str]]:
    """Find assembly_accessions shared by ≥ 2 strain_ids in a cohort.

    Used by `cmd_train` to decide between strain-id CV (safe when no
    duplicates) vs accession CV (required when duplicates exist; matches the
    2026-05-22 LESSON on duplicate-accession LOSO leakage). The cohort
    builder asserts uniqueness at build time (see `build_cohort`'s
    accession-uniqueness assertion), but a pre-existing parquet built before
    that assertion landed may still contain duplicates.

    Args:
        cohort: StrainCohort to scan.
        restrict_to_strain_ids: optional filter. If provided, only strains
            whose strain_id is in this list contribute to the scan. Useful
            for per-drug duplicate detection (e.g., `cmd_train` passes the
            drug-specific strain list to ignore duplicates that don't
            participate in the drug pool).

    Returns:
        Mapping `accession -> [strain_ids]` for accessions shared by ≥ 2
        strain_ids. Empty assembly_accession values are excluded (they
        represent legit-missing accessions). Returns empty dict if no
        duplicates exist.
    """
    restrict_set = set(restrict_to_strain_ids) if restrict_to_strain_ids is not None else None
    accession_to_strain_ids: dict[str, list[str]] = {}
    for s in cohort.strains:
        if restrict_set is not None and s.strain_id not in restrict_set:
            continue
        acc = (s.assembly_accession or "").strip()
        if not acc:
            continue
        accession_to_strain_ids.setdefault(acc, []).append(s.strain_id)
    return {acc: sorted(sids) for acc, sids in accession_to_strain_ids.items() if len(sids) > 1}


def save_cohort(cohort: StrainCohort, parquet_path: Path | str) -> Path:
    """Persist cohort metadata as a parquet file.

    Schema follows COHORT_COLUMNS plus per-drug binary label columns (e.g.,
    `ast_ciprofloxacin`, `ast_ceftriaxone`, `ast_tetracycline`). Missing labels
    are NaN.

    Wave 2.5 hardening C5b: also persists per-drug pool membership as boolean
    columns (e.g., `in_pool_ciprofloxacin`) and `in_three_drug_intersection`
    so `load_cohort` can round-trip the StrainCohort without recomputing.
    """
    out_path = Path(parquet_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    all_drugs = sorted({d for s in cohort.strains for d in s.ast_labels.keys()})
    pool_drugs = sorted(cohort.per_drug_strain_ids.keys())
    intersection_set = set(cohort.three_drug_intersection)

    for s in cohort.strains:
        row: dict[str, object] = {
            "strain_id": s.strain_id,
            "assembly_accession": s.assembly_accession,
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
        for drug in pool_drugs:
            row[f"in_pool_{drug}"] = s.strain_id in set(cohort.per_drug_strain_ids[drug])
        row["in_three_drug_intersection"] = s.strain_id in intersection_set
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_parquet(out_path, index=False)
    return out_path


def load_cohort(parquet_path: Path | str) -> StrainCohort:
    """Round-trip counterpart to `save_cohort`.

    Reconstructs the full StrainCohort: strains list (with assembly_accession,
    per-strain ast_labels) + per_drug_strain_ids + three_drug_intersection.

    Wave 2.5 hardening C5b — closes save/load asymmetry that previously
    forced Wave 3 to recompute `build_cohort` on every load.
    """
    in_path = Path(parquet_path)
    if not in_path.exists():
        raise FileNotFoundError(f"Cohort parquet not found at {in_path}")

    df = pd.read_parquet(in_path)

    # Reconstruct per-strain candidates
    ast_cols = [c for c in df.columns if c.startswith("ast_")]
    pool_cols = [c for c in df.columns if c.startswith("in_pool_")]
    drugs_from_pools = [c[len("in_pool_"):] for c in pool_cols]

    strains: list[CandidateStrain] = []
    per_drug_strain_ids: dict[str, list[str]] = {d: [] for d in drugs_from_pools}
    three_drug_intersection: list[str] = []

    for _, row in df.iterrows():
        ast_labels: dict[str, int] = {}
        for c in ast_cols:
            val = row[c]
            if pd.notna(val):
                ast_labels[c[len("ast_"):]] = int(val)
        plasmid_genes = tuple(
            g for g in str(row.get("plasmid_resistance_genes", "") or "").split(";") if g
        )
        chromosome_genes = tuple(
            g for g in str(row.get("chromosome_resistance_genes", "") or "").split(";") if g
        )
        strains.append(
            CandidateStrain(
                strain_id=str(row["strain_id"]),
                assembly_accession=str(row.get("assembly_accession", "") or ""),
                mlst=str(row.get("mlst", "") or ""),
                country=str(row.get("country", "") or ""),
                year=int(row.get("year", 0) or 0),
                contig_count=int(row.get("contig_count", 0) or 0),
                n50=int(row.get("n50", 0) or 0),
                ast_labels=ast_labels,
                plasmid_resistance_genes=plasmid_genes,
                chromosome_resistance_genes=chromosome_genes,
            )
        )
        for drug in drugs_from_pools:
            if bool(row.get(f"in_pool_{drug}", False)):
                per_drug_strain_ids[drug].append(str(row["strain_id"]))
        if bool(row.get("in_three_drug_intersection", False)):
            three_drug_intersection.append(str(row["strain_id"]))

    return StrainCohort(
        strains=strains,
        per_drug_strain_ids=per_drug_strain_ids,
        three_drug_intersection=sorted(three_drug_intersection),
    )


def candidates_from_bvbrc_ast(
    ast_table: pd.DataFrame,
    assembly_metadata: dict[str, dict[str, object]] | None = None,
) -> list[CandidateStrain]:
    """Build CandidateStrain list from Wave 1's `load_bvbrc_ast` output.

    Wave 2.5 hardening C5a — closes the wired path from Step 5 → Step 6 that
    previously left CandidateStrain construction as exercise to the caller.

    Args:
        ast_table: DataFrame returned by `dna_decode.data.ast_data.load_bvbrc_ast`.
            Must contain columns `strain_id`, `antibiotic`, `susceptibility_label`.
        assembly_metadata: optional dict[strain_id, metadata-dict] holding
            per-strain NCBI assembly facts: `assembly_accession`, `contig_count`,
            `n50`, `mlst`, `country`, `year`, plus optional
            `plasmid_resistance_genes` / `chromosome_resistance_genes` tuples.
            When missing, defaults to empty/zero values — strains will fail
            the assembly-quality filter in build_cohort unless metadata is
            supplied (intentional: phase-1 cohort construction requires
            real assembly QC).

    Returns:
        list[CandidateStrain], one per unique strain_id. Per-drug labels are
        aggregated by majority vote when a strain has multiple AST rows for
        the same drug (rare but possible across testing standards).
    """
    if "strain_id" not in ast_table.columns:
        raise ValueError("ast_table is missing required column 'strain_id'")

    assembly_metadata = assembly_metadata or {}

    from dna_decode.data.ast_data import SUSCEPTIBILITY_BINARY_MAP

    candidates: list[CandidateStrain] = []
    for strain_id, group in ast_table.groupby("strain_id", sort=True):
        # Aggregate per-drug labels — majority vote when multiple rows per drug
        per_drug: dict[str, int] = {}
        for drug, drug_rows in group.groupby("antibiotic"):
            labels = [
                SUSCEPTIBILITY_BINARY_MAP[lbl]
                for lbl in drug_rows["susceptibility_label"]
                if lbl in SUSCEPTIBILITY_BINARY_MAP
            ]
            if not labels:
                continue
            # Majority vote; ties → resistant (conservative for clinical use)
            per_drug[str(drug)] = int(round(sum(labels) / len(labels)))

        meta = assembly_metadata.get(str(strain_id), {})
        candidates.append(
            CandidateStrain(
                strain_id=str(strain_id),
                assembly_accession=str(meta.get("assembly_accession", "") or ""),
                mlst=str(meta.get("mlst", "") or ""),
                country=str(meta.get("country", "") or ""),
                year=int(meta.get("year", 0) or 0),
                contig_count=int(meta.get("contig_count", 0) or 0),
                n50=int(meta.get("n50", 0) or 0),
                ast_labels=per_drug,
                plasmid_resistance_genes=tuple(meta.get("plasmid_resistance_genes", ()) or ()),
                chromosome_resistance_genes=tuple(
                    meta.get("chromosome_resistance_genes", ()) or ()
                ),
            )
        )

    return candidates


def download_cohort_genomes(
    cohort: StrainCohort,
    cache_root: Path | str,
    max_workers: int = 4,
) -> dict[str, Path]:
    """Materialize genomes for every strain in the cohort.

    Wave 2.5 hardening C4: uses `CandidateStrain.assembly_accession` (NCBI
    Datasets accession like `GCF_000005845.2`) for the refseq fetch, NOT
    `strain_id` (which is the BV-BRC genome id). Strains without an
    assembly_accession raise CohortConstructionError — they cannot be fetched.

    Returns mapping strain_id -> per-accession cache directory. The map key
    stays as the BV-BRC strain_id so downstream code (cohort.parquet,
    embedding cache) uses one consistent ID space.

    Uses `ThreadPoolExecutor(max_workers=4)` by default (NCBI Datasets API
    politeness limit).
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from dna_decode.data import refseq

    out: dict[str, Path] = {}
    if not cohort.strains:
        return out

    missing_accession = [s.strain_id for s in cohort.strains if not s.assembly_accession]
    if missing_accession:
        raise CohortConstructionError(
            f"Cannot download: {len(missing_accession)} strains missing "
            f"assembly_accession (sample: {missing_accession[:3]}). Resolve via "
            f"NCBI assembly summary lookup before calling download_cohort_genomes."
        )

    def _fetch(strain: CandidateStrain) -> tuple[str, Path]:
        path = refseq.download_genome(strain.assembly_accession, cache_root)
        return strain.strain_id, path

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(_fetch, s): s.strain_id for s in cohort.strains}
        for fut in as_completed(futures):
            sid, path = fut.result()
            out[sid] = path

    return out
