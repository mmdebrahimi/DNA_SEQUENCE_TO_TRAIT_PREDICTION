"""Build a Gate B mini-cohort (12 strains, 6 R + 6 S, single drug) from a larger cohort parquet.

Per the /brainstorm Option C synthesis (2026-05-13): full 67-strain real-NT
populate took >2 hr wallclock and is at risk of mid-run disconnects on
external storage. A 12-strain cipro-balanced mini-cohort matches the
original Gate B framing ("12-20 strain mini-cohort dry-run"), finishes in
~30 min, and limits blast radius if the embedding storage drops mid-run.

Strain selection prioritizes high N50 + low contig_count (best assembly
quality) within each R/S class — preserves data-quality signal over
random sampling.

Usage:
    uv run python scripts/build_mini_cohort.py \
        --source data/processed/gate_b_cohort.parquet \
        --output data/processed/gate_b_mini_cohort.parquet \
        --drug ciprofloxacin \
        --per-class 6
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dna_decode.data.cohort import CandidateStrain, StrainCohort, load_cohort, save_cohort


def _assembly_quality_score(strain: CandidateStrain) -> tuple[int, int]:
    """Higher-quality assemblies sort first: low contig_count + high n50.

    Returns a tuple usable with sorted(reverse=False): smaller is better.
    Strain with contig_count=10, n50=300000 ranks above contig_count=200, n50=80000.
    """
    return (strain.contig_count, -strain.n50)


def select_mini_cohort(
    source: StrainCohort, drug: str, per_class: int
) -> StrainCohort:
    """Pick `per_class` highest-quality strains from each of R / S for `drug`."""
    drug_lower = drug.lower()
    if drug_lower not in source.per_drug_strain_ids:
        raise ValueError(
            f"drug {drug!r} not in source.per_drug_strain_ids "
            f"({list(source.per_drug_strain_ids.keys())})"
        )

    drug_strain_ids = set(source.per_drug_strain_ids[drug_lower])
    drug_strains = [s for s in source.strains if s.strain_id in drug_strain_ids]

    resistant = [s for s in drug_strains if s.ast_labels.get(drug_lower) == 1]
    susceptible = [s for s in drug_strains if s.ast_labels.get(drug_lower) == 0]

    resistant.sort(key=_assembly_quality_score)
    susceptible.sort(key=_assembly_quality_score)

    if len(resistant) < per_class:
        raise ValueError(
            f"only {len(resistant)} resistant strains for {drug!r}; need {per_class}"
        )
    if len(susceptible) < per_class:
        raise ValueError(
            f"only {len(susceptible)} susceptible strains for {drug!r}; need {per_class}"
        )

    picked = resistant[:per_class] + susceptible[:per_class]
    picked_ids = [s.strain_id for s in picked]

    return StrainCohort(
        strains=picked,
        per_drug_strain_ids={drug_lower: picked_ids},
        three_drug_intersection=picked_ids,  # single-drug context; treat as full intersection
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a balanced mini-cohort from a source cohort parquet."
    )
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--drug", default="ciprofloxacin", help="Drug for R/S balance + selection"
    )
    parser.add_argument(
        "--per-class", type=int, default=6, help="Strains per R / S class (mini = 6 → 12)"
    )
    args = parser.parse_args(argv)

    if not args.source.exists():
        print(f"[build_mini_cohort] source not found: {args.source}", file=sys.stderr)
        return 2

    source = load_cohort(args.source)
    print(f"[build_mini_cohort] source: {len(source)} strains")
    try:
        mini = select_mini_cohort(source, args.drug, args.per_class)
    except ValueError as e:
        print(f"[build_mini_cohort] {e}", file=sys.stderr)
        return 1

    save_cohort(mini, args.output)
    drug_lower = args.drug.lower()
    r_count = sum(1 for s in mini.strains if s.ast_labels.get(drug_lower) == 1)
    s_count = sum(1 for s in mini.strains if s.ast_labels.get(drug_lower) == 0)
    print(
        f"[build_mini_cohort] saved {args.output}: {len(mini)} strains "
        f"({r_count} R / {s_count} S for {args.drug})"
    )
    # Per-strain summary
    for s in mini.strains:
        print(
            f"  {s.strain_id}: contig_count={s.contig_count} n50={s.n50:,} "
            f"mlst={s.mlst!r} {args.drug}={s.ast_labels.get(drug_lower)}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
