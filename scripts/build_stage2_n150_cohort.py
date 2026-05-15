"""Phase A.0 of `plans/Stage2_N150_Prep_Plan.md` — build the Stage 2 N=150 cipro cohort.

`scripts/audit_cohort.py` only AUDITS existing cohorts; the actual cohort-construction
path is `dna_decode.data.cohort.candidates_from_bvbrc_ast` + `build_cohort` + `save_cohort`
called directly with explicit thresholds. This script wires that up with CLI flags.

Hard-fails (no silent imbalance, no silent under-target):
  - If the broth-microdilution-filtered cipro pool can't hit `--target-total` ± `--balance-slack`.
  - If R/S balance falls outside `(per_class - balance_slack, per_class + balance_slack)`.
  - If any selected strain lacks MLST (matches Stage 1's loud-MLST invariant).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dna_decode.data.ast_data import load_bvbrc_ast
from dna_decode.data.bvbrc_genome import load_bvbrc_genome_metadata
from dna_decode.data.cohort import (
    CohortSelectionCriteria,
    StrainCohort,
    _filter_by_assembly_quality,
    _mlst_balanced_selection,
    candidates_from_bvbrc_ast,
    save_cohort,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Stage 2 N=150 cipro cohort")
    parser.add_argument("--ast-tsv", type=Path, required=True,
                        help="Path to BV-BRC AST TSV/CSV (genome AMR phenotype table)")
    parser.add_argument("--assembly-metadata-csv", type=Path, required=True,
                        help="Path to BV-BRC genome metadata CSV (BVBRC_genome.csv)")
    parser.add_argument("--drug", default="ciprofloxacin",
                        help="Drug name (lowercase) for cohort construction")
    parser.add_argument("--target-total", type=int, default=150,
                        help="Target total strain count (default 150 for Stage 2)")
    parser.add_argument("--per-class", type=int, default=75,
                        help="Target R-count and S-count separately (default 75/75)")
    parser.add_argument("--balance-slack", type=int, default=10,
                        help="Acceptable per-class deviation (default +/-10, i.e. 65-85 each)")
    parser.add_argument("--n50-min", type=int, default=50_000,
                        help="Minimum N50 (repo default 50K; loosen for more candidates)")
    parser.add_argument("--contig-count-max", type=int, default=500,
                        help="Maximum contig count (repo default 500)")
    parser.add_argument("--output", type=Path,
                        default=Path("data/processed/stage2_n150_cipro_cohort.parquet"),
                        help="Output parquet path")
    args = parser.parse_args(argv)

    drug_lower = args.drug.lower()

    print(f"[build] loading AST TSV from {args.ast_tsv}")
    ast = load_bvbrc_ast(args.ast_tsv)
    print(f"[build] AST rows: {len(ast)}")

    print(f"[build] loading assembly metadata from {args.assembly_metadata_csv}")
    metadata = load_bvbrc_genome_metadata(args.assembly_metadata_csv)
    print(f"[build] metadata strains: {len(metadata)}")

    print(f"[build] constructing candidate list ({drug_lower}-labeled, AST table joined with metadata)")
    candidates = candidates_from_bvbrc_ast(ast_table=ast, assembly_metadata=metadata)
    print(f"[build] candidates (pre-MLST-filter): {len(candidates)}")

    # Pre-filter MLST-missing candidates so build_cohort's MLST-balanced selection
    # doesn't pick them (caught 2 such strains on first run; build_cohort doesn't
    # enforce MLST presence by default). Loud-fail at selection-time would force
    # re-run; pre-filter is cheaper.
    before = len(candidates)
    candidates = [c for c in candidates if getattr(c, "mlst", None) and str(c.mlst).strip() not in ("", "None")]
    print(f"[build] candidates after MLST-presence filter: {len(candidates)} (dropped {before - len(candidates)})")

    criteria = CohortSelectionCriteria(
        target_per_drug=args.target_total,
        three_drug_intersection_target=0,  # Stage 2 is single-drug
        assembly_contig_count_max=args.contig_count_max,
        assembly_n50_min=args.n50_min,
    )

    # Apply assembly-quality filter (same logic as build_cohort)
    quality_filtered = _filter_by_assembly_quality(candidates, criteria)
    print(f"[build] candidates after assembly-quality filter "
          f"(n50>={args.n50_min}, contig<={args.contig_count_max}): {len(quality_filtered)}")

    # Restrict to strains labeled for the target drug
    drug_labeled = [c for c in quality_filtered if drug_lower in c.ast_labels]
    print(f"[build] candidates labeled for {drug_lower}: {len(drug_labeled)}")

    # LABEL-STRATIFIED selection: pick per_class from R and per_class from S separately,
    # each with MLST-diversity balancing. Replaces build_cohort's diversity-only selection
    # which left available R strains on the table (Stage 2 R-ceiling was 49 with default
    # algo; ~77 R strains were actually available -- see scripts/diagnose_bvbrc_mlst_gaps.py
    # 2026-05-14 finding).
    r_pool = [c for c in drug_labeled if c.ast_labels.get(drug_lower) == 1]
    s_pool = [c for c in drug_labeled if c.ast_labels.get(drug_lower) == 0]
    print(f"[build] R pool: {len(r_pool)} candidates; S pool: {len(s_pool)} candidates")

    selected_r = _mlst_balanced_selection(r_pool, args.per_class)
    selected_s = _mlst_balanced_selection(s_pool, args.per_class)
    print(f"[build] selected: {len(selected_r)} R + {len(selected_s)} S = {len(selected_r) + len(selected_s)}")

    # Compose StrainCohort manually (bypassing build_cohort's drug-intersection logic
    # which is N/A for Stage 2's single-drug case)
    selected_strains = selected_r + selected_s
    strain_ids_list = [s.strain_id for s in selected_strains]
    cohort = StrainCohort(
        strains=selected_strains,
        per_drug_strain_ids={drug_lower: strain_ids_list},
        three_drug_intersection=[],
    )
    print(f"[build] cohort size: {len(cohort)} strains")

    # Balance + MLST loud-checks
    drug_strain_ids = cohort.per_drug_strain_ids.get(drug_lower, [])
    r_count = 0
    s_count = 0
    missing_mlst: list[str] = []
    for sid in drug_strain_ids:
        s = cohort.strain_by_id(sid)
        if s is None:
            continue
        lbl = s.ast_labels.get(drug_lower)
        if lbl == 1:
            r_count += 1
        elif lbl == 0:
            s_count += 1
        mlst = getattr(s, "mlst", None)
        if mlst is None or str(mlst).strip() in ("", "None"):
            missing_mlst.append(sid)

    print(f"[build] balance: {r_count}R / {s_count}S (target {args.per_class}/{args.per_class}, slack ±{args.balance_slack})")

    if missing_mlst:
        print(f"[build] ERROR: {len(missing_mlst)} selected strains missing MLST: {missing_mlst[:10]}",
              file=sys.stderr)
        return 2

    if r_count < args.per_class - args.balance_slack or r_count > args.per_class + args.balance_slack:
        print(f"[build] ERROR: R count {r_count} outside acceptable range "
              f"[{args.per_class - args.balance_slack}, {args.per_class + args.balance_slack}]; "
              f"escalate to user (loosen thresholds or accept imbalance)", file=sys.stderr)
        return 3
    if s_count < args.per_class - args.balance_slack or s_count > args.per_class + args.balance_slack:
        print(f"[build] ERROR: S count {s_count} outside acceptable range "
              f"[{args.per_class - args.balance_slack}, {args.per_class + args.balance_slack}]; "
              f"escalate to user", file=sys.stderr)
        return 3

    print(f"[build] saving to {args.output}")
    save_cohort(cohort, args.output)
    print(f"[build] DONE — cohort with {len(cohort)} strains ({r_count}R / {s_count}S) at {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
