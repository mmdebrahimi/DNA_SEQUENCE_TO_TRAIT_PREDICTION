"""Ceftriaxone MIC-tier audit for the gate-B cohort.

Reads the raw BV-BRC AST CSV, rejoins per-strain cef MIC rows onto the existing
cef cohort, and classifies each strain into a MIC confidence tier using the
shared breakpoints in `dna_decode.data.mic_tiers`.

Output:
- wiki/ceftriaxone_mic_audit_<date>.{md,json}
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from datetime import date as _date
from pathlib import Path
from statistics import median

import pandas as pd

from dna_decode.data.ast_data import load_bvbrc_ast
from dna_decode.data.cohort import load_cohort
from dna_decode.data.mic_tiers import breakpoints_for, classify_tier


DEFAULT_AST_CSV = Path(
    os.environ.get(
        "BVBRC_AST_TSV",
        "C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv",
    )
)


def _tier_detail(mics: list[float], breakpoints: dict) -> dict:
    if not mics:
        return {}
    med = median(mics)
    clsi_r = breakpoints["clsi_r"]
    clsi_s = breakpoints["clsi_s"]
    eucast_r = breakpoints.get("eucast_r")
    eucast_s = breakpoints.get("eucast_s")
    detail = {
        "median_mic": med,
        "min_mic": min(mics),
        "max_mic": max(mics),
        "n_mic_rows": len(mics),
        "clsi_call": "R" if med >= clsi_r else ("S" if med <= clsi_s else "I"),
        "distance_to_clsi_r": med / clsi_r if clsi_r else None,
        "distance_to_clsi_s": med / clsi_s if clsi_s else None,
    }
    if eucast_r is not None and eucast_s is not None:
        detail["eucast_call"] = (
            "R" if med >= eucast_r else ("S" if med <= eucast_s else "I")
        )
    else:
        detail["eucast_call"] = None
    return detail


def _cohort_ids_for_drug(cohort, drug: str, *, all_labeled: bool = False) -> tuple[set[str], str]:
    drug_lower = drug.lower()
    if all_labeled:
        return {s.strain_id for s in cohort.strains if drug_lower in s.ast_labels}, "labeled"
    pool_ids = cohort.per_drug_strain_ids.get(drug_lower, [])
    if pool_ids:
        return set(pool_ids), "pool"
    return {s.strain_id for s in cohort.strains if drug_lower in s.ast_labels}, "labeled"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cohort",
        type=Path,
        default=Path("data/processed/gate_b_cohort.parquet"),
    )
    parser.add_argument("--ast-csv", type=Path, default=DEFAULT_AST_CSV)
    parser.add_argument("--drug", default="ceftriaxone")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--all-labeled",
        action="store_true",
        help="Audit all labeled strains instead of the saved per-drug in_pool_<drug> subset.",
    )
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/{args.drug.lower()}_mic_audit_{_date.today().isoformat()}.md")

    if not args.ast_csv.exists():
        print(f"[cef_mic_audit] raw AST CSV not found: {args.ast_csv}", file=sys.stderr)
        return 2

    cohort = load_cohort(args.cohort)
    drug_lower = args.drug.lower()
    cohort_ids, cohort_scope = _cohort_ids_for_drug(cohort, args.drug, all_labeled=args.all_labeled)
    strains_with_label = [s for s in cohort.strains if s.strain_id in cohort_ids]
    label_by_id = {s.strain_id: int(s.ast_labels[drug_lower]) for s in strains_with_label}
    accession_by_id = {s.strain_id: s.assembly_accession for s in strains_with_label}
    mlst_by_id = {s.strain_id: s.mlst for s in strains_with_label}
    breakpoints = breakpoints_for(drug_lower)

    ast = load_bvbrc_ast(args.ast_csv)
    ast = ast[ast["antibiotic"] == drug_lower]
    ast = ast[ast["strain_id"].isin(cohort_ids)]

    per_strain_rows: dict[str, list[dict]] = defaultdict(list)
    for _, row in ast.iterrows():
        sid = row["strain_id"]
        per_strain_rows[sid].append(
            {
                "mic_value": row["mic_value"],
                "mic_units": row["mic_units"],
                "method": row["measurement_method"],
                "source": row["source"],
                "phenotype": row["susceptibility_label"],
            }
        )

    audit_rows: list[dict] = []
    r_tier_counts: dict[str, int] = defaultdict(int)
    s_tier_counts: dict[str, int] = defaultdict(int)
    for sid in sorted(cohort_ids):
        rows = per_strain_rows.get(sid, [])
        mics = [float(r["mic_value"]) for r in rows if pd.notna(r["mic_value"])]
        calls = {str(r["phenotype"]).upper() for r in rows if r["phenotype"]}
        tier = classify_tier(mics, calls, breakpoints)
        detail = _tier_detail(mics, breakpoints)
        label = label_by_id[sid]
        if label == 1:
            r_tier_counts[tier] += 1
        else:
            s_tier_counts[tier] += 1
        audit_rows.append(
            {
                "strain_id": sid,
                "accession": accession_by_id[sid],
                "mlst": mlst_by_id[sid],
                "cohort_binary_label": label,
                "cohort_binary_R_or_S": "R" if label == 1 else "S",
                "n_ast_rows": len(rows),
                "n_mic_rows": len(mics),
                "tier": tier,
                "detail": detail,
                "rows": rows,
            }
        )

    decisive_r = [
        r["strain_id"] for r in audit_rows if r["cohort_binary_label"] == 1 and r["tier"] == "HIGH_R"
    ]
    decisive_s = [
        r["strain_id"] for r in audit_rows if r["cohort_binary_label"] == 0 and r["tier"] == "HIGH_S"
    ]
    verdict = (
        "CLEAN"
        if (r_tier_counts.get("HIGH_R", 0) / max(1, sum(r_tier_counts.values()))) >= 0.7
        and r_tier_counts.get("NO_MIC", 0) <= 2
        else "MIXED"
        if (r_tier_counts.get("HIGH_R", 0) / max(1, sum(r_tier_counts.values()))) >= 0.4
        else "NOISY"
    )

    json_path = args.output.with_suffix(".json")
    payload = {
        "cohort_path": str(args.cohort),
        "ast_csv": str(args.ast_csv),
        "drug": args.drug,
        "breakpoints": breakpoints,
        "n_cohort": len(strains_with_label),
        "r_tier_counts": dict(r_tier_counts),
        "s_tier_counts": dict(s_tier_counts),
        "decisive_R_ids": decisive_r,
        "decisive_S_ids": decisive_s,
        "verdict": verdict,
        "per_strain": audit_rows,
    }
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    lines = [
        f"# Ceftriaxone AST/MIC audit ({_date.today().isoformat()})",
        "",
        f"- Cohort: `{args.cohort}`",
        f"- Source: `{args.ast_csv}`",
        f"- N {cohort_scope} strains: `{len(strains_with_label)}`",
        f"- Verdict: `{verdict}`",
        "",
        "| tier | R strains | S strains |",
        "|---|---:|---:|",
    ]
    for tier in sorted(set(r_tier_counts) | set(s_tier_counts)):
        lines.append(f"| {tier} | {r_tier_counts.get(tier, 0)} | {s_tier_counts.get(tier, 0)} |")
    lines.extend(
        [
            "",
            "| strain_id | accession | label | tier | median MIC | n_rows | n_mic | mlst |",
            "|---|---|---|---|---:|---:|---:|---|",
        ]
    )
    for row in sorted(audit_rows, key=lambda x: (x["cohort_binary_label"], x["tier"], x["strain_id"])):
        med = row["detail"].get("median_mic")
        med_str = f"{med:.3f}" if isinstance(med, (int, float)) else "-"
        lines.append(
            f"| {row['strain_id']} | {row['accession']} | {row['cohort_binary_R_or_S']} | "
            f"{row['tier']} | {med_str} | {row['n_ast_rows']} | {row['n_mic_rows']} | {row['mlst']} |"
        )
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[cef_mic_audit] wrote {args.output}")
    print(f"[cef_mic_audit] wrote {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
