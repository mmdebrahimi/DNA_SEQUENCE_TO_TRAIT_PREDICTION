"""Slim BV-BRC AST CSV to only the rows for the cef-pool cohort strains.

The full BVBRC_genome_amr.csv is ~432 MB — too large to ship via Gmail (25 MB cap).
Codex on Precision 7780 needs the cef-pool strains' MIC rows to run
`scripts/cef_mic_audit.py` for the audit-aware packet. This script extracts only
those rows so the result fits in a Gmail-able bundle.

Run from repo root:

    uv run python scripts/slim_ast_csv_for_bundle.py \
        --ast-csv "C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv" \
        --cohort data/processed/gate_b_cohort.parquet \
        --output "C:/Users/Farshad/AppData/Local/Temp/bvbrc_ast_slim_cohort_2026-05-26.csv"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ast-csv", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--include-all-cohort",
        action="store_true",
        default=True,
        help="Include all 67 cohort strains' rows (default), not just cef-pool.",
    )
    args = parser.parse_args(argv)

    if not args.ast_csv.exists():
        print(f"[slim] AST CSV missing: {args.ast_csv}", file=sys.stderr)
        return 2
    if not args.cohort.exists():
        print(f"[slim] cohort missing: {args.cohort}", file=sys.stderr)
        return 2

    cohort_df = pd.read_parquet(args.cohort)
    if args.include_all_cohort:
        strain_ids = set(cohort_df["strain_id"].astype(str).tolist())
        label = "all-67-cohort"
    else:
        strain_ids = set(
            cohort_df[cohort_df["in_pool_ceftriaxone"]]["strain_id"].astype(str).tolist()
        )
        label = "cef-pool"
    print(f"[slim] keeping {len(strain_ids)} {label} strain_ids")

    # Stream-read the 432 MB CSV in chunks; keep only rows where Genome ID matches.
    keeper_chunks: list[pd.DataFrame] = []
    rows_seen = 0
    for chunk in pd.read_csv(args.ast_csv, chunksize=50_000, dtype=str):
        rows_seen += len(chunk)
        keep = chunk[chunk["Genome ID"].isin(strain_ids)]
        if len(keep):
            keeper_chunks.append(keep)

    if not keeper_chunks:
        print("[slim] no matching rows found — wrong column or strain_id format?", file=sys.stderr)
        return 3
    out = pd.concat(keeper_chunks, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    in_size = args.ast_csv.stat().st_size
    out_size = args.output.stat().st_size
    print(f"[slim] scanned {rows_seen:,} rows; kept {len(out):,}")
    print(f"[slim] input:  {in_size:>15,} B  ({in_size/1024/1024:.1f} MB)")
    print(f"[slim] output: {out_size:>15,} B  ({out_size/1024/1024:.2f} MB)")
    print(f"[slim] wrote {args.output}")
    matched = out["Genome ID"].nunique()
    print(f"[slim] distinct Genome IDs matched: {matched} / {len(strain_ids)}")
    missing = strain_ids - set(out["Genome ID"].astype(str).tolist())
    if missing:
        print(f"[slim] WARNING: {len(missing)} cohort strain_ids had NO AST rows: {sorted(missing)[:5]}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
