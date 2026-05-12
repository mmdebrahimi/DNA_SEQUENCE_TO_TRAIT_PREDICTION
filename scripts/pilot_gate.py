"""CLI for Step 0.5 — Real-data pilot gate (HARD).

Usage:
    python -m scripts.pilot_gate \
        [--drugs cipro,ceftriaxone,tetracycline] \
        [--target-per-drug 150] \
        [--ast-tsv <local-bvbrc-ast-tsv-path>]

HARD-gate exit semantics:
    0 = GO; 1 = NO-GO; 2 = PilotGateError; 3 = NotImplementedError (no AST source).

Provide a BV-BRC AST TSV via --ast-tsv (or BVBRC_AST_TSV env var, or
config.bvbrc_ast.local_tsv_path). Without one, the live-API path raises
NotImplementedError → exit 3. Download an AST TSV from ftp.bvbrc.org.

Writes report to data/processed/pilot_report.md.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dna_decode.data.pilot import (
    CohortSelectionCriteria,
    PilotGateError,
    run_pilot_gate,
    write_pilot_report,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase 1 real-data pilot gate (Step 0.5). HARD gate: exit 0 only when "
        "per-drug counts meet target.",
    )
    parser.add_argument(
        "--drugs",
        default="ciprofloxacin,ceftriaxone,tetracycline",
        help="Comma-separated drug list (default: Phase 1 trio).",
    )
    parser.add_argument(
        "--target-per-drug",
        type=int,
        default=150,
        help="Minimum strains per drug after all filters (default: 150).",
    )
    parser.add_argument(
        "--three-drug-intersection-target",
        type=int,
        default=75,
        help="Minimum strains in 3-drug intersection (default: 75).",
    )
    parser.add_argument(
        "--config",
        default="config/datasources.yaml",
        help="Path to datasources config (default: config/datasources.yaml).",
    )
    parser.add_argument(
        "--output",
        default="data/processed/pilot_report.md",
        help="Output report path (default: data/processed/pilot_report.md).",
    )
    parser.add_argument(
        "--ast-tsv",
        default=None,
        help="Path to a local BV-BRC AST TSV. If omitted, falls back to "
        "BVBRC_AST_TSV env var or config.bvbrc_ast.local_tsv_path. If neither "
        "is set, the live-API path raises NotImplementedError (exit 3).",
    )
    args = parser.parse_args(argv)

    drugs = tuple(d.strip() for d in args.drugs.split(","))
    criteria = CohortSelectionCriteria(
        target_per_drug=args.target_per_drug,
        three_drug_intersection_target=args.three_drug_intersection_target,
    )

    try:
        report = run_pilot_gate(
            drugs=drugs,
            criteria=criteria,
            config_path=args.config,
            ast_tsv_path=args.ast_tsv,
        )
    except PilotGateError as e:
        print(f"PILOT GATE FAILED (cannot run): {e}", file=sys.stderr)
        return 2
    except NotImplementedError as e:
        print(
            f"PILOT GATE SCAFFOLD: real API integration deferred. {e}",
            file=sys.stderr,
        )
        print("Implement fetch_bvbrc_drug_counts + fetch_ncbi_assembly_quality first.", file=sys.stderr)
        return 3

    write_pilot_report(report, args.output)
    print(f"Pilot report written: {args.output}")
    print(f"Verdict: {report.go_no_go}")
    if report.go_no_go == "NO-GO":
        for reason in report.failure_reasons:
            print(f"  - {reason}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
