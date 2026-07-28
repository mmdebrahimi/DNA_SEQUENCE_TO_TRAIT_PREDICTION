"""`dna-decode metabolic` / `dna-metabolic` — carbon-source utilization from gene/operon presence.

    dna-decode metabolic --source lactose --genes lacZ,lacY
    dna-decode metabolic --source citrate --genes citD,citE,citF,citT --condition aerobic   # -> Cit- (anchor)
    dna-decode metabolic --source lactose --feature-table GCF_..._feature_table.txt.gz
    dna-decode metabolic --list

The DETERMINISTIC E. coli carbon-catabolism decoder: utilizes iff (catabolic enzymes present) AND (a
transporter present) AND (transporter expressed under the O2 condition). The uptake-gate is what a naive
AMR-style has-the-genes rule misses — the citrate aerobic Cit- anchor. Label-faithful (EcoCyc/textbook),
offline, no deps. Tier: KNOWLEDGE_BASELINE. NOT a clinical tool.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

_SCOPE = ("scope: E. coli carbon-catabolism, uptake-gated, presence-based; calls can/cannot DIRECTION not "
          "growth rate; reads gene presence not sequence integrity; validated vs measured K-12 phenotypes. "
          "NOT clinical.")


def _feature_table_symbols(path: str) -> set[str]:
    op = gzip.open if str(path).endswith(".gz") else open
    syms: set[str] = set()
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            isym = header.index("symbol")
        except ValueError:
            return syms
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > isym and p[0] == "CDS" and p[isym]:
                syms.add(p[isym])
    return syms


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode metabolic",
        description="E. coli carbon-source utilization from gene presence (deterministic, uptake-gated).",
        epilog=_SCOPE)
    ap.add_argument("--source", help="carbon source (e.g. lactose, L-arabinose, citrate, D-glucose)")
    ap.add_argument("--genes", help="comma-separated PRESENT gene symbols (e.g. lacZ,lacY)")
    ap.add_argument("--feature-table", help="NCBI feature_table.txt(.gz): read present symbols genome-wide")
    ap.add_argument("--condition", default="aerobic", choices=("aerobic", "anaerobic"),
                    help="O2 condition (default aerobic) — matters for uptake-gated sources like citrate")
    ap.add_argument("--list", action="store_true", help="list catalogued carbon sources + their genes")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.metabolic.carbon_catalog import (
        CARBON_SOURCES, MetabolicInputError, call_carbon_utilization, genes_for)

    if args.list:
        rows = {s: genes_for(s) for s in CARBON_SOURCES}
        if args.json:
            print(json.dumps({"substrates": rows, "scope": _SCOPE}, indent=2))
        else:
            for s, g in rows.items():
                print(f"{s}: {', '.join(g)}")
            print(_SCOPE)
        return 0

    if not args.source:
        print("error: give --source (with --genes or --feature-table), or --list", file=sys.stderr)
        return 2

    if args.feature_table:
        if not Path(args.feature_table).exists():
            print(f"error: cannot read --feature-table: {args.feature_table}", file=sys.stderr)
            return 2
        present = _feature_table_symbols(args.feature_table)
        if not present:
            print("error: no CDS rows with a 'symbol' column found in the feature table", file=sys.stderr)
            return 2
    elif args.genes is not None:
        present = {g.strip() for g in args.genes.split(",") if g.strip()}
    else:
        print("error: give --genes or --feature-table with --source", file=sys.stderr)
        return 2

    try:
        c = call_carbon_utilization(args.source, present, condition=args.condition)
    except MetabolicInputError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    d = c.as_dict(); d["scope"] = _SCOPE; d["tier"] = "KNOWLEDGE_BASELINE"
    if args.json:
        print(json.dumps(d, indent=2))
    else:
        verb = {"utilizes": "UTILIZES", "cannot_utilize": "CANNOT UTILIZE", "ABSTAIN": "ABSTAIN"}[c.capability]
        print(f"{c.substrate} ({c.condition}): {verb}  (confidence {c.confidence})")
        for n in c.notes:
            print(f"  - {n}")
        print(_SCOPE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
