"""Acquired-AMR-gene detection — in-package CLI (console `dna-resfinder`; also `dna-decode resfinder`).

Genome FASTA -> acquired AMR genes present (ResFinder allele DB via real blastn), grouped by antibiotic
class. An INDEPENDENT second opinion vs the AMRFinder-based `dna-amr` engine — use the concordance as a
cross-tool check. Deterministic curated-DB caller; offline-safe (status "unavailable", exit 3).

    dna-resfinder assembly.fna --sample-id MY_STRAIN          # all class DBs in data/resfinder_db/
    dna-resfinder X.fna --db data/resfinder_db/beta-lactam.fsa --drug-class beta-lactam   # single class
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from dna_decode.resfinder.runner import call_resistance_genes

DEFAULT_DB_DIR = "data/resfinder_db"
_DB_BASE = "https://bitbucket.org/genomicepidemiology/resfinder_db/raw/HEAD"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-resfinder",
                                 description="Deterministic acquired-AMR-gene detection (ResFinder via blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db", type=Path, default=None,
                    help="a single ResFinder class .fsa; omit to scan every *.fsa in --db-dir")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR, help=f"dir of ResFinder class .fsa (default {DEFAULT_DB_DIR})")
    ap.add_argument("--drug-class", default=None, help="label for a single --db (e.g. beta-lactam)")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--identity", type=float, default=90.0, help="min %% identity (ResFinder default 90)")
    ap.add_argument("--coverage", type=float, default=60.0, help="min %% coverage (default 60)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem

    # one class DB, or every *.fsa in the dir (drug_class = filename stem)
    if args.db is not None:
        targets = [(args.db, args.drug_class)]
    else:
        d = Path(args.db_dir)
        targets = [(p, p.stem) for p in sorted(d.glob("*.fsa"))] if d.exists() else []

    all_genes: list[dict] = []
    status = "ok"
    reason = None
    if not targets:
        status, reason = "unavailable", f"no ResFinder class .fsa found (db {args.db or args.db_dir})"
    for db, cls in targets:
        res = call_resistance_genes(args.fasta, db, drug_class=cls,
                                    identity_threshold=args.identity, coverage_threshold=args.coverage)
        if res["status"] != "ok":
            status, reason = "unavailable", res.get("reason")
            break
        all_genes.extend(res["genes"])
    all_genes.sort(key=lambda r: (-r["percent_coverage"], r["gene"]))

    rec = {
        "sample_id": sample_id, "trait": "acquired_amr_genes",
        "analysis_date": datetime.date.today().isoformat(), "schema": "resfinder-gene-call-v0",
        "status": status, "genes": [g["gene"] for g in all_genes],
        "drug_classes": sorted({g["drug_class"] for g in all_genes if g["drug_class"]}),
        "gene_detail": all_genes,
        "caller": {"name": "dna_decode-resfinder-blastn-v0", "method": "resfinder_blastn_v0",
                   "source": "ResFinder curated allele DB", "caller_is_independent_baseline": True},
        "caveat": ("INDEPENDENT acquired-gene caller (ResFinder DB) — use as a cross-tool concordance check "
                   "vs dna-amr (AMRFinder DB), NOT as its replacement. Detects ACQUIRED genes only (no "
                   "chromosomal point mutations / efflux-expression phenotypes). NOT a clinical tool."),
        "provenance": {"db": str(args.db or args.db_dir), "db_name": "resfinder",
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage,
                       "n_class_dbs": len(targets)},
    }
    if status != "ok":
        rec["reason"] = reason
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: acquired AMR genes (ResFinder)")
        if status != "ok":
            print(f"STATUS: {status} - {reason}")
            print(f"  download class DBs into {args.db_dir}/ from {_DB_BASE}/<class>.fsa")
        else:
            print(f"GENES: {len(all_genes)} acquired  | classes: {', '.join(rec['drug_classes']) or '(none)'}")
            for g in all_genes:
                print(f"  {g['gene']:22} [{g['drug_class']}]  {g['percent_identity']}% id / "
                      f"{g['percent_coverage']}% cov")
            if not all_genes:
                print("  (no acquired AMR genes detected at threshold)")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if status == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
