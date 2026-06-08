"""Biocide/disinfectant resistance — CLI (console `dna-disinfinder`; also `dna-decode disinfinder`).

Genome FASTA -> acquired biocide-resistance genes (DisinFinder allele DB via real blastn). Quaternary-
ammonium (qacA/qacB/...) + formaldehyde (formA) etc. — hospital infection-control relevant; these genes
often share plasmids with AMR (pair with dna-coloc). Deterministic curated-DB caller; offline-safe (exit 3).

    dna-disinfinder assembly.fna --sample-id MY_STRAIN
    dna-disinfinder X.fna --db data/disinfinder_db/disinfectants.fsa
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.disinfinder.runner import call_disinfectant_genes

DEFAULT_DB = "data/disinfinder_db/disinfectants.fsa"
_DB_URL = "https://bitbucket.org/genomicepidemiology/disinfinder_db/raw/HEAD/disinfectants.fsa"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-disinfinder",
                                 description="Deterministic biocide/disinfectant resistance genes (DisinFinder)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db", default=DEFAULT_DB, help=f"DisinFinder allele .fsa (default {DEFAULT_DB})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--identity", type=float, default=90.0, help="min %% identity (default 90)")
    ap.add_argument("--coverage", type=float, default=60.0, help="min %% coverage (default 60)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    db_sha = hashlib.sha256(Path(args.db).read_bytes()).hexdigest()[:16] if Path(args.db).exists() else None

    res = call_disinfectant_genes(args.fasta, args.db, identity_threshold=args.identity,
                                  coverage_threshold=args.coverage)
    rec = {
        "sample_id": sample_id, "trait": "biocide_resistance_genes",
        "analysis_date": datetime.date.today().isoformat(), "schema": "disinfectant-gene-call-v0",
        "status": res["status"], "genes": [g["gene"] for g in res.get("genes", [])],
        "gene_detail": res.get("genes", []),
        "caller": {"name": "dna_decode-disinfinder-blastn-v0", "method": res.get("method"),
                   "source": "DisinFinder curated allele DB", "caller_is_independent_baseline": False},
        "caveat": ("Acquired biocide/disinfectant-resistance genes (qac/form/...) via blastn over the "
                   "DisinFinder DB — faithful to the tool, not an independent baseline. Gene presence != "
                   "clinical disinfectant failure. Often plasmid-borne (pair with dna-coloc). NOT a clinical tool."),
        "provenance": {"db": str(args.db), "db_name": "disinfinder", "db_sha256_16": db_sha,
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage},
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: biocide/disinfectant resistance genes")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
            if res.get("tool") == "db":
                print(f"  download the DB: curl -sSL {_DB_URL} -o {args.db}")
        else:
            print(f"GENES: {len(res['genes'])} acquired")
            for g in res["genes"]:
                print(f"  {g['gene']:14} {g['percent_identity']}% id / {g['percent_coverage']}% cov")
            if not res["genes"]:
                print("  (no acquired biocide-resistance genes detected at threshold)")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
