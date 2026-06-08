"""Plasmid replicon typing — in-package CLI (console entry `dna-plasmid`; also `dna-decode plasmid`).

Genome FASTA -> the plasmid Inc replicons present (PlasmidFinder allele DB via real blastn) + provenance.
Deterministic curated-DB caller (sibling of dna-amr / dna-pathotype). Composes with dna-amr: which mobile
element likely carries the resistance. Offline-safe: no blastn / no DB -> status "unavailable" (exit 3),
never a crash.

    dna-plasmid path/to/assembly.fna --sample-id MY_STRAIN
    dna-plasmid X.fna --db data/plasmidfinder_db/enterobacteriales.fsa

NOT a clinical tool. The DB is the curated PlasmidFinder allele set; this caller is faithful to
PlasmidFinder's method (blastn over the allele DB, identity 95% / coverage 60%), not an independent baseline.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.plasmid.runner import call_replicons

DEFAULT_DB = "data/plasmidfinder_db/enterobacteriales.fsa"
_DB_URL = "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/raw/HEAD/enterobacteriales.fsa"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-plasmid",
                                 description="Deterministic plasmid Inc-replicon typing (PlasmidFinder via blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db", default=DEFAULT_DB, help=f"PlasmidFinder allele .fsa (default {DEFAULT_DB})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--identity", type=float, default=95.0, help="min %% identity (PlasmidFinder default 95)")
    ap.add_argument("--coverage", type=float, default=60.0, help="min %% coverage (PlasmidFinder default 60)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    db_sha = hashlib.sha256(Path(args.db).read_bytes()).hexdigest()[:16] if Path(args.db).exists() else None

    res = call_replicons(args.fasta, args.db, identity_threshold=args.identity,
                         coverage_threshold=args.coverage)
    rec = {
        "sample_id": sample_id, "trait": "plasmid_replicon",
        "analysis_date": datetime.date.today().isoformat(), "schema": "plasmid-replicon-call-v0",
        "status": res["status"],
        "replicons": [r["replicon"] for r in res.get("replicons", [])],
        "replicon_detail": res.get("replicons", []),
        "caller": {"name": "dna_decode-plasmidfinder-blastn-v0", "method": res.get("method"),
                   "source": "PlasmidFinder curated allele DB", "caller_is_independent_baseline": False},
        "caveat": ("Faithful to PlasmidFinder's own method (blastn over its allele DB); NOT an independent "
                   "baseline. Replicon presence != a complete plasmid; absence is bounded by the DB's "
                   "enterobacteriales scope. NOT a clinical tool."),
        "provenance": {"db": str(args.db), "db_name": "plasmidfinder_enterobacteriales", "db_sha256_16": db_sha,
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage},
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: plasmid replicon typing")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
            if res.get("tool") == "db":
                print(f"  download the DB: curl -sSL {_DB_URL} -o {args.db}")
        else:
            reps = res["replicons"]
            print(f"REPLICONS: {len(reps)} found")
            for r in reps:
                print(f"  {r['replicon']:24} {r['percent_identity']}% id / "
                      f"{r['percent_coverage']}% cov  ({r['best_allele']})")
            if not reps:
                print("  (no plasmid replicons detected at threshold)")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
