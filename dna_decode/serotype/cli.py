"""E. coli serotype typing — in-package CLI (console `dna-serotype`; also `dna-decode serotype`).

Genome FASTA -> O:H serotype (SerotypeFinder allele DB via real blastn) + provenance. Deterministic
curated-DB caller (sibling of dna-amr / dna-pathotype / dna-plasmid). Offline-safe (status "unavailable",
exit 3). The DB is two files (O_type.fsa + H_type.fsa) concatenated, downloaded on demand.

    dna-serotype assembly.fna --sample-id MY_STRAIN
    dna-serotype X.fna --db data/serotypefinder_db/serotypefinder.fsa
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.serotype.runner import call_serotype

DEFAULT_DB = "data/serotypefinder_db/serotypefinder.fsa"
_DB_BASE = "https://bitbucket.org/genomicepidemiology/serotypefinder_db/raw/HEAD"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-serotype",
                                 description="Deterministic E. coli O:H serotyping (SerotypeFinder via blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db", default=DEFAULT_DB, help=f"SerotypeFinder allele .fsa (O+H; default {DEFAULT_DB})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--identity", type=float, default=85.0, help="min %% identity (SerotypeFinder default 85)")
    ap.add_argument("--coverage", type=float, default=60.0, help="min %% coverage (default 60)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    db_sha = hashlib.sha256(Path(args.db).read_bytes()).hexdigest()[:16] if Path(args.db).exists() else None

    res = call_serotype(args.fasta, args.db, identity_threshold=args.identity, coverage_threshold=args.coverage)
    rec = {
        "sample_id": sample_id, "trait": "serotype",
        "analysis_date": datetime.date.today().isoformat(), "schema": "serotype-call-v0",
        "status": res["status"], "serotype": res.get("serotype"),
        "o_antigen": res.get("o_antigen"), "h_antigen": res.get("h_antigen"),
        "antigen_detail": res.get("antigens", []),
        "caller": {"name": "dna_decode-serotypefinder-blastn-v0", "method": res.get("method"),
                   "source": "SerotypeFinder curated allele DB", "caller_is_independent_baseline": False},
        "caveat": ("Faithful to SerotypeFinder's method (blastn over its O/H allele DB); NOT an independent "
                   "baseline. O?/H? = that antigen not resolved at threshold (novel/partial/absent locus). "
                   "E. coli only. NOT a clinical tool."),
        "provenance": {"db": str(args.db), "db_name": "serotypefinder", "db_sha256_16": db_sha,
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage},
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: E. coli serotype")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
            if res.get("tool") == "db":
                print(f"  build the DB: curl -sSL {_DB_BASE}/O_type.fsa {_DB_BASE}/H_type.fsa > {args.db}")
        else:
            print(f"SEROTYPE: {res['serotype'] or '(none resolved)'}")
            for a in res.get("antigens", []):
                print(f"  {a['antigen']:5} via {a['gene']:5} {a['percent_identity']}% id / "
                      f"{a['percent_coverage']}% cov  ({a['best_allele']})")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
