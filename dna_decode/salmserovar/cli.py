"""Salmonella enterica serovar typing -- in-package CLI (console `dna-salmserovar`; also `dna-decode salmserovar`).

Genome FASTA -> Kauffmann-White antigenic formula (O:H1:H2) + serovar (SeqSero2-style antigen allele DB via
real blastn) + provenance. Deterministic curated-DB caller (sibling of dna-serotype / dna-ktype).
Offline-safe (status "unavailable", exit 3). DB = a dir with salmonella_antigens.fasta + serovar_table.tsv.

    dna-salmserovar assembly.fna --sample-id MY_STRAIN
    dna-salmserovar X.fna --db-dir data/salmserovar_db
    # DB derives from the SeqSero2 database (github denglab/SeqSero2); see wiki/salm_serovar_report_card.md.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.salmserovar.runner import call_serovar

DEFAULT_DB_DIR = "data/salmserovar_db"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-salmserovar",
                                 description="Deterministic Salmonella serovar typing (antigenic formula via blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR,
                    help=f"dir with salmonella_antigens.fasta + serovar_table.tsv (default {DEFAULT_DB_DIR})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--identity", type=float, default=90.0, help="min %% identity (default 90)")
    ap.add_argument("--coverage", type=float, default=80.0, help="min %% coverage (default 80)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    antigens = Path(args.db_dir) / "salmonella_antigens.fasta"
    db_sha = hashlib.sha256(antigens.read_bytes()).hexdigest()[:16] if antigens.exists() else None

    res = call_serovar(args.fasta, args.db_dir, identity_threshold=args.identity, coverage_threshold=args.coverage)
    rec = {
        "sample_id": sample_id, "trait": "serovar", "organism": "Salmonella enterica",
        "analysis_date": datetime.date.today().isoformat(), "schema": "serovar-call-v0",
        "status": res["status"], "serovar": res.get("serovar"),
        "antigenic_formula": res.get("antigenic_formula"),
        "o_antigen": res.get("o_antigen"), "h1_antigen": res.get("h1_antigen"), "h2_antigen": res.get("h2_antigen"),
        "antigen_detail": res.get("antigens", []),
        "caller": {"name": "dna_decode-salmserovar-blastn-v0", "method": res.get("method"),
                   "source": "SeqSero2 antigen DB + Kauffmann-White-Le Minor scheme",
                   "caller_is_independent_baseline": False},
        "caveat": ("Faithful to the SeqSero2 / Kauffmann-White method (blastn over the antigen allele DB + "
                   "formula lookup); NOT an independent baseline. Serovar reported only when the O:H1:H2 "
                   "formula resolves uniquely (else formula-only, like O?/H?). Salmonella enterica only. "
                   "NOT a clinical tool."),
        "provenance": {"db_dir": str(args.db_dir), "db_name": "seqsero2_antigen", "antigens_sha256_16": db_sha,
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage},
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")

    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: Salmonella serovar")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
        else:
            print(f"SEROVAR: {res.get('serovar') or '(formula unresolved)'}  "
                  f"[formula {res.get('antigenic_formula')}]")
            for a in res.get("antigens", []):
                print(f"  {a['axis']:3} {a['antigen']:8} {a['percent_identity']}% id / "
                      f"{a['percent_coverage']}% cov  ({a['best_allele']})")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
