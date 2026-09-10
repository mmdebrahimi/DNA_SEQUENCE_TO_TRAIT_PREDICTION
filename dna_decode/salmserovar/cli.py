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

from dna_decode.salmserovar.runner import (
    SEROVAR_COVERAGE_THRESHOLD,
    SEROVAR_IDENTITY_THRESHOLD,
    call_serovar,
)

DEFAULT_DB_DIR = "data/salmserovar_db"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-salmserovar",
                                 description="Deterministic Salmonella serovar typing (antigenic formula via blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR,
                    help=f"dir with salmonella_antigens.fasta + serovar_table.tsv (default {DEFAULT_DB_DIR})")
    ap.add_argument("--sample-id", default=None)
    # DERIVED, never restated. This CLI shipped `--coverage default=80.0` from the cell's first commit
    # and was NOT updated when SEROVAR_COVERAGE_THRESHOLD was lowered 80 -> 40 on 2026-09-04 against a
    # pre-registered bar. Because `main` passes `args.coverage` explicitly, the shipped entry point kept
    # overriding the validated value with the old one: every validation ran through `call_serovar`
    # (which uses the constant), so nothing exercised the seam. Defaults now come FROM the constants.
    ap.add_argument("--identity", type=float, default=SEROVAR_IDENTITY_THRESHOLD,
                    help=f"min %% identity (default {SEROVAR_IDENTITY_THRESHOLD:g})")
    ap.add_argument("--coverage", type=float, default=SEROVAR_COVERAGE_THRESHOLD,
                    help=f"min %% coverage (default {SEROVAR_COVERAGE_THRESHOLD:g})")
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
        # WHICH RULE produced the O antigen is part of the call, not decoration: a differential-marker
        # decision, a branch default and an ordinary best-allele pick are three different kinds of
        # evidence, and reporting only the antigen string collapses them.
        "o_antigen_rule": res.get("o_antigen_rule"),
        "o_antigen_best_hit": res.get("o_antigen_best_hit"),
        # Best hit PER AXIS -- for O this is often NOT the call (each row carries `is_call`).
        "antigen_detail": res.get("antigens", []),
        "caller": {"name": "dna_decode-salmserovar-blastn-v0", "method": res.get("method"),
                   "source": "SeqSero2 antigen DB + Kauffmann-White-Le Minor scheme",
                   "caller_is_independent_baseline": False},
        # This string PRINTS on every human-readable run, so a false claim here is a shipped falsehood,
        # not a stale doc. It said "reported only when the formula resolves uniquely" until 2026-09-09 --
        # untrue, because the table keeps an arbitrary winner on the 195 contested formulas
        # (ambiguous_policy="first") and there is an O+H1 phase-incomplete fallback on top.
        "caveat": ("Faithful to the SeqSero2 / Kauffmann-White method (blastn over the antigen allele DB + "
                   "formula lookup); NOT an independent baseline. A serovar is reported when the formula "
                   "RESOLVES -- which is not the same as being unique: 195 formulas are shared by >1 "
                   "Kauffmann-White serovar and the table records an arbitrary winner for those, and an "
                   "O+H1 fallback covers phase-incomplete genomes. Which of those produced a given call "
                   "is NOT yet disclosed. No serovar (formula-only, like O?/H?) means it did not resolve "
                   "at all. Salmonella enterica only. NOT a clinical tool."),
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
                # An evidence row that is NOT the call must say so. The O procedure routinely calls an
                # antigen that has no allele of its own (plain O9 is reached from a `wbaV` hit), so an
                # unmarked `O 9,46` line beside a `9:...` formula reads as a contradiction.
                label = "" if a.get("is_call", True) else "  <- best hit, NOT the call"
                print(f"  {a['axis']:3} {a['antigen']:8} {a['percent_identity']}% id / "
                      f"{a['percent_coverage']}% cov  ({a['best_allele']}){label}")
            if res.get("o_antigen_rule"):
                print(f"  O call: {res.get('o_antigen') or 'unresolved'}"
                      f"   [rule: {res['o_antigen_rule']}]")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
