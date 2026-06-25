"""Streptococcus pneumoniae capsular serotype typing -- in-package CLI (console `dna-pneumo-serotype`;
also `dna-decode pneumoserotype`).

Genome FASTA -> capsular serotype (cps-reference DB via real blastn) + provenance. Deterministic curated-DB
caller (sibling of dna-serotype / dna-ktype). Offline-safe (status "unavailable", exit 3). DB = a dir with
cps_references.fasta (per-serotype cps reference sequences; derives from PneumoCaT / SeroBA references).

    dna-pneumo-serotype assembly.fna --sample-id MY_STRAIN
    dna-pneumo-serotype X.fna --db-dir data/pneumoserotype_db
    # DB derives from PneumoCaT (github phe-bioinformatics/PneumoCaT) / SeroBA; see wiki/pneumo_serotype_report_card.md.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.pneumoserotype.runner import call_pneumo_serotype

DEFAULT_DB_DIR = "data/pneumoserotype_db"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-pneumo-serotype",
                                 description="Deterministic S. pneumoniae capsular serotyping (cps-reference via blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR,
                    help=f"dir with cps_references.fasta (default {DEFAULT_DB_DIR})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--identity", type=float, default=90.0, help="min %% identity (default 90)")
    ap.add_argument("--coverage", type=float, default=70.0, help="min %% coverage (default 70)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    cps = Path(args.db_dir) / "cps_references.fasta"
    db_sha = hashlib.sha256(cps.read_bytes()).hexdigest()[:16] if cps.exists() else None

    res = call_pneumo_serotype(args.fasta, args.db_dir, identity_threshold=args.identity,
                               coverage_threshold=args.coverage)
    rec = {
        "sample_id": sample_id, "trait": "pneumo_serotype", "organism": "Streptococcus pneumoniae",
        "analysis_date": datetime.date.today().isoformat(), "schema": "pneumo-serotype-call-v0",
        "status": res["status"], "serotype": res.get("serotype"), "best_reference": res.get("best_reference"),
        "percent_identity": res.get("percent_identity"), "percent_coverage": res.get("percent_coverage"),
        "caller": {"name": "dna_decode-cps-reference-blastn-v0", "method": res.get("method"),
                   "source": "PneumoCaT / SeroBA cps reference set",
                   "caller_is_independent_baseline": False},
        "caveat": ("Faithful to the cps-reference typing method (blastn best-match); NOT an independent "
                   "baseline. Single-best-reference v0 resolves SEROGROUP reliably; within-serogroup pairs "
                   "(6A/6B, 19A/19F) need allele-level logic the full tools add. Published in-silico-vs-"
                   "Quellung ceiling ~89% (GPS pipeline). S. pneumoniae only. NOT a clinical tool."),
        "provenance": {"db_dir": str(args.db_dir), "db_name": "pneumo_cps_reference", "cps_sha256_16": db_sha,
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage},
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")

    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: S. pneumoniae capsular serotype")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
        else:
            st = res.get("serotype")
            print(f"SEROTYPE: {st or '(none resolved)'}  "
                  f"(best ref {res.get('best_reference')}, {res.get('percent_identity')}% id / "
                  f"{res.get('percent_coverage')}% cov)")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
