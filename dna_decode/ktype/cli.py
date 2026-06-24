"""Klebsiella K-antigen (capsule) typing -- in-package CLI (console `dna-ktype`; also `dna-decode ktype`).

Genome FASTA -> predicted Klebsiella K-locus (KL) type via the wzi allele scheme (real blastn) + provenance.
Deterministic curated-DB caller (sibling of dna-serotype / dna-amr). Offline-safe (status 'unavailable',
exit 3). DB = a dir with wzi.fasta + wzi.txt (the BIGSdb Pasteur wzi scheme as bundled by Kleborate),
downloaded on demand.

    dna-ktype assembly.fna --sample-id MY_STRAIN
    dna-ktype X.fna --db-dir data/ktype_db
    # fetch the DB once (Kleborate-bundled BIGSdb wzi scheme):
    #   B=https://raw.githubusercontent.com/klebgenomics/Kleborate/main/kleborate/modules/klebsiella_pneumo_complex__wzi/data
    #   curl -L $B/wzi.fasta -o data/ktype_db/wzi.fasta ; curl -L $B/wzi.txt -o data/ktype_db/wzi.txt
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.ktype.runner import call_ktype

DEFAULT_DB_DIR = "data/ktype_db"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-ktype",
                                 description="Deterministic Klebsiella K-antigen (capsule) typing via wzi (blastn)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR, help=f"dir with wzi.fasta + wzi.txt (default {DEFAULT_DB_DIR})")
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
    wzi_fasta = Path(args.db_dir) / "wzi.fasta"
    db_sha = hashlib.sha256(wzi_fasta.read_bytes()).hexdigest()[:16] if wzi_fasta.exists() else None

    res = call_ktype(args.fasta, args.db_dir, identity_threshold=args.identity, coverage_threshold=args.coverage)
    rec = {
        "sample_id": sample_id, "trait": "ktype", "organism": "Klebsiella",
        "analysis_date": datetime.date.today().isoformat(), "schema": "ktype-call-v0",
        "status": res["status"], "predicted_k": res.get("predicted_k"), "kl_type": res.get("kl_type"),
        "wzi_allele": res.get("wzi_allele"),
        "percent_identity": res.get("percent_identity"), "percent_coverage": res.get("percent_coverage"),
        "caller": {"name": "dna_decode-wzi-blastn-v0", "method": res.get("method"),
                   "source": "BIGSdb Pasteur wzi scheme (Kleborate-bundled)",
                   "caller_is_independent_baseline": False},
        "caveat": ("Faithful to the Kleborate/BIGSdb wzi typing method (blastn over the wzi allele DB); NOT "
                   "an independent baseline. wzi->K-type is ~94% predictive and NOT one-to-one (isolates "
                   "with distinct K-types can share a wzi allele; Brisse 2013 JCM). Full-locus typing "
                   "(Kaptive) is more accurate; this is the single-gene v0. Klebsiella only. NOT a clinical tool."),
        "provenance": {"db_dir": str(args.db_dir), "db_name": "wzi_bigsdb_pasteur", "wzi_fasta_sha256_16": db_sha,
                       "identity_threshold": args.identity, "coverage_threshold": args.coverage},
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")

    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: K-antigen (capsule)  organism: Klebsiella")
        if res["status"] != "ok":
            print(f"CALL: [{res['status']}]  ({res.get('reason')})")
        else:
            k = res.get("predicted_k")
            print(f"CALL: {k or 'K?'}  (wzi allele {res.get('wzi_allele')}, "
                  f"{res.get('percent_identity')}% id / {res.get('percent_coverage')}% cov)")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
