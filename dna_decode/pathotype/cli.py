"""v0 pathotype compatibility-resolver CLI.

    python -m dna_decode.pathotype <assembly.fasta> [--gff3 G] [--db DB] [--out J] [--sample-id ID]

FASTA-in (+ optional GFF3, accepted for forward-compat; v0 detection is assembly-only)
-> audit-grade provenance JSON (ledger v0 Output Contract schema) + a human-readable
summary. Emits a `*_COMPATIBLE` call with virulence-cluster provenance and abstention.
NOT a clinical predictor (see markers.py compatibility-resolver framing).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from dna_decode.pathotype.detect import (
    CONFIDENT_COV, PARTIAL_COV, assembly_qc, build_vf_index, detect, parse_fasta,
)
from dna_decode.pathotype.markers import RULES_VERSION
from dna_decode.pathotype.resolve import resolve_call

DEFAULT_DB = "data/virulencefinder_db/virulence_ecoli.fsa"


def analyze(fasta: str, db: str, sample_id: str | None, gff3: str | None = None,
            analysis_date: str | None = None) -> dict:
    contigs = parse_fasta(fasta)
    qc = assembly_qc(contigs)
    vf_index = build_vf_index(db)
    det = detect(contigs, vf_index)
    call = resolve_call(det["cluster_profile"],
                        partial_clusters=det["partial_clusters"],
                        qc_pass=(qc["qc_verdict"] != "FAIL"))
    db_sha = hashlib.sha256(Path(db).read_bytes()).hexdigest()[:16]
    return {
        "sample_id": sample_id or Path(fasta).stem,
        "analysis_date": analysis_date or datetime.date.today().isoformat(),
        "schema": "pathotype-provenance-v0",
        "assembly": {"path": str(fasta), **qc, "gff3": gff3},
        "caller": {
            "name": "dna_decode-kmer-seed-v0",
            "software_version": RULES_VERSION,
            "db_name": "virulencefinder_ecoli",
            "db_version": "bitbucket-HEAD",
            "db_commit_or_checksum": f"sha256:{db_sha}",
            "parameters": {"k": 15, "confident_coverage": CONFIDENT_COV,
                           "partial_coverage": PARTIAL_COV},
            "method": "kmer_seed_coverage_v0",
            "caller_is_independent_baseline": False,
        },
        "marker_hits": det["marker_hits"],
        "cluster_profile": {c: v for c, v in det["cluster_profile"].items() if v},
        "derived_call": call,
    }


def _summary(rec: dict) -> str:
    c = rec["derived_call"]
    clusters = ", ".join(sorted(rec["cluster_profile"])) or "(none)"
    lines = [
        f"sample: {rec['sample_id']}",
        f"assembly QC: {rec['assembly']['qc_verdict']} "
        f"({rec['assembly']['n_contigs']} contigs, {rec['assembly']['total_bp']:,} bp, N50 {rec['assembly']['n50']:,})",
        f"CALL: {c['primary']}  [{c['confidence_tier']} | {c['external_validity']}]",
    ]
    if c["secondary"]:
        lines.append(f"  secondary: {', '.join(c['secondary'])}")
    lines += [f"  rule: {c['rule_id']} ({c['rule_version']})",
              f"  reason: {c['reason']}",
              f"  clusters present: {clusters}",
              f"  driven by: " + "; ".join(
                  f"{h['cluster']}/{h['gene']} ({h['percent_coverage']}% cov, {h['hit_status']}, "
                  f"{h['contig']}:{h['start']})" for h in rec["marker_hits"] if h["hit_status"] == "CONFIDENT") or "  driven by: (no confident hits)"]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna_decode.pathotype",
                                 description="v0 E. coli pathotype compatibility resolver (FASTA -> provenance JSON)")
    ap.add_argument("fasta", help="genome assembly FASTA")
    ap.add_argument("--gff3", default=None, help="optional GFF3 (accepted; v0 is assembly-only)")
    ap.add_argument("--db", default=DEFAULT_DB, help=f"VirulenceFinder E. coli .fsa (default {DEFAULT_DB})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--out", default=None, help="write provenance JSON here (default: alongside summary on stdout)")
    ap.add_argument("--json-only", action="store_true", help="print only the JSON")
    args = ap.parse_args(argv)

    if not Path(args.db).exists():
        print(f"ERROR: marker DB not found: {args.db}\n"
              f"  fetch: curl -L -o {args.db} "
              f"https://bitbucket.org/genomicepidemiology/virulencefinder_db/raw/HEAD/virulence_ecoli.fsa",
              file=sys.stderr)
        return 2
    rec = analyze(args.fasta, args.db, args.sample_id, args.gff3)
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(_summary(rec))
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
        else:
            print("\n" + json.dumps(rec, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
