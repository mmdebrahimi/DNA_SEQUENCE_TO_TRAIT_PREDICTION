"""Chromosomal AMR point-mutation decoder — CLI (console `dna-pointfinder`; also `dna-decode pointfinder`).

Genome FASTA -> catalogued chromosomal resistance point mutations (PointFinder ref genes via blastn +
codon-position lookup) + the resistances they confer. An INDEPENDENT point-mutation caller (complements
dna-amr's AMRFinder POINT + dna-resfinder's acquired genes). v0 scope: E. coli fluoroquinolone QRDR
(gyrA/parC/gyrB/parE). Offline-safe (status 'unavailable', exit 3). DB downloaded on demand.

    dna-pointfinder assembly.fna --sample-id MY_STRAIN
    dna-pointfinder X.fna --db-dir data/pointfinder_db/escherichia_coli --genes gyrA,parC,gyrB,parE
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from dna_decode.pointfinder.runner import call_point_mutations, parse_overview

DEFAULT_DB_DIR = "data/pointfinder_db/escherichia_coli"
DEFAULT_GENES = "gyrA,parC,gyrB,parE"
_DB_BASE = "https://bitbucket.org/genomicepidemiology/pointfinder_db/raw/HEAD/escherichia_coli"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-pointfinder",
                                 description="Deterministic chromosomal AMR point-mutation caller (PointFinder)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR, help=f"PointFinder species dir (default {DEFAULT_DB_DIR})")
    ap.add_argument("--genes", default=DEFAULT_GENES, help=f"comma gene list (default {DEFAULT_GENES})")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    d = Path(args.db_dir)
    overview_path = d / "resistens-overview.txt"
    gene_refs = {g: d / f"{g}.fsa" for g in (x.strip() for x in args.genes.split(",")) if g}
    have = overview_path.exists() and any(p.exists() for p in gene_refs.values())
    if not have:
        rec = {"sample_id": sample_id, "trait": "chromosomal_point_mutations", "status": "unavailable",
               "schema": "pointfinder-call-v0",
               "reason": f"PointFinder DB not found at {args.db_dir} (need resistens-overview.txt + <gene>.fsa)"}
        if args.out:
            Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
        if args.json_only:
            print(json.dumps(rec, indent=2))
        else:
            print(f"STATUS: unavailable - {rec['reason']}")
            print(f"  download: curl -sSL {_DB_BASE}/resistens-overview.txt -o {overview_path}  "
                  f"(+ {args.genes} .fsa from {_DB_BASE}/)")
        return 3

    overview = parse_overview(overview_path)
    res = call_point_mutations(args.fasta, {g: p for g, p in gene_refs.items() if p.exists()}, overview)
    rec = {
        "sample_id": sample_id, "trait": "chromosomal_point_mutations",
        "analysis_date": datetime.date.today().isoformat(), "schema": "pointfinder-call-v0",
        "status": res["status"],
        "mutations": [m["mutation"] for m in res.get("mutations", [])],
        "resistances": res.get("resistances", []),
        "mutation_detail": res.get("mutations", []),
        "genes_aligned": res.get("genes_aligned", []),
        "caller": {"name": "dna_decode-pointfinder-blastn-v0", "method": res.get("method"),
                   "source": "PointFinder curated resistens-overview + ref genes",
                   "caller_is_independent_baseline": True},
        "caveat": ("INDEPENDENT chromosomal point-mutation caller (PointFinder DB) — complements dna-amr "
                   "(AMRFinder POINT) + dna-resfinder (acquired only). v0: E. coli QRDR (gyrA/parC/gyrB/parE); "
                   "epistasis (Required_mut) recorded but NOT enforced. NOT a clinical tool."),
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: chromosomal AMR point mutations")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
        else:
            print(f"MUTATIONS: {len(rec['mutations'])}  | resistances: {', '.join(rec['resistances']) or '(none)'}")
            for m in res.get("mutations", []):
                print(f"  {m['gene']} {m['mutation']:8} -> {', '.join(m['resistances'])}")
            if not rec["mutations"]:
                print(f"  (no catalogued point mutations; genes aligned: {', '.join(rec['genes_aligned']) or 'none'})")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
