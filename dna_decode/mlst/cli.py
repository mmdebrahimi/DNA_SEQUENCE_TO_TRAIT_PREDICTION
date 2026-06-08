"""MLST sequence typing — CLI (console `dna-mlst`; also `dna-decode mlst`).

Genome FASTA -> 7-gene Sequence Type (ST) via exact per-locus allele matching + the PubMLST profiles table.
v0 scope: E. coli Achtman (adk/fumC/gyrB/icd/mdh/purA/recA). Offline-safe (status 'unavailable', exit 3).
`--fetch-db` pulls the scheme (7 loci FASTAs + profiles.tsv) from PubMLST into --db-dir.

    dna-mlst --fetch-db                       # one-time: download the E. coli Achtman scheme
    dna-mlst assembly.fna --sample-id MY_STRAIN
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
import urllib.request
from pathlib import Path

from dna_decode.mlst.runner import call_mlst

DEFAULT_DB_DIR = "data/mlst_db/ecoli_achtman"
LOCI = ["adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA"]
_PUB = "https://rest.pubmlst.org/db/pubmlst_ecoli_achtman_seqdef"
_SCHEME = "4"


def _fetch_db(db_dir: Path) -> int:
    db_dir.mkdir(parents=True, exist_ok=True)
    try:
        for loc in LOCI:
            url = f"{_PUB}/loci/{loc}/alleles_fasta"
            with urllib.request.urlopen(url, timeout=120) as r:
                (db_dir / f"{loc}.fasta").write_bytes(r.read())
            print(f"  fetched {loc}.fasta")
        with urllib.request.urlopen(f"{_PUB}/schemes/{_SCHEME}/profiles_csv", timeout=120) as r:
            (db_dir / "profiles.tsv").write_bytes(r.read())
        print(f"  fetched profiles.tsv -> {db_dir}")
        return 0
    except Exception as e:
        print(f"ERROR: PubMLST fetch failed ({type(e).__name__}: {e})", file=sys.stderr)
        return 3


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-mlst", description="Deterministic MLST sequence typing (PubMLST)")
    ap.add_argument("fasta", type=Path, nargs="?", help="genome assembly FASTA (omit with --fetch-db)")
    ap.add_argument("--db-dir", default=DEFAULT_DB_DIR, help=f"scheme dir (loci .fasta + profiles.tsv; default {DEFAULT_DB_DIR})")
    ap.add_argument("--loci", default=",".join(LOCI), help="comma locus list")
    ap.add_argument("--fetch-db", action="store_true", help="download the E. coli Achtman scheme from PubMLST")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    db = Path(args.db_dir)
    if args.fetch_db:
        return _fetch_db(db)
    if not args.fasta or not args.fasta.exists():
        print("ERROR: genome FASTA required (or run --fetch-db first)", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem
    loci = [x.strip() for x in args.loci.split(",") if x.strip()]
    locus_fastas = {loc: db / f"{loc}.fasta" for loc in loci}
    profiles = db / "profiles.tsv"
    if not profiles.exists() or not any(p.exists() for p in locus_fastas.values()):
        rec = {"sample_id": sample_id, "trait": "mlst", "status": "unavailable", "schema": "mlst-call-v0",
               "reason": f"MLST scheme not found at {args.db_dir} — run: dna-mlst --fetch-db --db-dir {args.db_dir}"}
        if args.out:
            Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
        print(json.dumps(rec, indent=2) if args.json_only else f"STATUS: unavailable - {rec['reason']}")
        return 3

    res = call_mlst(args.fasta, {l: p for l, p in locus_fastas.items() if p.exists()}, profiles)
    rec = {
        "sample_id": sample_id, "trait": "mlst",
        "analysis_date": datetime.date.today().isoformat(), "schema": "mlst-call-v0",
        "status": res["status"],
        "sequence_type": res.get("st"), "clonal_complex": res.get("clonal_complex"),
        "profile": res.get("profile", {}), "complete": res.get("complete"), "novel": res.get("novel"),
        "scheme": "ecoli_achtman", "scheme_loci": res.get("scheme_loci", []),
        "caller": {"name": "dna_decode-mlst-pubmlst-v0", "method": res.get("method"),
                   "source": "PubMLST E. coli Achtman scheme", "caller_is_independent_baseline": False},
        "caveat": ("Exact per-locus allele match (blastn 100/100) + PubMLST profile lookup. A locus with no "
                   "exact allele -> null (novel/partial allele or assembly gap) -> ST not called. 'novel' = "
                   "full profile not yet in the ST table. v0: E. coli Achtman only. NOT a clinical tool."),
    }
    if res["status"] != "ok":
        rec["reason"] = res.get("reason")
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  trait: MLST (E. coli Achtman)")
        if res["status"] != "ok":
            print(f"STATUS: {res['status']} - {res.get('reason')}")
        else:
            st = res.get("st")
            label = f"ST{st}" if st else ("novel ST (full profile, not in table)" if res.get("novel")
                                          else "ST not called (incomplete profile)")
            print(f"  {label}" + (f"  [{res['clonal_complex']}]" if res.get("clonal_complex") else ""))
            print("  profile: " + "  ".join(f"{loc}={res['profile'].get(loc)}" for loc in res["scheme_loci"]))
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
