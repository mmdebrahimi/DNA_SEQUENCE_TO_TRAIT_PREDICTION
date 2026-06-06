"""Deterministic, interpretable AMR decoder CLI (the Phase-2-decision shipping artifact).

Genome FASTA (or a cached AMRFinder run) -> R/S call per drug + the exact curated resistance
determinants that drove it + provenance. Mechanism-feature decoding, NOT embeddings (per
`plans/AMR_embedding_niche_decision_2026-06-05.md`: embeddings have no E. coli AMR niche; this is the
honest interpretable tool the decision points to). The sibling of `dna_decode/pathotype` (dna-pathotype).

    # from a cached AMRFinder run dir:
    python -m scripts.amr_predict --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_xxx.x
    # from a genome (runs AMRFinder via Docker; needs the DB at data/amrfinder_db):
    python -m scripts.amr_predict --drug ciprofloxacin --genome-fasta X.fna --sample-id X

NOT a clinical decision tool. The call is a transparent rule over AMRFinder's curated database; op-chars
(cipro N=147): acc 0.85 / sens 0.96 / spec 0.75 — sensitivity-favoring (see amr_rules caveat).
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.mic_tiers import supported_drugs
from dna_decode.eval.amr_rules import call_resistance


def _run_amrfinder_for_genome(fasta: Path, sample_id: str, out_root: Path, db: Path) -> Path:
    """Run AMRFinder on a single genome via the existing drug_mechanism_audit runner; return its run dir."""
    from scripts.drug_mechanism_audit import _run_amrfinder
    out_dir = out_root / (sample_id or fasta.stem)
    out_dir.mkdir(parents=True, exist_ok=True)
    # _run_amrfinder mounts the resolved DB (handles the latest-symlink/Docker-D: fixes) + writes main.tsv
    import scripts.drug_mechanism_audit as dma
    if db:
        dma.AMRFINDER_DB = str(db)
    _run_amrfinder(fasta, out_dir)
    return out_dir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="amr_predict",
                                 description="Deterministic AMR R/S decoder from AMRFinder curated determinants")
    ap.add_argument("--drug", required=True, choices=supported_drugs())
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--amrfinder-run", type=Path, help="existing AMRFinder run dir (contains main.tsv)")
    src.add_argument("--genome-fasta", type=Path, help="genome FASTA (runs AMRFinder via Docker)")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--amrfinder-db", type=Path, default=Path("data/amrfinder_db"),
                    help="AMRFinder DB root (Docker-readable; default data/amrfinder_db)")
    ap.add_argument("--out-root", type=Path, default=Path("data/amrfinder_runs"))
    ap.add_argument("--resistance-threshold", type=int, default=2,
                    help="min #curated determinants for an R call (default 2, cipro/QRDR-validated; "
                         "use 1 for acquired-gene-dominant drugs e.g. cef beta-lactamases)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if args.amrfinder_run:
        run_dir = args.amrfinder_run
        sample_id = args.sample_id or run_dir.name
    else:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            run_dir = _run_amrfinder_for_genome(args.genome_fasta, sample_id, args.out_root, args.amrfinder_db)
        except Exception as e:  # Docker/AMRFinder failure → honest error, no fabricated call
            print(f"ERROR: AMRFinder run failed ({type(e).__name__}: {e}). "
                  f"Need Docker + a Docker-readable DB at {args.amrfinder_db}.", file=sys.stderr)
            return 3

    call = call_resistance(run_dir / "main.tsv", args.drug, args.resistance_threshold)
    rec = {
        "sample_id": sample_id,
        "drug": args.drug,
        "analysis_date": datetime.date.today().isoformat(),
        "schema": "amr-mechanism-call-v0",
        "prediction": call["prediction"],
        "confidence": call["confidence"],
        "n_determinants": call["n_determinants"],
        "determinants": call["determinants"],
        "caller": {"name": "dna_decode-amr-rules-v0", "rule": call["rule"],
                   "source": "AMRFinderPlus curated main.tsv", "caller_is_independent_baseline": False},
        "caveat": call["caveat"],
        "provenance": {"amrfinder_run": str(run_dir)},
    }
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        d = call["determinants"]
        print(f"sample: {sample_id}  drug: {args.drug}")
        print(f"CALL: {call['prediction']}  [{call['confidence']} | {call['n_determinants']} determinant(s)]")
        for x in d:
            print(f"  driven by: {x['symbol']}  ({x['subclass'] or x['class']}, {x['pct_identity']}% id)")
        if not d:
            print("  driven by: (no curated resistance determinants for this drug)")
        print(f"  {call['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if call["prediction"] != "INDETERMINATE" else 4


if __name__ == "__main__":
    raise SystemExit(main())
