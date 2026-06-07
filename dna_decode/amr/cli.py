"""Deterministic AMR mechanism decoder — in-package CLI (console entry `dna-amr`).

Genome FASTA (or a cached AMRFinder run) -> R/S call per drug + the curated resistance determinants that
drove it + provenance. Mechanism-feature decoding, NOT embeddings (per
`plans/AMR_embedding_niche_decision_2026-06-05.md`). Sibling of `dna_decode.pathotype.cli` (dna-pathotype).

In-package so it ships in the wheel. Cached-run mode is pure (reads main.tsv via amr_rules — no Docker).
Genome mode lazily imports the AMRFinder Docker runner from `scripts/` (repo-only; needs Docker + a
Docker-readable DB) and errors cleanly if unavailable — so the console entry installs + imports without
the scripts/ dir.

    dna-amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_xxx.x
    dna-amr --drug ciprofloxacin --genome-fasta X.fna --sample-id X      # needs Docker + data/amrfinder_db

NOT a clinical decision tool. cipro N=147 op-chars (threshold=2): acc 0.939 / sens 0.931 / spec 0.947.
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from dna_decode.data.mic_tiers import supported_drugs
from dna_decode.eval.amr_rules import AMRFINDER_IMAGE_PINNED, call_resistance


def _run_amrfinder_for_genome(fasta: Path, sample_id: str, out_root: Path, db: Path,
                              organism: str = "Escherichia") -> Path:
    """Genome mode: lazily import the repo's AMRFinder Docker runner (not in the wheel).

    `organism` selects AMRFinder's `-O` (organism-specific point-mutation detection — gyrA/parC QRDR calls
    are organism-specific, so a Klebsiella genome MUST use 'Klebsiella_pneumoniae' or its QRDR is missed)."""
    try:
        import scripts.drug_mechanism_audit as dma  # repo-only; needs Docker + DB
        from scripts.drug_mechanism_audit import _run_amrfinder
    except ImportError as e:
        raise RuntimeError(
            "genome mode needs the repo's AMRFinder runner (scripts/drug_mechanism_audit) + Docker + a "
            "Docker-readable DB at --amrfinder-db; not available in a wheel install. Use --amrfinder-run "
            f"with a precomputed run instead. ({e})"
        ) from e
    out_dir = out_root / (sample_id or fasta.stem)
    out_dir.mkdir(parents=True, exist_ok=True)
    if db:
        dma.AMRFINDER_DB = str(db)
    _run_amrfinder(fasta, out_dir, organism=organism)
    return out_dir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-amr",
                                 description="Deterministic AMR R/S decoder from AMRFinder curated determinants")
    ap.add_argument("--drug", required=True, choices=supported_drugs())
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--amrfinder-run", type=Path, help="existing AMRFinder run dir (contains main.tsv)")
    src.add_argument("--genome-fasta", type=Path, help="genome FASTA (runs AMRFinder via Docker)")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--organism", default="Escherichia",
                    help="AMRFinder -O organism for genome mode (organism-specific QRDR point-mutation "
                         "detection). E.g. Escherichia (default), Klebsiella_pneumoniae. Validated "
                         "cross-organism for E. coli + K. pneumoniae.")
    ap.add_argument("--amrfinder-db", type=Path, default=Path("data/amrfinder_db"),
                    help="AMRFinder DB root (Docker-readable; default data/amrfinder_db)")
    ap.add_argument("--out-root", type=Path, default=Path("data/amrfinder_runs"))
    ap.add_argument("--resistance-threshold", type=int, default=None,
                    help="min #curated determinants for an R call. Default: per-drug validated config "
                         "(cipro=2 QRDR; cef=1 + extended-spectrum refinement; tet/gent=1). Pass an int to override.")
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
            run_dir = _run_amrfinder_for_genome(args.genome_fasta, sample_id, args.out_root,
                                                args.amrfinder_db, organism=args.organism)
        except Exception as e:
            print(f"ERROR: AMRFinder run failed ({type(e).__name__}: {e}).", file=sys.stderr)
            return 3

    call = call_resistance(run_dir / "main.tsv", args.drug, args.resistance_threshold)
    rec = {
        "sample_id": sample_id, "drug": args.drug,
        "analysis_date": datetime.date.today().isoformat(), "schema": "amr-mechanism-call-v1",
        "prediction": call["prediction"], "confidence": call["confidence"],
        "n_determinants": call["n_determinants"], "determinants": call["determinants"],
        "resistance_threshold": call.get("resistance_threshold"),
        "undetectable_mechanisms": call.get("undetectable_mechanisms", []),
        "caller": {"name": "dna_decode-amr-rules-v1", "rule": call["rule"],
                   "source": "AMRFinderPlus curated main.tsv", "caller_is_independent_baseline": False},
        "caveat": call["caveat"],
        "provenance": {"amrfinder_run": str(run_dir), "amrfinder_image": AMRFINDER_IMAGE_PINNED,
                       "amrfinder_organism": args.organism},
    }
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  drug: {args.drug}")
        print(f"CALL: {call['prediction']}  [{call['confidence']} | {call['n_determinants']} determinant(s)]")
        for x in call["determinants"]:
            print(f"  driven by: {x['symbol']}  ({x['subclass'] or x['class']}, {x['pct_identity']}% id)")
        if not call["determinants"]:
            print("  driven by: (no curated resistance determinants for this drug)")
        print(f"  {call['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if call["prediction"] != "INDETERMINATE" else 4


if __name__ == "__main__":
    raise SystemExit(main())
