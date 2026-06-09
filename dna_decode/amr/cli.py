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

from dna_decode.data.fungal_amr import (
    call_from_observed_substitutions,
    supported_fungal_drugs,
)
from dna_decode.data.mic_tiers import supported_drugs
from dna_decode.eval.amr_rules import AMRFINDER_IMAGE_PINNED, call_resistance

# Fungal target-site path (BLAST ERG11/FKS1 -> catalog), NOT AMRFinder (no AMRFinder-for-fungi). Routed by
# drug: fluconazole/voriconazole/caspofungin/micafungin -> fungal engine. G1-validated on C. auris
# (2026-06-08, wiki/fungal_ep7_g1_closeout): method transfers, sens 1.0 across clades; label-limited spec.
_DEFAULT_ERG11_REF = Path(__file__).resolve().parent.parent.parent / "data" / "fungal_ref" / "Cauris_ERG11_cds.fna"


def _parse_observed(observed: str) -> dict[str, set[str]]:
    """'ERG11:Y132F,ERG11:K143R,FKS1:S639F' -> {'ERG11': {'Y132F','K143R'}, 'FKS1': {'S639F'}}."""
    out: dict[str, set[str]] = {}
    for tok in (t.strip() for t in observed.split(",") if t.strip()):
        gene, _, sub = tok.partition(":")
        if not sub:
            raise ValueError(f"bad --observed token {tok!r}; expected GENE:SUBSTITUTION (e.g. ERG11:Y132F)")
        out.setdefault(gene.strip(), set()).add(sub.strip())
    return out


def _fungal_record(call, sample_id: str, drug: str, organism: str, provenance: dict) -> dict:
    """Map a FungalCall onto the uniform amr-mechanism-call-v1 record (same shape as the bacterial path)."""
    dets = [{"symbol": d.split(":", 1)[0], "subclass": d.split(":", 1)[1] if ":" in d else "",
             "class": "TARGET_SITE_MUTATION", "pct_identity": None} for d in call.determinants]
    return {
        "sample_id": sample_id, "drug": drug,
        "analysis_date": datetime.date.today().isoformat(), "schema": "amr-mechanism-call-v1",
        "prediction": call.prediction,
        "confidence": "high" if call.prediction == "R" else ("n/a" if call.prediction == "INDETERMINATE" else "screen"),
        "n_determinants": len(call.determinants), "determinants": dets,
        "resistance_threshold": 1,
        "undetectable_mechanisms": list(call.undetectable_mechanisms),
        "caller": {"name": "dna_decode-fungal-target-mutation-v0", "rule": call.rule,
                   "source": "hand-curated fungal catalog (no AMRFinder-for-fungi)",
                   "caller_is_independent_baseline": False},
        "caveat": call.caveat,
        "provenance": {**provenance, "organism": organism},
    }


def _fungal_main(args) -> int:
    """Fungal target-site decoder branch (routed when --drug is a fungal drug)."""
    # the --organism default 'Escherichia' is the bacterial AMRFinder default; on the fungal path, relabel
    # to the validated fungal organism unless the user explicitly set one.
    if args.organism == "Escherichia":
        args.organism = "Candida_auris"
    if args.observed is not None:
        obs = _parse_observed(args.observed)
        call = call_from_observed_substitutions(args.drug, obs)
        sample_id = args.sample_id or "observed"
        prov = {"mode": "observed-substitutions", "observed": args.observed}
    elif args.genome_fasta is not None:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            from scripts.fungal_erg11_caller import call_erg11  # repo-only; needs BLAST+
        except ImportError as e:
            print(f"ERROR: fungal genome mode needs scripts/fungal_erg11_caller + BLAST+ ({e}). "
                  "Use --observed with known substitutions for a wheel-only call.", file=sys.stderr)
            return 3
        call = call_erg11(str(args.genome_fasta), str(args.erg11_ref), args.drug)
        prov = {"mode": "blast-erg11", "erg11_ref": str(args.erg11_ref)}
    else:
        print("ERROR: fungal drug needs --genome-fasta OR --observed GENE:SUB[,...]", file=sys.stderr)
        return 2

    rec = _fungal_record(call, sample_id, args.drug, args.organism, prov)
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  drug: {args.drug}  organism: {args.organism}")
        print(f"CALL: {call.prediction}  [{rec['confidence']} | {len(call.determinants)} determinant(s)]")
        for d in call.determinants:
            print(f"  driven by: {d}")
        if not call.determinants:
            print("  driven by: (no catalogued target-site resistance mutation)")
        if call.undetectable_mechanisms:
            print(f"  blind spots (an S call can't rule out): {', '.join(call.undetectable_mechanisms)}")
        print(f"  {call.caveat}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if call.prediction != "INDETERMINATE" else 4


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
    ap.add_argument("--drug", required=True, choices=sorted(set(supported_drugs()) | set(supported_fungal_drugs())),
                    metavar="DRUG", help="bacterial (AMRFinder engine) or fungal (BLAST-ERG11 engine) drug")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--amrfinder-run", type=Path, help="[bacterial] existing AMRFinder run dir (main.tsv)")
    src.add_argument("--genome-fasta", type=Path, help="genome FASTA (bacterial: AMRFinder via Docker; "
                                                       "fungal: BLAST ERG11 via scripts/fungal_erg11_caller)")
    src.add_argument("--observed", default=None, help="[fungal] known substitutions GENE:SUB[,...] "
                                                      "(e.g. ERG11:Y132F) — pure, no BLAST")
    ap.add_argument("--erg11-ref", type=Path, default=_DEFAULT_ERG11_REF,
                    help="[fungal] in-frame ERG11 CDS reference FASTA (default: committed C. auris allele)")
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

    # Route by drug: fungal drugs -> BLAST-ERG11 target-site engine; bacterial -> AMRFinder engine.
    if args.drug in supported_fungal_drugs():
        if args.amrfinder_run is not None:
            print("ERROR: --amrfinder-run is bacterial-only; fungal drugs use --genome-fasta or --observed",
                  file=sys.stderr)
            return 2
        return _fungal_main(args)

    if args.observed is not None:
        print("ERROR: --observed is fungal-only; bacterial drugs use --amrfinder-run or --genome-fasta",
              file=sys.stderr)
        return 2

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

    # Forward --organism so a calibrated registry entry (opt-in) is used when the organism is known. The
    # default 'Escherichia' has no registry entry -> DRUG_RULE default (unchanged); an explicit
    # Campylobacter/Klebsiella/Salmonella resolves its independent-cohort-validated config, and an
    # EXPRESSION_FLOOR organism (Acinetobacter/Pseudomonas carbapenem) returns prediction 'ABSTAIN'.
    call = call_resistance(run_dir / "main.tsv", args.drug, args.resistance_threshold,
                           organism=args.organism)
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
        print(f"sample: {sample_id}  drug: {args.drug}  organism: {args.organism}")
        if call["prediction"] == "ABSTAIN":
            print("CALL: ABSTAIN  [gene-presence cannot decode this organism×drug]")
            print(f"  {call['caveat']}")
        else:
            nd = call["n_determinants"]
            print(f"CALL: {call['prediction']}  [{call['confidence']} | {nd} determinant(s)]")
            for x in call["determinants"]:
                print(f"  driven by: {x['symbol']}  ({x['subclass'] or x['class']}, {x['pct_identity']}% id)")
            if not call["determinants"]:
                print("  driven by: (no curated resistance determinants for this drug)")
            print(f"  {call['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return {"INDETERMINATE": 4, "ABSTAIN": 5}.get(call["prediction"], 0)


if __name__ == "__main__":
    raise SystemExit(main())
