"""`dna-decode forward` — the forward variant-effect decoder as a first-class CLI command.

The forward direction of the decoder: given a genetic EDIT (a point mutation on a protein), predict the change
in MOLECULAR phenotype (Regime B — enzyme fitness/stability), validated per-variant against measured Deep
Mutational Scanning (the one place the project's label wall does not bind). This is the "edit -> effect"
complement to the deterministic determinant->R/S AMR decoder (`dna-decode amr`), which this does NOT touch.

    dna-decode forward --mutation M69L --protein-seq MSIQHFRVALIPFFAAFCLPVFA...   # BLOSUM62, instant
    dna-decode forward --mutation S83L --protein-fasta gyrA.faa --protein gyrA
    dna-decode forward --mutation G12D --protein-seq <seq> --json

v0 CLI is BLOSUM62-only — deterministic, wheel-only, no network, no GPU (the same offline-safe posture as the
blastn decoders). The learned methods (ESM2-650M / AlphaMissense / ESM-IF) beat BLOSUM everywhere but need a
precomputed score table (they run the model ONCE per protein); they stay in the Python API
(`dna_decode.forward.predict_effect(..., method="esm2", esm_table=...)`; see `dna_decode/forward/README.md`).

HONEST SCOPE: this predicts MOLECULAR fitness rank (the DMS-validated quantity), NOT clinical antibiotic
resistance — the antagonistic-selection direction fails for likelihood/severity scorers; use `dna-decode amr`
(the frozen determinant catalogue) for R/S. NOT a clinical tool.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _read_fasta_seq(path: str) -> str:
    seq = "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))
    return "".join(c for c in seq.upper() if c.isalpha())


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode forward",
        description="Forward variant-effect decode: a point mutation -> predicted molecular-phenotype change "
                    "(Regime B, DMS-validated). BLOSUM62 (deterministic) in v0; learned methods via the API.",
        epilog="scope: molecular fitness rank, NOT clinical resistance (use `dna-decode amr` for R/S).",
    )
    ap.add_argument("--mutation", required=True, help="point mutation, e.g. M69L (WT-pos-ALT, 1-based)")
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--protein-seq", help="protein amino-acid sequence (verifies the WT residue at the position)")
    src.add_argument("--protein-fasta", help="path to a protein FASTA whose sequence to use")
    ap.add_argument("--protein", default="protein", help="protein name/label (default: protein)")
    ap.add_argument("--method", default="blosum62", choices=["blosum62"],
                    help="v0 CLI = blosum62 only (deterministic, offline). Learned methods (esm2/alphamissense/"
                         "esm_if) beat it but need a precomputed table -> use the Python API.")
    ap.add_argument("--regime", default="B_molecular", choices=["B_molecular", "C_organismal"],
                    help="B_molecular = enzyme fitness/stability (default); C_organismal = abstain (closed negative)")
    ap.add_argument("--json", action="store_true", help="emit the ForwardPrediction as JSON")
    args = ap.parse_args(argv)

    from dna_decode.forward import predict_effect  # lazy: keep the top-level CLI import light

    seq = ""
    if args.protein_fasta:
        try:
            seq = _read_fasta_seq(args.protein_fasta)
        except OSError as e:
            print(f"error: cannot read --protein-fasta {args.protein_fasta}: {e}", file=sys.stderr)
            return 2
    elif args.protein_seq:
        seq = "".join(c for c in args.protein_seq.upper() if c.isalpha())

    try:
        pred = predict_effect(seq, args.mutation, protein=args.protein, method=args.method, regime=args.regime)
    except ValueError as e:
        # WT/coordinate mismatch or malformed mutation — fail loudly, never a silent wrong call.
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(pred.as_dict(), indent=2))
        return 0

    d = pred.as_dict()
    score = d["raw_score"]
    print("forward variant-effect decode")
    print(f"  protein: {d['protein']}   mutation: {d['mutation']} ({d['wt']}->{d['alt']} at {d['pos']})")
    print(f"  method: {d['method']}   regime: {d['regime']}")
    print(f"  predicted_effect: {d['predicted_effect']}   "
          f"raw_score: {'nan' if score != score else f'{score:.3f}'}   confidence: {d['confidence']}")
    print(f"  phenotype axis: {d['phenotype_axis']}")
    if d["notes"]:
        for n in d["notes"]:
            print(f"  note: {n}")
    if not seq:
        print("  note: no protein sequence supplied — WT residue NOT verified against a reference")
    print("  [BLOSUM62 = deterministic substitution severity; learned methods (ESM2/AlphaMissense/ESM-IF) "
          "beat it but are API-only — see dna_decode/forward/README.md]")
    print("  [scope: molecular fitness rank, NOT clinical resistance — use `dna-decode amr` for R/S]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
