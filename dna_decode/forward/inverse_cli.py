"""`dna-decode inverse` — effect -> edit. Propose the edits at a target percentile of predicted damage.

The INVERSE of `dna-decode forward`. Uses the DMS-validated forward oracle as label-free ground truth: no
phenotype label is consulted, which is the move that dodges this project's binding constraint.

    dna-decode inverse --protein-seq MSIQ... --target-percentile 0.05          # the 5% most damaging
    dna-decode inverse --protein-fasta tem1.faa --target-percentile 0.5 --top-k 10
    dna-decode inverse --protein-fasta tem1.faa --cds-fasta blatem.fna         # single-nt-reachable only
    dna-decode inverse --protein-seq MSIQ... --target-percentile 0.9 --json

IT RANKS. IT DOES NOT DOSE. It proposes *edits near the top of the damaging tail*; it can never say
*fold-change 4.2*. That limit is measured, not modest -- the magnitude version needs a calibrator fit on the
TARGET protein's own DMS (which would make the inverse unnecessary), and its conformal interval is
uninformative even where it brackets. See `dna_decode/forward/inverse.py` + `wiki/forward_inverse_*`.

Measured: beats an exact no-oracle null on 4/4 proteins across 4 kingdoms (~2-5 percentile points at
top-5), graded against measured wet-lab DMS. The learned oracle beats plain BLOSUM62 on only 3/4, so the
BLOSUM62 default is often the right answer rather than a fallback.

SCOPE: Regime B (molecular fitness) ONLY -- never clinical resistance (use `dna-decode amr`, the frozen
determinant catalogue; this scorer class is below chance on antagonistically-selected resistance). NOT a
clinical tool; research use only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .inverse import EVIDENCE, UNSUPPORTED_CLAIMS, InverseError, propose_edits


def _read_fasta(path: str) -> str:
    seq = "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))
    return "".join(c for c in seq.upper() if c.isalpha())


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode inverse",
        description="Inverse design: propose the edits at a target percentile of PREDICTED molecular damage "
                    "(Regime B). Ranks, never doses.",
        epilog="scope: molecular fitness RANK, not clinical resistance (use `dna-decode amr` for R/S).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--protein-seq", help="the protein sequence")
    src.add_argument("--protein-fasta", help="FASTA holding the protein sequence")
    ap.add_argument("--target-percentile", type=float, required=True,
                    help="0.0 = most damaging predicted edit, 1.0 = most tolerated")
    ap.add_argument("--top-k", type=int, default=5,
                    help="how many edits to propose (default 5). Assay ALL of them and keep the best -- "
                         "single-shot is ~4x worse than best-of-5")
    ap.add_argument("--method", choices=["blosum62"], default="blosum62",
                    help="blosum62 (deterministic, wheel-only). The learned methods need a precomputed "
                         "score table and stay in the Python API -- and earn their keep on only 3/4 proteins")
    ap.add_argument("--cds-fasta", help="the protein's CDS -> restrict proposals to single-nt-reachable "
                                        "edits (real genome editing). Without it, all substitutions.")
    ap.add_argument("--no-diverse", action="store_true",
                    help="allow several proposals at the SAME residue (plain window). Default is one edit "
                         "per residue -- measured free for ESM and BETTER for BLOSUM (its ties are coarse)")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    a = ap.parse_args(argv)

    try:
        seq = a.protein_seq.strip().upper() if a.protein_seq else _read_fasta(a.protein_fasta)
        cds = _read_fasta(a.cds_fasta) if a.cds_fasta else None
        res = propose_edits(seq, a.target_percentile, top_k=a.top_k, method=a.method, cds=cds,
                            diverse=not a.no_diverse)
    except InverseError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    if a.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    print(f"Inverse design — edits at percentile {res.target_percentile:.2f} of PREDICTED damage")
    print(f"  method: {res.method}   candidates: {res.n_candidates} ({res.candidate_space})")
    print(f"  {'mutation':>10s} {'score':>8s} {'pctile':>7s}   codon")
    for p in res.proposals:
        print(f"  {p.mutation:>10s} {p.score:8.2f} {p.score_percentile:7.3f}   {p.codon_change or '-'}")
    for n in res.notes:
        print(f"  note: {n}")
    print(f"\n  RANK ONLY — this does NOT support:")
    for c in UNSUPPORTED_CLAIMS:
        print(f"    - {c}")
    print(f"  evidence: ESM (learned) beats a no-oracle null on {EVIDENCE['esm_rank_inverse_beats_null_at_scale']}. "
          f"The shipped blosum62 DEFAULT beats it {EVIDENCE['shipped_blosum62_default_beats_null_at_scale']}."
          f"\n  -> run scripts/forward_inverse_deployable.py on a DMS assay for YOUR protein before trusting "
          f"the default; see {EVIDENCE['artifact']}")
    print("  research use only — not a clinical tool.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
