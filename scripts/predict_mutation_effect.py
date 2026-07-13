"""CLI for the rung-2 mutation-effect predictor (2026-07-13).

Point at a protein + a single point mutation -> the deterministic sequence change + a molecular-effect RANK
score (ESM2-650M damage_llr) with honesty rails. v0 requires a SUPPLIED sequence (no gene/accession network
lookup — avoids strain/numbering ambiguity + keeps it reproducible/offline).

  # supply the sequence inline or from a FASTA-ish file (first non-header line used)
  uv run python scripts/predict_mutation_effect.py --seq MAQ...GV --mutation A10V
  uv run python scripts/predict_mutation_effect.py --seq-file prot.fasta --mutation A10V

Face-validity anchor (the honest accuracy number): scripts/protein_effect_facevalidity.py ->
ARGR_ECOLI stability Spearman ~0.50 (matches ProteinGym ESM2-650M ~0.52 stability). Frozen decoder
surface untouched.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.protein_effect import gene_lookup as G  # noqa: E402
from dna_decode.protein_effect import predictor as P  # noqa: E402

CACHE_DIR = REPO / "data" / "processed" / "protein_effect_cache"


def _read_seq(seq: str | None, seq_file: str | None) -> str:
    if seq:
        return seq.strip().upper()
    lines = [l.strip() for l in open(seq_file, encoding="utf-8") if l.strip() and not l.startswith(">")]
    return "".join(lines).upper()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--seq", help="protein sequence inline")
    g.add_argument("--seq-file", help="FASTA-ish file; first non-header lines concatenated")
    g.add_argument("--gene", help="resolve an E. coli K-12 gene name via UniProt (e.g. gyrA); adds a network dep")
    ap.add_argument("--organism-id", default=G.ECOLI_K12["organism_id"],
                    help="UniProt organism_id for --gene (default 83333 = E. coli K-12)")
    ap.add_argument("--mutation", required=True, help="single substitution, e.g. A123V (1-based)")
    ap.add_argument("--cache", default=None, help="masked-marginal cache JSON (reuse across mutations of one protein)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    gene_rec = None
    if a.gene:
        gene_rec = G.fetch_protein_sequence(a.gene, organism_id=a.organism_id)   # raises on ambiguity/unreachable
        seq = gene_rec["sequence"].upper()
        print(f"[gene] {a.gene} -> {gene_rec['accession']} ({gene_rec['organism']}), {gene_rec['length']} aa "
              f"[{'reviewed' if gene_rec['reviewed'] else 'UNREVIEWED'}]")
    else:
        seq = _read_seq(a.seq, a.seq_file)
    P.apply_edit(seq, a.mutation)                    # validate the edit up-front (clear error if bad)
    import hashlib
    cache = a.cache or (CACHE_DIR / (f"{a.organism_id}_{a.gene.lower()}.json" if a.gene
                                     else f"seq_{hashlib.sha256(seq.encode()).hexdigest()[:12]}.json"))
    logp = P.masked_marginals(seq, cache_path=cache, progress=True)
    out = P.predict(seq, a.mutation, logp)
    if gene_rec is not None:
        out["gene_resolution"] = {k: gene_rec[k] for k in ("gene", "accession", "organism", "reviewed",
                                                           "protein_name", "provenance")}
    txt = json.dumps(out, indent=2)
    if a.out:
        a.out.write_text(txt, encoding="utf-8")
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
