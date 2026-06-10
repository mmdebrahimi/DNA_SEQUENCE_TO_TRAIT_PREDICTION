"""Fungal ERG11/FKS1 target-site mutation caller (EP-7 step 2 → Gate G0).

Calls azole/echinocandin resistance from a fungal genome by BLASTing a target-gene CDS reference against
the assembly, translating the aligned subject region, and checking catalogued resistance substitutions
(`dna_decode/data/fungal_amr.py`). The eukaryotic analogue of the AMRFinder→`amr_rules` step — but since
there is NO AMRFinder-for-fungi, the reference allele + catalog are supplied here.

Uses `blastn` (CDS-vs-genome). C. auris ERG11 is intronless, so the HSP is colinear and codon-mapping is
direct (gap-aware: a position interrupted by an indel is reported uncalled, never mis-translated). `tblastn`
would be cleaner (protein-vs-genome) but is not installed; `blastn` + `makeblastdb` are
(`C:/Users/Farshad/ncbi-blast/bin`). Offline-safe: absent BLAST → INDETERMINATE with a reason (mirrors the
pathotype `vf_runner` degrade pattern), so tests stay green without the binaries.

G0 = this machinery is validated against a KNOWN planted mutation (the bundled test). G0-completion (EP-7
step 3) swaps the synthetic reference for the real C. auris ERG11 allele + validates on a real
resistant genome (e.g. a Y132F / VF125AL isolate) to confirm the catalog numbering matches the reference.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from dna_decode.data.fungal_amr import FungalCall, call_from_observed_substitutions

_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_BLAST_DIRS = [os.environ.get("NCBI_BLAST_BIN", ""), "", str(Path.home() / "ncbi-blast" / "bin"),
               "C:/Users/Farshad/ncbi-blast/bin"]


def _find(tool: str) -> str | None:
    for d in _BLAST_DIRS:
        cand = shutil.which(tool, path=d) if d else shutil.which(tool)
        if cand:
            return cand
    return None


def _translate(seq: str) -> str:
    seq = seq.upper().replace("-", "")
    return "".join(_CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def observed_substitutions(genome_fasta: str, erg11_cds_ref_fasta: str,
                           gene: str = "ERG11") -> dict[str, set[str]] | None:
    """BLAST the in-frame CDS reference vs the genome; return {gene: {observed substitutions}}.

    INTRON-AWARE (multi-HSP): the query (CDS reference) is contiguous, so codon position P always maps to
    CDS nt (3P-2..3P) regardless of how the SUBJECT genome splits the gene into exons. For an
    intron-containing gene, blastn returns one HSP per exon (colinear within an exon, separated by introns
    in the subject); we STITCH the query-position → subject-nucleotide map across all HSPs on the gene's
    contig. A codon spanning an exon-exon boundary (its 3 nts in two different HSPs) is therefore still
    translated correctly — the 3 coding nts are contiguous in the query though physically separated by an
    intron in the genome. Intronless genes (ERG11, K13) are the single-HSP special case → identical to the
    prior behavior (guarded by their existing tests).

    Returns None if BLAST is unavailable (caller maps that to INDETERMINATE). The reference FASTA MUST be
    the in-frame CDS (starts at ATG); protein position P ↔ CDS nt (3P-2..3P), 1-based.
    """
    blastn, makedb = _find("blastn"), _find("makeblastdb")
    if not blastn or not makedb:
        return None
    ref_seq = _read_single_fasta(erg11_cds_ref_fasta)
    ref_prot = _translate(ref_seq)
    with tempfile.TemporaryDirectory() as td:
        dbpath = str(Path(td) / "genome")
        subprocess.run([makedb, "-in", genome_fasta, "-dbtype", "nucl", "-out", dbpath],
                       check=True, capture_output=True)
        out = subprocess.run(
            [blastn, "-query", erg11_cds_ref_fasta, "-db", dbpath,
             "-outfmt", "6 sseqid qstart qend sstart send qseq sseq bitscore", "-max_target_seqs", "10"],
            check=True, capture_output=True, text=True).stdout
    if not out.strip():
        return {gene: set()}
    rows = [ln.split("\t") for ln in out.strip().splitlines() if ln.strip()]
    # Restrict to the single best CONTIG (highest total HSP bitscore) — the gene's genomic location —
    # so paralog/repeat HSPs on other contigs don't pollute the map.
    by_contig: dict[str, list] = {}
    for r in rows:
        by_contig.setdefault(r[0], []).append(r)
    best_contig = max(by_contig, key=lambda c: sum(float(r[7]) for r in by_contig[c]))
    hsps = by_contig[best_contig]
    # Stitch query-CDS 1-based position -> subject aligned nucleotide across all HSPs (gap-aware).
    # Highest-bitscore HSP wins on any query-position overlap (a weaker HSP never overwrites).
    qpos_to_snt: dict[int, str] = {}
    for r in sorted(hsps, key=lambda r: -float(r[7])):
        qstart, qend, qseq, sseq = int(r[1]), int(r[2]), r[5], r[6]
        if qstart > qend:        # query (CDS) should be plus-strand; skip exotic orientations
            continue
        qp = qstart
        for qc, sc in zip(qseq, sseq):
            if qc == "-":        # insertion in subject relative to ref CDS — advance subject only
                continue
            if qp not in qpos_to_snt:
                qpos_to_snt[qp] = sc     # sc may be '-' (deletion); handled below
            qp += 1
    subs: set[str] = set()
    n_codons = len(ref_prot)
    for p in range(1, n_codons + 1):
        # codon callable iff all 3 query positions are covered by SOME HSP (intron-spanning OK)
        nts = [qpos_to_snt.get(3 * p - 2), qpos_to_snt.get(3 * p - 1), qpos_to_snt.get(3 * p)]
        if any(n is None or n == "-" for n in nts):
            continue             # uncovered or indel-interrupted codon → uncalled (never mis-translate)
        q_res = _CODON.get("".join(nts).upper(), "X")
        ref_res = ref_prot[p - 1]
        if q_res not in ("X", "*") and q_res != ref_res:
            subs.add(f"{ref_res}{p}{q_res}")
    return {gene: subs}


def call_erg11(genome_fasta: str, erg11_cds_ref_fasta: str, drug: str = "fluconazole",
               gene: str = "ERG11") -> FungalCall:
    """Full call: BLAST → observed substitutions → deterministic R/S vs the fungal catalog."""
    obs = observed_substitutions(genome_fasta, erg11_cds_ref_fasta, gene=gene)
    if obs is None:
        return FungalCall("INDETERMINATE", drug, [], [], "fungal_erg11_blastn_v0",
                          "blastn/makeblastdb not found — install BLAST+ or use --eval-only with cached calls")
    return call_from_observed_substitutions(drug, obs)


def _read_single_fasta(path: str) -> str:
    seq = []
    for ln in Path(path).read_text().splitlines():
        if not ln.startswith(">"):
            seq.append(ln.strip())
    return "".join(seq)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--erg11-ref", required=True, help="in-frame ERG11 CDS reference FASTA")
    ap.add_argument("--drug", default="fluconazole")
    a = ap.parse_args()
    c = call_erg11(a.genome, a.erg11_ref, a.drug)
    print(f"CALL: {c.prediction} [{c.drug}]  determinants={c.determinants}")
    if c.undetectable_mechanisms:
        print(f"  blind spots: {c.undetectable_mechanisms}")
    print(f"  {c.caveat}")
