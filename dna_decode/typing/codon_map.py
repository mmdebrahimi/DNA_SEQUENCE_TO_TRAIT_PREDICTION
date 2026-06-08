"""Codon translation + gap-aware codon→subject-amino-acid mapping from a blastn HSP.

Shared by point-mutation callers (PointFinder chromosomal AMR; the fungal ERG11 caller uses the same logic).
Given a blastn HSP of an in-frame CDS reference (query) vs a genome (subject), recover the SUBJECT amino acid
at each 1-based reference codon position — gap-aware: a codon interrupted by an indel is left UNCALLED
(never mis-translated). The reference must be the in-frame CDS (codon P ↔ query nt 3P-2..3P, 1-based).
"""
from __future__ import annotations

CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str) -> str:
    s = seq.upper().replace("-", "")
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def subject_aa_by_codon(qseq: str, sseq: str, qstart: int, n_codons: int) -> dict[int, str]:
    """{1-based reference codon position -> subject amino acid} for the aligned region (gap-aware).

    qseq/sseq are the aligned query(=ref CDS)/subject strings from a blastn HSP (`-outfmt '... qseq sseq'`),
    qstart the 1-based query start. Codons with an indel-interrupted or out-of-HSP position are omitted.
    """
    # map query-CDS 1-based nt position -> subject aligned nt (gap-aware: skip subject insertions)
    qpos_to_snt: dict[int, str] = {}
    qp = qstart
    qend = qstart - 1
    for qc, sc in zip(qseq, sseq):
        if qc == "-":            # insertion in subject relative to ref -> advance subject only
            continue
        qpos_to_snt[qp] = sc     # sc may be '-' (deletion); handled below
        qend = qp
        qp += 1
    out: dict[int, str] = {}
    for p in range(1, n_codons + 1):
        c1, c2, c3 = 3 * p - 2, 3 * p - 1, 3 * p
        if not (qstart <= c1 and c3 <= qend):
            continue
        nts = [qpos_to_snt.get(c1), qpos_to_snt.get(c2), qpos_to_snt.get(c3)]
        if any(n is None or n == "-" for n in nts):
            continue             # indel-interrupted codon -> uncalled
        aa = CODON.get("".join(nts).upper(), "X")
        if aa not in ("X",):
            out[p] = aa
    return out
