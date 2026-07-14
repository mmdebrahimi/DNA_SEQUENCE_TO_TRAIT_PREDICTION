"""Genome-level forward edit: a single nucleotide change in a coding sequence (CDS) -> the codon it hits ->
the amino-acid change -> the protein-level phenotype predictor (variant_effect.predict_effect).

This lifts the forward cell's INPUT from a protein point-mutation ('M69L') to a real genome edit
(CDS position + ref base + alt base), classifying it as SILENT (synonymous), NONSENSE (premature stop),
or MISSENSE, and routing missense/nonsense to the Regime-B predictor. Reference base is verified against
the CDS (a mismatch fails LOUDLY — the coordinate-integrity discipline), and if a protein sequence is
supplied the translated codon must match it.
"""
from __future__ import annotations

from dataclasses import dataclass

from .variant_effect import ForwardPrediction, predict_effect

# Standard genetic code (NCBI table 1; internal-codon translation is identical to bacterial table 11 —
# they differ only in alternative START codons, which do not affect a point-substitution's residue call).
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


def translate_codon(codon: str) -> str:
    """Standard-code translation of one 3-nt codon -> one-letter AA ('*' = stop). Raises on a bad codon."""
    c = codon.upper()
    if len(c) != 3 or any(b not in "ACGT" for b in c):
        raise ValueError(f"not a valid DNA codon: {codon!r}")
    return _CODON[c]


@dataclass
class GenomeEditPrediction:
    nt_pos: int                 # 1-based CDS position of the edited base
    ref_base: str
    alt_base: str
    aa_pos: int                 # 1-based residue the codon encodes
    wt_aa: str
    alt_aa: str
    wt_codon: str
    alt_codon: str
    consequence: str            # "silent" | "missense" | "nonsense"
    aa_mutation: str | None     # e.g. "M69L" (None for silent)
    protein_prediction: ForwardPrediction | None   # None for silent (no protein change)

    def as_dict(self) -> dict:
        d = {
            "nt_pos": self.nt_pos, "ref_base": self.ref_base, "alt_base": self.alt_base,
            "aa_pos": self.aa_pos, "wt_aa": self.wt_aa, "alt_aa": self.alt_aa,
            "wt_codon": self.wt_codon, "alt_codon": self.alt_codon, "consequence": self.consequence,
            "aa_mutation": self.aa_mutation,
            "protein_prediction": (self.protein_prediction.as_dict() if self.protein_prediction else None),
        }
        return d


def cds_point_edit(cds: str, nt_pos: int, ref_base: str, alt_base: str) -> dict:
    """Resolve a 1-based CDS base substitution to its codon consequence. Verifies ref_base against the CDS."""
    if nt_pos < 1 or nt_pos > len(cds):
        raise ValueError(f"nt_pos {nt_pos} out of range for CDS length {len(cds)}")
    ref_base, alt_base = ref_base.upper(), alt_base.upper()
    if alt_base not in "ACGT" or ref_base not in "ACGT":
        raise ValueError(f"ref/alt must be single DNA bases; got {ref_base!r}->{alt_base!r}")
    idx = nt_pos - 1
    have = cds[idx].upper()
    if have != ref_base:
        raise ValueError(f"REF mismatch at CDS pos {nt_pos}: sequence has {have!r}, edit asserts {ref_base!r} "
                         f"(coordinate error — refusing)")
    codon_no = idx // 3                       # 0-based codon index
    within = idx % 3                          # 0..2 position within codon
    cstart = codon_no * 3
    wt_codon = cds[cstart:cstart + 3].upper()
    if len(wt_codon) != 3:
        raise ValueError(f"edit at nt_pos {nt_pos} falls in an incomplete terminal codon (CDS not a multiple of 3?)")
    alt_codon = wt_codon[:within] + alt_base + wt_codon[within + 1:]
    return {
        "aa_pos": codon_no + 1, "within_codon": within,
        "wt_codon": wt_codon, "alt_codon": alt_codon,
        "wt_aa": translate_codon(wt_codon), "alt_aa": translate_codon(alt_codon),
    }


def predict_genome_edit(cds: str, nt_pos: int, ref_base: str, alt_base: str, *,
                        protein_seq: str | None = None, protein: str = "protein",
                        phenotype_axis: str = "molecular fitness (DMS-measured)",
                        method: str = "blosum62", esm_table: dict | None = None) -> GenomeEditPrediction:
    """Full genome-edit -> phenotype path: CDS base substitution -> codon -> AA change -> Regime-B predictor.

    - SILENT (synonymous): no protein change -> protein_prediction=None (predicted neutral at the protein level).
    - NONSENSE / MISSENSE: build the AA mutation and delegate to predict_effect (which re-verifies the WT
      residue against `protein_seq` if supplied — double coordinate check: CDS ref base AND translated WT AA).
    """
    info = cds_point_edit(cds, nt_pos, ref_base, alt_base)
    aa_pos, wt_aa, alt_aa = info["aa_pos"], info["wt_aa"], info["alt_aa"]

    if wt_aa == alt_aa:
        return GenomeEditPrediction(nt_pos, ref_base.upper(), alt_base.upper(), aa_pos, wt_aa, alt_aa,
                                    info["wt_codon"], info["alt_codon"], "silent", None, None)

    consequence = "nonsense" if alt_aa == "*" else "missense"
    aa_mut = f"{wt_aa}{aa_pos}{'*' if alt_aa == '*' else alt_aa}"
    pred = predict_effect(protein_seq or "", aa_mut, protein=protein, phenotype_axis=phenotype_axis,
                          method=method, esm_table=esm_table)
    return GenomeEditPrediction(nt_pos, ref_base.upper(), alt_base.upper(), aa_pos, wt_aa, alt_aa,
                                info["wt_codon"], info["alt_codon"], consequence, aa_mut, pred)
