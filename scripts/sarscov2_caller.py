"""SARS-CoV-2 Mpro (3CLpro) target-site caller — genome-FASTA mode for the coronavirus antiviral cell.

Calls Mpro-inhibitor resistance (nirmatrelvir / ensitrelvir / lufotrelvir) from a SARS-CoV-2 genome (or an
ORF1ab / Mpro segment) by BLASTing the in-frame Mpro CDS reference vs the assembly, translating the aligned
subject region, and checking the catalogued major Mpro substitutions (`dna_decode/data/sarscov2_amr.py`).
Mirrors the HIV RT + influenza NA + fungal ERG11 callers exactly and REUSES the same gene-generic
BLAST+codon-mapping machinery (`observed_substitutions`) — the only new pieces are the Mpro catalog (the
data module) + this thin wrapper + the committed reference.

REFERENCE / NUMBERING: the shipped reference is the SARS-CoV-2 Mpro (nsp5) in-frame CDS (Wuhan-Hu-1
NC_045512.2:10055-10972, 306 codons) at `data/sarscov2_ref/SARSCoV2_Mpro_NC045512_cds.fna`. Mpro residue 1 =
Ser (N-terminus SGFRKM...), so protein position P <-> CDS nt (3P-2..3P), 1-based. The reference matches the
Mpro WT at every catalogued position (asserted by `tests/test_sarscov2_caller.py`; catalytic dyad H41/C145
+ E166 verified). Mpro is encoded contiguously within ORF1ab (before the -1 frameshift), so the CDS-vs-genome
HSP is colinear and codon-mapping is direct (intronless, like the NA / HIV-RT / bacterial CDS path).

Offline-safe: absent BLAST+ -> INDETERMINATE with a reason (same degrade contract as the other callers).
"""
from __future__ import annotations

from dna_decode.data.sarscov2_amr import SARSCoV2Call, call_sarscov2_observed
# Reuse the PROVEN gene-generic BLAST + codon-mapping from the fungal caller (DRY — BLASTs any in-frame CDS
# reference vs a genome and codon-maps the best HSP; the same mapper the HIV/NA cells use).
from scripts.fungal_erg11_caller import observed_substitutions


def call_sarscov2_target(genome_fasta: str, cds_ref_fasta: str, drug: str, gene: str = "Mpro") -> SARSCoV2Call:
    """Full call: BLAST the in-frame Mpro CDS reference vs the genome -> observed substitutions -> R/S via
    the Mpro dispatcher. `cds_ref_fasta` MUST be the in-frame Mpro CDS (the committed Wuhan-Hu-1 reference
    ships as the default); protein position P <-> CDS nt (3P-2..3P), 1-based."""
    obs = observed_substitutions(genome_fasta, cds_ref_fasta, gene=gene)
    if obs is None:
        return SARSCoV2Call("INDETERMINATE", drug, [], [], "sarscov2_mpro_blastn_v0",
                            "blastn/makeblastdb not found — install BLAST+ to call from a genome")
    return call_sarscov2_observed(drug, obs)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--mpro-ref", required=True, help="in-frame SARS-CoV-2 Mpro CDS reference FASTA")
    ap.add_argument("--drug", default="nirmatrelvir")
    a = ap.parse_args()
    c = call_sarscov2_target(a.genome, a.mpro_ref, a.drug)
    print(f"CALL: {c.prediction} [{c.drug}]  determinants={c.determinants}")
    if c.undetectable_mechanisms:
        print(f"  blind spots: {c.undetectable_mechanisms}")
    print(f"  {c.caveat}")
