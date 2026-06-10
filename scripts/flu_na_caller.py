"""Influenza A neuraminidase (NA) inhibitor-resistance caller — the viral kingdom jump (4th kingdom).

Calls NAI (oseltamivir / peramivir / zanamivir) resistance from an influenza A genome/segment by BLASTing
the NA CDS reference vs the assembly, translating the aligned subject region, and checking the
CDC/WHO-AVWG-recognized NA markers (`dna_decode/data/antiviral_amr.py`). Mirrors the K13 + fungal ERG11
callers exactly — and REUSES the same gene-generic BLAST+codon-mapping machinery (`observed_substitutions`)
rather than re-implementing it. The only new pieces are the NA catalog (the data module) + this thin wrapper.

NA is encoded on a non-spliced influenza segment, so the CDS-vs-genome HSP is colinear and codon-mapping is
direct (intronless, like K13 / bacterial CDS; unlike intron-containing pfcrt). Numbering is N1 (the shipped
reference NC_026434.1 has WT His at 275 = the H275Y marker). Offline-safe: absent BLAST+ → INDETERMINATE
with a reason (same degrade contract as the other callers), so tests stay green without the binaries.
"""
from __future__ import annotations

from dna_decode.data.antiviral_amr import AntiviralCall, call_from_observed_substitutions
# Reuse the PROVEN gene-generic BLAST + codon-mapping from the fungal caller (DRY — it BLASTs any in-frame
# CDS reference vs a genome and codon-maps the best HSP; not ERG11/K13-specific despite the lineage).
from scripts.fungal_erg11_caller import observed_substitutions


def call_neuraminidase(genome_fasta: str, na_cds_ref_fasta: str, drug: str = "oseltamivir",
                       gene: str = "NA") -> AntiviralCall:
    """Full call: BLAST NA-CDS-ref vs genome → observed substitutions → R/S vs the recognized NAI catalog.

    `na_cds_ref_fasta` MUST be the in-frame NA CDS (starts at ATG); protein position P ↔ CDS nt
    (3P-2..3P), 1-based, WT residue at 275 = Histidine on the N1 reference."""
    obs = observed_substitutions(genome_fasta, na_cds_ref_fasta, gene=gene)
    if obs is None:
        return AntiviralCall("INDETERMINATE", drug, [], [], "antiviral_na_blastn_v0",
                             "blastn/makeblastdb not found — install BLAST+ to call from a genome")
    return call_from_observed_substitutions(drug, obs)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--na-ref", required=True, help="in-frame N1 NA CDS reference FASTA")
    ap.add_argument("--drug", default="oseltamivir")
    a = ap.parse_args()
    c = call_neuraminidase(a.genome, a.na_ref, a.drug)
    print(f"CALL: {c.prediction} [{c.drug}]  determinants={c.determinants}")
    if c.undetectable_mechanisms:
        print(f"  blind spots: {c.undetectable_mechanisms}")
    print(f"  {c.caveat}")
