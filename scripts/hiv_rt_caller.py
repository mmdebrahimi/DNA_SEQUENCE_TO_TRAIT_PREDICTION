"""HIV-1 reverse-transcriptase (RT) target-site caller — genome-FASTA mode for the HIV viral cell.

Calls NNRTI / NRTI resistance from an HIV-1 genome (or a pol/RT segment) by BLASTing the in-frame RT CDS
reference vs the assembly, translating the aligned subject region, and checking the established major DRMs
(`dna_decode/data/hiv_amr.py`). Mirrors the influenza NA + fungal ERG11 + antimalarial K13 callers exactly,
and REUSES the same gene-generic BLAST+codon-mapping machinery (`observed_substitutions`) rather than
re-implementing it — the only new pieces are the HIV RT catalog (the data module, already shipped) + this
thin wrapper + the committed reference.

REFERENCE / NUMBERING: the shipped reference is the HXB2 RT p66 in-frame CDS (NCBI K03455.1:2550-4229,
560 codons) at `data/hiv_ref/HIV1_RT_HXB2_cds.fna`. RT residue 1 = mature-RT Pro (N-terminus PISP...), so
protein position P <-> CDS nt (3P-2..3P), 1-based, in consensus-B RT numbering. HXB2 matches consensus-B WT
at EVERY catalogued NNRTI+NRTI DRM position (asserted by `tests/test_hiv_rt_caller.py::test_reference_*`),
which is what makes a substitution call comparable to the Stanford-sourced catalog. RT is encoded on the
unspliced pol gene, so the CDS-vs-genome HSP is colinear and codon-mapping is direct (intronless, like NA /
K13 / bacterial CDS).

ROUTING: NNRTI drugs -> `call_from_observed_substitutions`; NRTI drugs -> `call_nrti_from_observed` (the v0
position-based call). `is_nrti` picks the branch (the CLI already knows the class from `--drug`).

Offline-safe: absent BLAST+ -> INDETERMINATE with a reason (same degrade contract as the other callers), so
tests stay green without the binaries.
"""
from __future__ import annotations

from dna_decode.data.hiv_amr import HIVCall, call_hiv_observed
# Reuse the PROVEN gene-generic BLAST + codon-mapping from the fungal caller (DRY — it BLASTs any in-frame
# CDS reference vs a genome and codon-maps the best HSP; not ERG11/NA/K13-specific despite the lineage).
from scripts.fungal_erg11_caller import observed_substitutions


def call_hiv_target(genome_fasta: str, cds_ref_fasta: str, drug: str, gene: str) -> HIVCall:
    """Full call: BLAST the in-frame CDS reference for `gene` (RT/PR/IN/CA) vs the genome -> observed
    substitutions -> R/S via the unified HIV dispatcher (NNRTI/NRTI/PI/INSTI/CAI, routed by drug).

    `cds_ref_fasta` MUST be the in-frame CDS for `gene` (HXB2 K03455.1 references ship as defaults); protein
    position P <-> CDS nt (3P-2..3P), 1-based, consensus-B numbering."""
    obs = observed_substitutions(genome_fasta, cds_ref_fasta, gene=gene)
    if obs is None:
        return HIVCall("INDETERMINATE", drug, [], [], "hiv_target_blastn_v0",
                       "blastn/makeblastdb not found — install BLAST+ to call from a genome")
    return call_hiv_observed(drug, obs)


def call_hiv_rt(genome_fasta: str, rt_cds_ref_fasta: str, drug: str = "efavirenz",
                is_nrti: bool = False, gene: str = "RT") -> HIVCall:
    """Back-compat RT wrapper (NNRTI/NRTI). `is_nrti` is vestigial — the unified dispatcher routes by drug."""
    return call_hiv_target(genome_fasta, rt_cds_ref_fasta, drug, gene)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--rt-ref", required=True, help="in-frame HIV-1 RT CDS reference FASTA")
    ap.add_argument("--drug", default="efavirenz")
    ap.add_argument("--nrti", action="store_true", help="route to the NRTI position-based call")
    a = ap.parse_args()
    c = call_hiv_rt(a.genome, a.rt_ref, a.drug, is_nrti=a.nrti)
    print(f"CALL: {c.prediction} [{c.drug}]  determinants={c.determinants}")
    if c.undetectable_mechanisms:
        print(f"  blind spots: {c.undetectable_mechanisms}")
    print(f"  {c.caveat}")
