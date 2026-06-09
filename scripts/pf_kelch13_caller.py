"""Plasmodium falciparum Pfkelch13 (K13) artemisinin-resistance caller — the protozoan kingdom jump.

Calls artemisinin partial resistance from a P. falciparum genome by BLASTing the K13 CDS reference vs the
assembly, translating the aligned subject region, and checking the WHO-validated K13 propeller markers
(`dna_decode/data/antimalarial_amr.py`). Mirrors the fungal ERG11 caller exactly — and REUSES its
gene-generic BLAST+codon-mapping machinery (`observed_substitutions`, gap-aware, plus-strand HSP) rather
than re-implementing it. The only new pieces are the K13 catalog (the data module) + this thin wrapper.

K13 (PF3D7_1343700) is intronless across the propeller domain, so the CDS-vs-genome HSP is colinear and
codon-mapping is direct. Offline-safe: absent BLAST+ → INDETERMINATE with a reason (same degrade contract
as the fungal caller), so tests stay green without the binaries.

G0 (machinery) = validated against a planted C580Y (the bundled synthetic test). G0-completion = the real
3D7 K13 reference + a real C580Y isolate (the committed `data/antimalarial_ref/` fixtures, when present).
"""
from __future__ import annotations

from dna_decode.data.antimalarial_amr import AntimalarialCall, call_from_observed_substitutions
# Reuse the PROVEN gene-generic BLAST + codon-mapping from the fungal caller (DRY — it is not ERG11-specific
# despite the parameter name; it BLASTs any in-frame CDS reference vs a genome and codon-maps the best HSP).
from scripts.fungal_erg11_caller import observed_substitutions


def call_kelch13(genome_fasta: str, k13_cds_ref_fasta: str, drug: str = "artemisinin",
                 gene: str = "K13") -> AntimalarialCall:
    """Full call: BLAST K13-CDS-ref vs genome → observed substitutions → R/S vs the WHO-validated catalog.

    `k13_cds_ref_fasta` MUST be the in-frame K13 CDS (starts at ATG); protein position P ↔ CDS nt
    (3P-2..3P), 1-based, WT residue at 580 = Cysteine on the 3D7 reference."""
    obs = observed_substitutions(genome_fasta, k13_cds_ref_fasta, gene=gene)
    if obs is None:
        return AntimalarialCall("INDETERMINATE", drug, [], [], "antimalarial_k13_blastn_v0",
                                "blastn/makeblastdb not found — install BLAST+ to call from a genome")
    return call_from_observed_substitutions(drug, obs)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--k13-ref", required=True, help="in-frame K13 CDS reference FASTA (3D7)")
    ap.add_argument("--drug", default="artemisinin")
    a = ap.parse_args()
    c = call_kelch13(a.genome, a.k13_ref, a.drug)
    print(f"CALL: {c.prediction} [{c.drug}]  determinants={c.determinants}")
    if c.undetectable_mechanisms:
        print(f"  blind spots: {c.undetectable_mechanisms}")
    print(f"  {c.caveat}")
