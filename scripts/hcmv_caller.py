"""HCMV target-site caller — genome-FASTA mode for the herpesvirus antiviral cell (UL97/UL54/UL56).

Calls ganciclovir/valganciclovir/cidofovir/foscarnet (UL97+UL54) or letermovir (UL56) resistance from an HCMV
genome (or a UL97/UL54/UL56 gene segment) by BLASTing each committed in-frame CDS reference vs the assembly,
codon-mapping the aligned region, and checking the catalogued mutations (`dna_decode/data/hcmv_amr.py`).
Mirrors the SARS-CoV-2 / HIV / influenza-NA callers and REUSES the same gene-generic BLAST+codon-mapping
machinery (`observed_substitutions`) — the only new piece is the MULTI-GENE loop (a drug is scored against its
target gene(s); ganciclovir needs both UL97 and UL54).

REFERENCE / NUMBERING: the shipped references are the three HCMV in-frame CDS from strain **Merlin RefSeq
NC_006273.2** — UL97 141798..143921 (fwd, 707 aa), UL54 complement(78194..81922) (1242 aa), UL56
complement(84752..87304) (850 aa) — committed at `data/hcmv_ref/HCMV_{UL97,UL54,UL56}_Merlin_cds.fna` (NCBI
already reverse-complements the complement-strand fetches, so each file is the sense-strand CDS). Merlin
numbering matches the AD169/Chou catalog convention at EVERY catalogued RESISTANCE position (asserted by
`tests/test_hcmv_caller.py`). The ONE documented strain difference is the BENIGN polymorphic site UL54:897
(Merlin encodes L897 where the AD169 catalog references S897) — it affects no resistance call and is a pinned
known-exception in the integrity test. Merlin is preferred over the high-passage AD169 lab strain (cleaner WT).

Offline-safe: absent BLAST+ -> INDETERMINATE with a reason (same degrade contract as the other callers).
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.data.hcmv_amr import HCMVCall, call_hcmv_observed, genes_for_hcmv_drug
# Reuse the PROVEN gene-generic BLAST + codon-mapping (BLASTs any in-frame CDS reference vs a genome and
# codon-maps the best HSP; the same mapper the HIV/NA/SARS-CoV-2 cells use).
from scripts.fungal_erg11_caller import observed_substitutions

_REF_DIR = Path(__file__).resolve().parent.parent / "data" / "hcmv_ref"


def default_ref_for(gene: str) -> Path:
    """Committed Merlin NC_006273.2 in-frame CDS reference for a catalogued gene."""
    return _REF_DIR / f"HCMV_{gene}_Merlin_cds.fna"


def call_hcmv_target(genome_fasta: str, drug: str, *, ref_map: dict[str, str] | None = None) -> HCMVCall:
    """Full call: for each of the DRUG's target gene(s), BLAST the in-frame CDS reference vs the genome ->
    observed substitutions; merge across genes -> R/S via the HCMV dispatcher. `ref_map` overrides the
    committed defaults per gene. Returns INDETERMINATE if BLAST+ is absent for the first gene attempted."""
    genes = genes_for_hcmv_drug(drug)
    if not genes:
        return HCMVCall("INDETERMINATE", drug, [], [], "hcmv_target_site_blastn_v0",
                        f"no HCMV catalog for {drug!r}")
    observed: dict[str, set[str]] = {}
    for gene in genes:
        ref = (ref_map or {}).get(gene) or str(default_ref_for(gene))
        obs = observed_substitutions(genome_fasta, ref, gene=gene)
        if obs is None:
            return HCMVCall("INDETERMINATE", drug, [], [], "hcmv_target_site_blastn_v0",
                            "blastn/makeblastdb not found -- install BLAST+ to call from a genome")
        observed.update(obs)   # obs is {gene: {subs}}
    return call_hcmv_observed(drug, observed)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--drug", default="ganciclovir")
    a = ap.parse_args()
    c = call_hcmv_target(a.genome, a.drug)
    print(f"CALL: {c.prediction} [{c.drug}]  determinants={c.determinants}")
    if c.undetectable_mechanisms:
        print(f"  blind spots: {c.undetectable_mechanisms}")
    print(f"  {c.caveat}")
