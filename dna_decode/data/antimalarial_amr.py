"""Antimalarial (artemisinin) resistance determinant catalog — the THIRD kingdom (protozoan).

After the bacterial AMR decoder (AMRFinder → amr_rules) and the fungal kingdom jump (C. auris ERG11 →
fungal_amr), this extends the SAME proven deterministic target-site method to *Plasmodium falciparum*
artemisinin partial resistance via the **Pfkelch13 (K13)** propeller domain — the WHO-recognized molecular
marker. Like fungi, there is NO AMRFinder-equivalent for Plasmodium, so the determinant catalog is
HAND-CURATED here and the call is a BLAST(K13-CDS-vs-genome) → translate → check-codon step (reusing the
fungal caller's gene-generic `observed_substitutions`).

KEY DIFFERENCE from AMR/fungal MIC: artemisinin partial resistance is NOT defined by an MIC breakpoint —
it is a delayed parasite-clearance phenotype (clearance half-life >=5 h / RSA0-3h survival >=1%). The
**validated K13 marker IS the genotypic resistance call** (WHO uses K13 status as the surveillance
definition). So there is no `mic_to_phenotype` here; a catalogued K13 propeller substitution = ART-R.

Determinants = the WHO-VALIDATED Pfkelch13 artemisinin-partial-resistance markers (WHO "Report on
antimalarial drug efficacy, resistance and response", validated list; Ariey 2014 + WWARN K13 surveyor).
C580Y is the canonical/most-widespread (Greater Mekong); R561H (Rwanda/East Africa), A675V/R622I (East
Africa/Horn) are the African-emergence markers. Numbering is 1-based on the 726-aa 3D7 K13 reference
(PF3D7_1343700), WT residue at 580 = Cysteine (so the marker is C580Y).
"""
from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Resistance-determinant catalog: drug -> gene -> set of resistance amino-acid substitutions.
# Substitution form "<wt><pos><mut>" (1-letter, 1-based on the 3D7 K13 reference).
# ---------------------------------------------------------------------------

# WHO-VALIDATED Pfkelch13 artemisinin partial-resistance markers (propeller-domain substitutions).
_K13_VALIDATED = {
    "F446I", "N458Y", "M476I", "Y493H", "R539T", "I543T", "P553L",
    "R561H", "P574L", "C580Y", "R622I", "A675V",
}

ANTIMALARIAL_RESISTANCE_MUTATIONS: dict[str, dict[str, set[str]]] = {
    # The artemisinin component of every ACT (and artesunate monotherapy) shares the K13 target.
    "artemisinin": {"K13": set(_K13_VALIDATED)},
    "artesunate": {"K13": set(_K13_VALIDATED)},
    "dihydroartemisinin": {"K13": set(_K13_VALIDATED)},
    # Chloroquine: pfcrt K76T is THE validated molecular marker — the necessary (near-sufficient)
    # determinant of P. falciparum chloroquine resistance (WHO/WWARN; the CVIET/SVMNT haplotypes all carry
    # the 76T substitution). Direction is unambiguous (K76T -> CQ-R), unlike the pfmdr1 partner-drug markers
    # (N86Y selects OPPOSITELY for amodiaquine vs lumefantrine) — those are deliberately NOT catalogued here.
    "chloroquine": {"pfcrt": {"K76T"}},
}

# Single target gene per drug (the catalog has one gene per drug in v0). Used to route the caller +
# the genome-mode intron guard (pfcrt is intron-containing; K13 is intronless).
def gene_for_drug(drug: str) -> str | None:
    genes = ANTIMALARIAL_RESISTANCE_MUTATIONS.get(drug.lower())
    return next(iter(genes), None) if genes else None

# Genes whose genomic locus is INTRONLESS, so the colinear single-HSP BLAST codon-mapper is valid in
# genome mode. K13 is intronless (genome mode OK); pfcrt has 13 exons (GenBank seqs are ~2471bp genomic) —
# genome mode for pfcrt needs intron-aware multi-HSP stitching (DEFERRED), so chloroquine is --observed-only
# for now (the catalog/--observed path needs no reference + no alignment, so it is unaffected by introns).
INTRONLESS_GENES = {"K13"}

# Mechanisms a K13-only target-site scan CANNOT see (the protozoan analogue of the bacterial efflux /
# fungal aneuploidy blind spot). A SUSCEPTIBLE (no-K13-marker) call cannot rule these out.
ANTIMALARIAL_UNDETECTABLE_MECHANISMS = sorted({
    "non_K13_background_mutations",        # pffd/pfarps10/pfmdr2/pfcrt background that modulates K13-R
    "pfcoronin_artemisinin_resistance",    # K13-independent ART-R (in-vitro selected)
    "partner_drug_resistance",             # pfmdr1 CNV (lumefantrine/mefloquine), pfcrt (CQ),
                                           # plasmepsin2/3 CNV (piperaquine) — ACT failure not via K13
    "novel_uncatalogued_K13_propeller_substitution",
})

ANTIMALARIAL_DRUG_CLASS = {
    "artemisinin": "ARTEMISININ", "artesunate": "ARTEMISININ", "dihydroartemisinin": "ARTEMISININ",
    "chloroquine": "4-AMINOQUINOLINE",
}


def supported_antimalarial_drugs() -> list[str]:
    return sorted(ANTIMALARIAL_RESISTANCE_MUTATIONS)


def resistance_mutations_for(drug: str) -> dict[str, set[str]]:
    """gene -> resistance-substitution set for an antimalarial drug. Raises KeyError on unknown drug."""
    d = drug.lower()
    if d not in ANTIMALARIAL_RESISTANCE_MUTATIONS:
        raise KeyError(f"no antimalarial resistance catalog for drug {drug!r}; "
                       f"configured: {supported_antimalarial_drugs()}")
    return ANTIMALARIAL_RESISTANCE_MUTATIONS[d]


def is_resistance_mutation(drug: str, gene: str, substitution: str) -> bool:
    """True iff <substitution> (e.g. 'C580Y') in <gene> is a validated ART-R marker for <drug>."""
    return substitution in ANTIMALARIAL_RESISTANCE_MUTATIONS.get(drug.lower(), {}).get(gene, set())


@dataclass(frozen=True)
class AntimalarialCall:
    """A deterministic antimalarial R/S call (shape mirrors FungalCall / amr_rules output)."""
    prediction: str            # R / S / INDETERMINATE
    drug: str
    determinants: list[str]    # e.g. ["K13:C580Y"]
    undetectable_mechanisms: list[str]
    rule: str
    caveat: str


def call_from_observed_substitutions(drug: str, observed: dict[str, set[str]]) -> AntimalarialCall:
    """Deterministic artemisinin R/S call from a genome's observed K13 substitutions.

    `observed` = {gene: {substitutions present}} (from the BLAST caller). Rule (v0): R iff the genome
    carries >=1 WHO-validated K13 propeller marker. An S call surfaces ANTIMALARIAL_UNDETECTABLE_MECHANISMS
    — it means 'no validated K13 marker found', NOT 'definitely artemisinin-sensitive'.
    """
    if drug.lower() not in ANTIMALARIAL_RESISTANCE_MUTATIONS:
        return AntimalarialCall("INDETERMINATE", drug, [], [], "antimalarial_k13_target_mutation_v0",
                                f"no catalog for {drug!r}")
    hits = []
    for gene, subs in resistance_mutations_for(drug).items():
        for s in observed.get(gene, set()):
            if s in subs:
                hits.append(f"{gene}:{s}")
    pred = "R" if hits else "S"
    undetectable = ANTIMALARIAL_UNDETECTABLE_MECHANISMS if pred == "S" else []
    gene = gene_for_drug(drug) or "target"
    phen = ("Artemisinin partial resistance is a clearance phenotype, not an MIC — the validated K13 "
            "marker IS the genotypic call. " if gene == "K13"
            else f"The validated {gene} marker IS the genotypic resistance call. ")
    caveat = (f"deterministic {gene} target-site call ({ANTIMALARIAL_DRUG_CLASS.get(drug.lower(), '?')}). "
              "Hand-curated WHO/WWARN-validated marker catalog (no AMRFinder-equivalent for Plasmodium). "
              + phen
              + ("An S call cannot rule out non-target / partner-drug resistance "
                 f"({', '.join(ANTIMALARIAL_UNDETECTABLE_MECHANISMS)})." if pred == "S" else ""))
    return AntimalarialCall(pred, drug, sorted(hits), undetectable, "antimalarial_k13_target_mutation_v0",
                            caveat)
