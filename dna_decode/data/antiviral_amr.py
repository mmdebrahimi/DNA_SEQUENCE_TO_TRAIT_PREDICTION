"""Antiviral (influenza neuraminidase-inhibitor) resistance determinant catalog — the FOURTH kingdom (viral).

After the bacterial AMR decoder (AMRFinder -> amr_rules), the fungal kingdom jump (C. auris ERG11 ->
fungal_amr), and the protozoan kingdom (P. falciparum K13/pfcrt -> antimalarial_amr), this extends the SAME
proven deterministic target-site method one domain further: to *influenza A virus* neuraminidase-inhibitor
(NAI) resistance via the **neuraminidase (NA)** active-site framework residues — the WHO/CDC-recognized
molecular markers. Like fungi + Plasmodium, there is NO AMRFinder-equivalent for influenza, so the
determinant catalog is HAND-CURATED here and the call is a BLAST(NA-CDS-vs-genome) -> translate ->
check-codon step (reusing the fungal caller's gene-generic `observed_substitutions`).

KEY DIFFERENCE from bacterial/fungal MIC: NAI resistance is NOT defined by a growth MIC — it is a phenotypic
**IC50 fold-change** in a neuraminidase-inhibition (NI) enzyme assay (WHO AVWG "highly/normal/reduced
inhibition" tiers). The **validated NA marker IS the genotypic resistance call** (CDC/WHO surveillance uses
NA marker status). So there is no `mic_to_phenotype` here; a catalogued NA active-site substitution = R.

NUMBERING: N1 numbering (NOT the N2 convention sometimes seen in the literature — H275Y here = "H274Y" in
N2 numbering). The reference is RefSeq NC_026434.1 (A/California/07/2009(H1N1)pdm09 NA, 470 aa, WT His at
275), so the catalog numbering is internally consistent with the shipped reference. NA is encoded on a
non-spliced influenza segment -> INTRONLESS (genome-mode single-HSP codon-mapping is valid).

Markers = CDC/WHO-AVWG-recognized influenza A(H1N1) NAI resistance substitutions (N1 numbering). H275Y is
the canonical, globally dominant oseltamivir/peramivir marker (highly reduced inhibition); I223R + S247N are
recognized accessory oseltamivir markers; Q136K + E119G are zanamivir-affecting framework markers. Direction
is unambiguous (each substitution REDUCES inhibitor binding).
"""
from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Resistance-determinant catalog: drug -> gene -> set of resistance amino-acid substitutions.
# Substitution form "<wt><pos><mut>" (1-letter, 1-based on the N1 NA reference, WT His at 275).
# ---------------------------------------------------------------------------

ANTIVIRAL_RESISTANCE_MUTATIONS: dict[str, dict[str, set[str]]] = {
    # Oseltamivir (Tamiflu): H275Y is THE dominant N1 marker (highly reduced inhibition); I223R + S247N
    # are recognized accessory/secondary oseltamivir markers in N1.
    "oseltamivir": {"NA": {"H275Y", "I223R", "S247N"}},
    # Peramivir shares the same active-site framework; H275Y confers cross-resistance.
    "peramivir": {"NA": {"H275Y", "I223R"}},
    # Zanamivir (Relenza) binds differently; H275Y barely affects it. The N1 zanamivir-affecting framework
    # markers are Q136K + E119G (and I223R contributes a mild reduction).
    "zanamivir": {"NA": {"Q136K", "E119G", "I223R"}},
}

# WT framework residues at the catalogued positions (N1 numbering on NC_026434.1) — documentation only.
_N1_WT = {119: "E", 136: "Q", 223: "I", 247: "S", 275: "H"}

ANTIVIRAL_DRUG_CLASS = {
    "oseltamivir": "NEURAMINIDASE_INHIBITOR",
    "peramivir": "NEURAMINIDASE_INHIBITOR",
    "zanamivir": "NEURAMINIDASE_INHIBITOR",
}


# Single target gene per drug (NA for every NAI in v0). Used to route the caller + pick the reference.
def gene_for_drug(drug: str) -> str | None:
    genes = ANTIVIRAL_RESISTANCE_MUTATIONS.get(drug.lower())
    return next(iter(genes), None) if genes else None


# NA is on a non-spliced influenza segment -> intronless, so the colinear single-HSP BLAST codon-mapper is
# valid in genome mode (the same as bacterial CDS / K13; unlike intron-containing pfcrt).
INTRONLESS_GENES = {"NA"}

# Mechanisms an NA-only target-site scan CANNOT see (the viral analogue of the bacterial efflux / fungal
# aneuploidy / protozoan partner-drug blind spot). A SUSCEPTIBLE (no-marker) call cannot rule these out.
ANTIVIRAL_UNDETECTABLE_MECHANISMS = sorted({
    "novel_uncatalogued_NA_active_site_substitution",
    "permissive_secondary_NA_mutations",      # framework-stabilizing background that enables a primary marker
    "PA_PB2_baloxavir_class_resistance",       # cap-endonuclease (I38X) — a DIFFERENT drug class, not NAI
    "M2_adamantane_resistance",                # S31N etc. — a different target/drug class entirely
    "reassortment_or_mixed_population_minor_variant",
})


def supported_antiviral_drugs() -> list[str]:
    return sorted(ANTIVIRAL_RESISTANCE_MUTATIONS)


def resistance_mutations_for(drug: str) -> dict[str, set[str]]:
    """gene -> resistance-substitution set for an antiviral drug. Raises KeyError on unknown drug."""
    d = drug.lower()
    if d not in ANTIVIRAL_RESISTANCE_MUTATIONS:
        raise KeyError(f"no antiviral resistance catalog for drug {drug!r}; "
                       f"configured: {supported_antiviral_drugs()}")
    return ANTIVIRAL_RESISTANCE_MUTATIONS[d]


def is_resistance_mutation(drug: str, gene: str, substitution: str) -> bool:
    """True iff <substitution> (e.g. 'H275Y') in <gene> is a recognized NAI-R marker for <drug>."""
    return substitution in ANTIVIRAL_RESISTANCE_MUTATIONS.get(drug.lower(), {}).get(gene, set())


@dataclass(frozen=True)
class AntiviralCall:
    """A deterministic antiviral R/S call (shape mirrors AntimalarialCall / FungalCall / amr_rules)."""
    prediction: str            # R / S / INDETERMINATE
    drug: str
    determinants: list[str]    # e.g. ["NA:H275Y"]
    undetectable_mechanisms: list[str]
    rule: str
    caveat: str


def call_from_observed_substitutions(drug: str, observed: dict[str, set[str]]) -> AntiviralCall:
    """Deterministic NAI R/S call from a virus's observed NA substitutions.

    `observed` = {gene: {substitutions present}} (from the BLAST caller). Rule (v0): R iff the NA carries
    >=1 recognized NAI-resistance marker for the drug. An S call surfaces ANTIVIRAL_UNDETECTABLE_MECHANISMS
    — it means 'no recognized NA marker found', NOT 'definitely inhibitor-sensitive'.
    """
    if drug.lower() not in ANTIVIRAL_RESISTANCE_MUTATIONS:
        return AntiviralCall("INDETERMINATE", drug, [], [], "antiviral_na_target_mutation_v0",
                             f"no catalog for {drug!r}")
    hits = []
    for gene, subs in resistance_mutations_for(drug).items():
        for s in observed.get(gene, set()):
            if s in subs:
                hits.append(f"{gene}:{s}")
    pred = "R" if hits else "S"
    undetectable = ANTIVIRAL_UNDETECTABLE_MECHANISMS if pred == "S" else []
    caveat = (f"deterministic NA target-site call ({ANTIVIRAL_DRUG_CLASS.get(drug.lower(), '?')}, N1 "
              "numbering). Hand-curated CDC/WHO-AVWG-recognized marker catalog (no AMRFinder-equivalent "
              "for influenza). NAI resistance is an NI-assay IC50 fold-change, not a growth MIC — the "
              "validated NA marker IS the genotypic call. "
              + ("An S call cannot rule out non-NA / different-drug-class resistance "
                 f"({', '.join(ANTIVIRAL_UNDETECTABLE_MECHANISMS)})." if pred == "S" else ""))
    return AntiviralCall(pred, drug, sorted(hits), undetectable, "antiviral_na_target_mutation_v0", caveat)
