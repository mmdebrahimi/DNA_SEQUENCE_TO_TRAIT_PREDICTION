"""Fungal azole/echinocandin resistance determinant catalog (EP-7 — eukaryotic kingdom jump).

The first eukaryotic substrate (roadmap Phase 6) is FUNGAL AMR — Candida auris azole resistance — because
it is the kingdom jump that reuses the project's PROVEN deterministic-determinant method: same phenotype
class (AMR / MIC), documented target-site mutations, NO foundation model / GPU / money needed (a BLAST →
known-mutation scan, the eukaryotic analogue of `amr_rules.py` over AMRFinder).

KEY DIFFERENCE from bacterial AMR: there is NO AMRFinder-equivalent for fungi ("genome-based predictors of
AMR are not available for fungal pathogens" — Frontiers Microbiol 2023, PMC10157239). So the determinant
catalog must be HAND-CURATED (this module) and the calling step is custom (BLAST the gene allele vs the
genome → translate → check the codon against the resistance-substitution set). BLAST+ is installed natively
(`C:/Users/Farshad/ncbi-blast/bin`, per the pathotype work).

This module is the substrate-independent core (domain knowledge — valid regardless of cohort). The
BLAST-based caller + the WGS+MIC validation cohort (from paper supplementaries: 188 S. Africa /
350 India, near-perfect ERG11↔MIC linkage) are EP-7 steps 2-3 (see plans/EP7_Fungal_AMR_Substrate.md).

Determinants sourced from the 2026-06-07 substrate survey:
- Lockhart et al. (clade-specific ERG11): Y132F (Venezuela/India/Pakistan), K143R (India/Pakistan),
  F126T (South Africa) — `research_outputs/eukaryotic-multimodal-substrate-feasibility-2026-06-07.md`.
- S. Africa outbreak (PMC10521600): clade III VF125AL; clade IV K177R/N335S/E343D.
- Multi-locus caveat (PMC8092288): ERG11 alone insufficient — TAC1b + Cdr1 efflux + aneuploidy also drive
  azole-R (the fungal analogue of the bacterial efflux blind spot → `undetectable_mechanisms`).
"""
from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Resistance-determinant catalog: drug -> gene -> set of resistance amino-acid substitutions.
# Substitution form: "<wt><pos><mut>" (1-letter, 1-based), e.g. "Y132F". Detected by translating the
# BLAST-aligned gene and checking the residue at <pos>.
# ---------------------------------------------------------------------------

FUNGAL_RESISTANCE_MUTATIONS: dict[str, dict[str, set[str]]] = {
    # Azoles (fluconazole / voriconazole): target = ERG11 (lanosterol 14-alpha-demethylase, Cyp51A).
    "fluconazole": {
        "ERG11": {
            # C. auris clade-specific azole-R substitutions (Lockhart 2017 + S. Africa 2023)
            "Y132F", "K143R", "F126T", "F126L", "VF125AL",
            "K177R", "N335S", "E343D",
            # C. albicans canonical ERG11 azole-R hotspots (cross-species, well documented)
            "Y132H", "G464S", "S405F", "G448E", "F145L", "D446E",
        },
    },
    "voriconazole": {
        # voriconazole shares the ERG11 target; same substitution set drives cross-azole resistance
        "ERG11": {
            "Y132F", "K143R", "F126T", "F126L", "Y132H", "G464S", "G448E",
        },
    },
    # Echinocandins (caspofungin / micafungin / anidulafungin): target = FKS1 (beta-1,3-glucan synthase),
    # hot-spot region HS1 (~S639 in C. auris numbering).
    "caspofungin": {
        "FKS1": {"S639F", "S639P", "S639Y", "F635del", "R1354H"},
    },
    "micafungin": {
        "FKS1": {"S639F", "S639P", "S639Y", "F635del", "R1354H"},
    },
}

# Lineage / clade-background markers: catalogued substitutions that CO-SEGREGATE with a resistant clade
# but do NOT independently discriminate R from S. An R call driven ONLY by these is LOW-confidence
# ("lineage-associated, mechanism-unconfirmed") rather than a mechanism-attributable R.
# C. auris clade IV ERG11 haplotype K177R/N335S/E343D: the AR-Bank re-validation (2026-07-20) found the
# IDENTICAL haplotype in a fluconazole-R isolate (SAMN11570381) AND a fluconazole-S isolate (SAMN10139552)
# -- same genotype, opposite phenotype -> +1 TP and +1 FP -> zero net discrimination. The originating
# S. Africa outbreak paper (PMC10521600) reported these in clade IV R isolates but had no clade IV S
# control; SAMN10139552 is that control. They stay in the catalog (preserve sensitivity -- clade IV is
# often genuinely R) but a haplotype-only call is flagged LOW confidence. See
# wiki/ar_bank_caur_powered_validation_*.  The fungal analogue of the documented QRDR-vs-lineage confound.
FUNGAL_LINEAGE_MARKERS: dict[str, dict[str, set[str]]] = {
    "fluconazole": {"ERG11": {"K177R", "N335S", "E343D"}},
}


def lineage_markers_for(drug: str) -> dict[str, set[str]]:
    """gene -> lineage/clade-background substitution set (non-discriminative) for a fungal drug."""
    return FUNGAL_LINEAGE_MARKERS.get(drug.lower(), {})


# Mechanisms NOT detectable by a target-site-mutation scan (the fungal analogue of the bacterial
# efflux/porin/regulatory blind spot). A SUSCEPTIBLE call cannot rule these out.
FUNGAL_UNDETECTABLE_MECHANISMS = sorted({
    "TAC1b_efflux_regulator",   # upregulates CDR1/CDR2 efflux (azole-R without ERG11 mutation)
    "CDR1_efflux_overexpression",
    "MDR1_efflux_overexpression",
    "ERG11_copy_number_aneuploidy",   # chr5/ERG11 duplication raises target dosage
    "ERG3_loss_of_function",
})

# Drug class for reporting
FUNGAL_DRUG_CLASS = {
    "fluconazole": "AZOLE", "voriconazole": "AZOLE",
    "caspofungin": "ECHINOCANDIN", "micafungin": "ECHINOCANDIN",
}

# CDC *tentative* MIC breakpoints for C. auris resistance (ug/mL). There is NO formal CLSI/EUCAST
# species-specific breakpoint for C. auris; CDC publishes tentative ones widely used in genomic-epi.
# fluconazole >=32 = R is the canonical azole cutoff (e.g. the S.Africa outbreak: 181/188 isolates
# with MIC>32 carried ERG11 mutations). Source: CDC C. auris antifungal-susceptibility guidance.
CAURIS_TENTATIVE_R_MIC = {
    "fluconazole": 32.0,
    "caspofungin": 2.0,
    "micafungin": 4.0,
}


def mic_to_phenotype(drug: str, mic: float) -> str | None:
    """Map an MIC (ug/mL) to 'R'/'S' via the CDC tentative C. auris breakpoint.

    Returns None when no tentative breakpoint is configured for the drug (caller treats as unlabelable).
    R iff mic >= breakpoint (CDC uses >= for the resistant call).
    """
    bp = CAURIS_TENTATIVE_R_MIC.get(drug.lower())
    if bp is None:
        return None
    return "R" if float(mic) >= bp else "S"


def supported_fungal_drugs() -> list[str]:
    return sorted(FUNGAL_RESISTANCE_MUTATIONS)


def resistance_mutations_for(drug: str) -> dict[str, set[str]]:
    """gene -> resistance-substitution set for a fungal drug. Raises KeyError on unknown drug."""
    d = drug.lower()
    if d not in FUNGAL_RESISTANCE_MUTATIONS:
        raise KeyError(f"no fungal resistance catalog for drug {drug!r}; "
                       f"configured: {supported_fungal_drugs()}")
    return FUNGAL_RESISTANCE_MUTATIONS[d]


def is_resistance_mutation(drug: str, gene: str, substitution: str) -> bool:
    """True iff <substitution> (e.g. 'Y132F') in <gene> is a known resistance determinant for <drug>."""
    return substitution in FUNGAL_RESISTANCE_MUTATIONS.get(drug.lower(), {}).get(gene, set())


@dataclass(frozen=True)
class FungalCall:
    """A deterministic fungal R/S call (shape mirrors amr_rules.call_resistance output)."""
    prediction: str            # R / S / INDETERMINATE
    drug: str
    determinants: list[str]    # e.g. ["ERG11:Y132F"]
    undetectable_mechanisms: list[str]
    rule: str
    caveat: str
    confidence: str = "NA"     # HIGH (causal marker) / LOW_LINEAGE_ONLY (clade-background only) / NA (S/indet)
    lineage_only_determinants: tuple[str, ...] = ()   # the clade-background hits that drove a LOW call


def call_from_observed_substitutions(drug: str, observed: dict[str, set[str]]) -> FungalCall:
    """Deterministic azole/echinocandin R/S call from a genome's observed target-gene substitutions.

    `observed` = {gene: {substitutions present in this genome}} (produced by the EP-7 BLAST caller — step 2).
    Rule (v0): R iff the genome carries >=1 catalogued resistance substitution in a target gene. Mirrors
    the acquired-determinant bacterial rule. An S call surfaces FUNGAL_UNDETECTABLE_MECHANISMS (efflux /
    aneuploidy) — it means 'no target-site resistance mutation found', NOT 'definitely susceptible'.
    """
    if drug.lower() not in FUNGAL_RESISTANCE_MUTATIONS:
        return FungalCall("INDETERMINATE", drug, [], [], "fungal_target_mutation_v0",
                          f"no catalog for {drug!r}")
    lineage = lineage_markers_for(drug)
    hits, causal_hits, lineage_hits = [], [], []
    for gene, subs in resistance_mutations_for(drug).items():
        lin_g = lineage.get(gene, set())
        for s in observed.get(gene, set()):
            if s in subs:
                tag = f"{gene}:{s}"
                hits.append(tag)
                (lineage_hits if s in lin_g else causal_hits).append(tag)
    pred = "R" if hits else "S"
    # Confidence: causal marker present -> HIGH; R driven ONLY by clade-background haplotype -> LOW.
    if pred == "S":
        confidence = "NA"
    elif causal_hits:
        confidence = "HIGH"
    else:
        confidence = "LOW_LINEAGE_ONLY"
    undetectable = FUNGAL_UNDETECTABLE_MECHANISMS if pred == "S" else []
    caveat = (f"deterministic target-site call ({FUNGAL_DRUG_CLASS.get(drug.lower(), '?')}). "
              "Hand-curated catalog (no AMRFinder-equivalent for fungi). "
              + ("An S call cannot rule out efflux/aneuploidy-mediated resistance "
                 f"({', '.join(FUNGAL_UNDETECTABLE_MECHANISMS)})." if pred == "S" else "")
              + (f"R driven only by clade-background marker(s) {sorted(lineage_hits)} "
                 "(lineage-associated, mechanism-unconfirmed; non-discriminative within the clade)."
                 if confidence == "LOW_LINEAGE_ONLY" else ""))
    return FungalCall(pred, drug, sorted(hits), undetectable, "fungal_target_mutation_v0", caveat,
                      confidence=confidence, lineage_only_determinants=tuple(sorted(lineage_hits)))
