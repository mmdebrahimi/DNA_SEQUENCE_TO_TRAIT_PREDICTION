"""MIC-tier classification + per-drug breakpoints (CLSI 2024 + EUCAST 14.0).

Shared module lifted from `scripts/cipro_mic_audit.py::_confidence_tier` +
`scripts/bvbrc_strict_mic_4drug_census.py::classify_strict_mic`. Both call
sites previously duplicated the tier-classification logic with cipro
breakpoints hardcoded vs parameterized; this module is the single source
of truth.

E. coli breakpoint catalog covers the Phase 2-relevant drugs:

  - Ciprofloxacin  (fluoroquinolone, concentrated-signal QRDR mechanism)
  - Ceftriaxone    (β-lactam, concentrated-signal plasmid β-lactamase mechanism)
  - Tetracycline   (distributed mobile-element efflux mechanism)
  - Gentamicin     (aminoglycoside, 4th-mechanism-class candidate)

Add drugs to `DRUG_BREAKPOINTS` to extend (not callers — callers should
discover the catalog via `breakpoints_for(drug)`).
"""
from __future__ import annotations

from math import isnan
from statistics import median
from typing import Optional


# ---------------------------------------------------------------------------
# Per-drug breakpoint catalog (CLSI 2024 + EUCAST 14.0, E. coli)
# ---------------------------------------------------------------------------
#
# Convention: breakpoints expressed as the lowest MIC integer at which a call
# becomes R or S. For example, "CLSI R>=2" means MIC values 2, 4, 8, 16, ...
# are CLSI-R; MIC 1 is intermediate.
#
# EUCAST values may be None when no E. coli breakpoints exist for the drug
# (e.g., tetracycline uses ECOFF only); CLSI is always present.

DRUG_BREAKPOINTS: dict[str, dict[str, Optional[float]]] = {
    "ciprofloxacin": {
        "clsi_r": 2.0, "clsi_s": 0.5,
        "eucast_r": 1.0, "eucast_s": 0.25,
    },
    "ceftriaxone": {
        "clsi_r": 4.0, "clsi_s": 1.0,
        "eucast_r": 2.0, "eucast_s": 1.0,
    },
    "tetracycline": {
        "clsi_r": 16.0, "clsi_s": 4.0,
        "eucast_r": None, "eucast_s": None,
    },
    "gentamicin": {
        "clsi_r": 16.0, "clsi_s": 4.0,
        "eucast_r": 4.0, "eucast_s": 2.0,
    },
    "meropenem": {
        # CLSI 2024 + EUCAST 14.0 Enterobacterales (carbapenem)
        "clsi_r": 4.0, "clsi_s": 1.0,
        "eucast_r": 8.0, "eucast_s": 2.0,
    },
}


class UnknownDrugError(KeyError):
    """Raised when a drug name is not in the breakpoint catalog."""


def breakpoints_for(drug: str) -> dict[str, Optional[float]]:
    """Return the breakpoint dict for a drug. Raises UnknownDrugError if absent."""
    bp = DRUG_BREAKPOINTS.get(drug.lower())
    if bp is None:
        raise UnknownDrugError(
            f"No breakpoints configured for drug: {drug!r}. "
            f"Configured drugs: {sorted(DRUG_BREAKPOINTS.keys())}"
        )
    return bp


def supported_drugs() -> list[str]:
    """List drug names with breakpoint entries."""
    return sorted(DRUG_BREAKPOINTS.keys())


# ---------------------------------------------------------------------------
# Tier classification
# ---------------------------------------------------------------------------


def classify_tier(
    mics: list[float],
    distinct_calls: set[str],
    breakpoints: dict[str, Optional[float]],
) -> str:
    """Classify a strain by MIC + label coherence.

    Returns one of:
      HIGH_R       — median MIC >= 4 * CLSI-R breakpoint (definite R, 4x safety margin)
      DECISIVE_R   — 2 * CLSI_R < median MIC < 4 * CLSI_R (clearly R, no safety margin)
      HIGH_S       — median MIC <= CLSI-S / 4 (definite S, 4x safety margin)
      DECISIVE_S   — CLSI_S / 4 < median MIC < CLSI_S / 2 (clearly S, no safety margin)
      BORDERLINE   — CLSI_S / 2 <= median MIC <= 2 * CLSI_R (gray zone)
      AMBIGUOUS    — CLSI vs EUCAST disagree on R/S call (skipped when EUCAST=None)
      CONFLICT     — multiple AST rows disagree on R/S (per-strain not per-row)
      NO_MIC       — no numeric MIC values

    Strict-MIC pass = HIGH_R or HIGH_S only.
    Relaxed-MIC pass = HIGH_R + DECISIVE_R + HIGH_S + DECISIVE_S.

    Compares against `proba` (or here, `med`) directly rather than `abs(...)`
    to keep float-precision behavior consistent with the cipro_mic_audit
    precursor (no boundary drift at MIC == 2 * CLSI_R, etc.).
    """
    r_calls = distinct_calls & {"R", "RESISTANT"}
    s_calls = distinct_calls & {"S", "SUSCEPTIBLE"}
    if r_calls and s_calls:
        return "CONFLICT"
    valid_mics = [m for m in mics if m is not None and not isnan(m)]
    if not valid_mics:
        return "NO_MIC"
    med = median(valid_mics)

    clsi_r = breakpoints["clsi_r"]
    clsi_s = breakpoints["clsi_s"]
    eucast_r = breakpoints.get("eucast_r")
    eucast_s = breakpoints.get("eucast_s")

    clsi_call = "R" if med >= clsi_r else ("S" if med <= clsi_s else "I")
    if eucast_r is not None and eucast_s is not None:
        eucast_call = "R" if med >= eucast_r else ("S" if med <= eucast_s else "I")
        if clsi_call != eucast_call:
            return "AMBIGUOUS"

    if med >= 4 * clsi_r:
        return "HIGH_R"
    if med <= clsi_s / 4:
        return "HIGH_S"
    if clsi_s / 2 <= med <= 2 * clsi_r:
        return "BORDERLINE"
    if med > 2 * clsi_r:
        return "DECISIVE_R"
    return "DECISIVE_S"


# ---------------------------------------------------------------------------
# Per-drug mechanism catalog (the "audit framework relaxes upstream filtering"
# insight from /review 2026-05-18: mechanism classes are drug-specific; the
# cipro_* audit scripts had this hardcoded. Pulling here so cef + tet + gent
# audits can share the same logic.)
# ---------------------------------------------------------------------------
#
# Per-drug AMRFinder Class filter (what gets kept from the mutations.tsv +
# main.tsv outputs). MULTIDRUG is included for drugs where regulatory hits
# show up under MULTIDRUG class (acrR, marR, etc.) — see CLAUDE.md gotcha
# "AMRFinder cipro-relevant Class filter: keep MULTIDRUG".

DRUG_AMRFINDER_CLASSES: dict[str, frozenset[str]] = {
    "ciprofloxacin": frozenset({"QUINOLONE", "FLUOROQUINOLONE", "MULTIDRUG"}),
    "ceftriaxone": frozenset({"BETA-LACTAM", "CARBAPENEM", "CEPHALOSPORIN", "MULTIDRUG"}),
    "tetracycline": frozenset({"TETRACYCLINE", "MULTIDRUG"}),
    "gentamicin": frozenset({"AMINOGLYCOSIDE", "MULTIDRUG"}),
    "meropenem": frozenset({"BETA-LACTAM", "CARBAPENEM", "MULTIDRUG"}),
}


def amrfinder_classes_for(drug: str) -> frozenset[str]:
    """Return the AMRFinder Class filter for a drug. Raises UnknownDrugError if absent."""
    cls = DRUG_AMRFINDER_CLASSES.get(drug.lower())
    if cls is None:
        raise UnknownDrugError(
            f"No AMRFinder class filter configured for drug: {drug!r}. "
            f"Configured drugs: {sorted(DRUG_AMRFINDER_CLASSES.keys())}"
        )
    return cls


# ---------------------------------------------------------------------------
# Per-drug mechanism catalog (loci → mechanism class mapping)
# ---------------------------------------------------------------------------
#
# Sourced primarily from CARD + AMRFinder + textbook reviews. Cipro catalog
# matches `scripts/cipro_mechanism_audit.py::CIPRO_LOCI_BY_MECHANISM`
# 1-for-1 (preserved as the canonical Phase 1 catalog).
#
# Cross-drug shared mechanisms (efflux pump components, porin loss,
# regulatory frameshifts) are tracked separately in CO_RESISTANCE_MECHANISMS
# below — these affect multiple drug classes and are conventionally classified
# as co-resistance modifiers rather than primary determinants.

DRUG_LOCI_BY_MECHANISM: dict[str, dict[str, set[str]]] = {
    "ciprofloxacin": {
        "QRDR_target_alteration": {"gyrA", "gyrB", "parC", "parE"},
        "plasmid_protect_modify": {
            "qnrA", "qnrB", "qnrC", "qnrD", "qnrS",
            "aac(6')-Ib-cr", "aac(6')-Ib", "aac6-Ib-cr",
        },
        "efflux": {"acrA", "acrB", "tolC", "oqxA", "oqxB", "mdfA", "mdtK"},
        "porin_loss": {"ompC", "ompF"},
        "regulatory": {"marR", "marA", "marB", "soxR", "soxS", "acrR"},
    },
    "ceftriaxone": {
        "acquired_beta_lactamase": {
            # ESBL + AmpC + carbapenemase families that hydrolyze 3rd-gen cephalosporins
            "blaCTX-M", "blaCMY", "blaTEM", "blaSHV", "blaOXA",
            "blaKPC", "blaNDM", "blaIMP", "blaVIM", "blaACT", "blaDHA", "blaMOX",
        },
        "ampC_hyperproduction": {"ampC"},
        "efflux": {"acrA", "acrB", "tolC", "mdtK"},
        "porin_loss": {"ompC", "ompF"},
        "regulatory": {"marR", "marA", "soxR", "soxS", "acrR", "ampR"},
    },
    "tetracycline": {
        "tet_efflux": {
            "tetA", "tetB", "tetC", "tetD", "tetE", "tetG", "tetH",
            "tetJ", "tetK", "tetL", "tetY", "tetZ", "tet39", "tet42",
        },
        "tet_ribosomal_protection": {
            "tetM", "tetO", "tetQ", "tetS", "tetT", "tetW", "tetBP",
        },
        "tet_enzymatic": {"tetX"},  # monooxygenase that inactivates tet
        "efflux": {"acrA", "acrB", "tolC"},
        "porin_loss": {"ompC", "ompF"},
        "regulatory": {"marR", "marA", "tetR", "acrR"},
    },
    "meropenem": {
        "carbapenemase_acquired": {
            # the carbapenemase families (AMRFinder Subclass CARBAPENEM) that hydrolyze meropenem
            "blaKPC", "blaNDM", "blaOXA-48", "blaOXA-181", "blaOXA-232",
            "blaVIM", "blaIMP", "blaGES", "blaSPM", "blaSME", "blaIMI",
        },
        # ESBL/AmpC + porin loss can raise meropenem MIC but are NOT carbapenemases — modifiers here
        "esbl_ampc_plus_porin": {"blaCTX-M", "blaCMY", "blaSHV", "ampC"},
        "efflux": {"acrA", "acrB", "tolC"},
        "porin_loss": {"ompC", "ompF", "ompK35", "ompK36"},
        "regulatory": {"marR", "marA", "soxR", "soxS", "acrR"},
    },
    "gentamicin": {
        "aminoglycoside_modifying_enzymes": {
            # Acetyltransferases (aac), phosphotransferases (aph),
            # nucleotidyltransferases (ant). The aac(6')-Ib-cr cross-listed
            # with cipro plasmid_protect_modify is included here too since
            # it's structurally a gent-acetylating enzyme.
            "aac(3)-IIa", "aac(3)-IId", "aac(3)-IV", "aac(3)-VIa",
            "aac(6')-Ib", "aac(6')-Iy", "aac(6')-Ib-cr",
            "ant(2'')-Ia", "ant(4')-Ia",
            "aph(2'')-Ia", "aph(3')-Ia", "aph(3')-IIa", "aph(3')-VI",
        },
        "16S_rRNA_methyltransferase": {
            "armA", "rmtA", "rmtB", "rmtC", "rmtD", "rmtE", "rmtF", "rmtG", "npmA",
        },
        "efflux": {"acrA", "acrB", "tolC"},
        "porin_loss": {"ompC", "ompF"},
        "regulatory": {"marR", "marA", "acrR"},
    },
}

# Per-drug primary mechanisms (vs co-resistance modifiers). A strain with a
# primary mechanism is the "textbook explanation" for resistance; efflux +
# regulatory + porin_loss are real but considered modifiers.
DRUG_PRIMARY_MECHANISMS: dict[str, frozenset[str]] = {
    "ciprofloxacin": frozenset({"QRDR_target_alteration", "plasmid_protect_modify"}),
    "ceftriaxone": frozenset({"acquired_beta_lactamase", "ampC_hyperproduction"}),
    "tetracycline": frozenset({"tet_efflux", "tet_ribosomal_protection", "tet_enzymatic"}),
    "gentamicin": frozenset({"aminoglycoside_modifying_enzymes", "16S_rRNA_methyltransferase"}),
    "meropenem": frozenset({"carbapenemase_acquired"}),
}

# Cross-drug shared co-resistance modifiers (the audit framework classifies
# these as "OPAQUE" when present without a primary mechanism)
CO_RESISTANCE_MECHANISMS: frozenset[str] = frozenset({"efflux", "regulatory", "porin_loss"})


def loci_by_mechanism_for(drug: str) -> dict[str, set[str]]:
    """Return the mechanism → loci dict for a drug. Raises UnknownDrugError if absent."""
    catalog = DRUG_LOCI_BY_MECHANISM.get(drug.lower())
    if catalog is None:
        raise UnknownDrugError(
            f"No loci-by-mechanism catalog configured for drug: {drug!r}. "
            f"Configured drugs: {sorted(DRUG_LOCI_BY_MECHANISM.keys())}"
        )
    return catalog


def primary_mechanisms_for(drug: str) -> frozenset[str]:
    """Return the primary-mechanism set for a drug. Raises UnknownDrugError if absent."""
    prim = DRUG_PRIMARY_MECHANISMS.get(drug.lower())
    if prim is None:
        raise UnknownDrugError(
            f"No primary mechanisms configured for drug: {drug!r}. "
            f"Configured drugs: {sorted(DRUG_PRIMARY_MECHANISMS.keys())}"
        )
    return prim


def classify_gene_symbol(drug: str, symbol: str) -> str:
    """Return the mechanism class for a gene symbol under the given drug.

    Returns "" when the symbol is not in the drug's catalog. Tolerant
    prefix-match (e.g., `qnrB19` → `qnrB`; `gyrA_S83L` → `gyrA`).

    Mirrors `scripts/cipro_mechanism_audit.py::_classify_symbol` but
    parameterized over the drug.
    """
    if not symbol:
        return ""
    catalog = loci_by_mechanism_for(drug)
    sym_norm = symbol.split("_")[0].strip()  # gyrA_S83L -> gyrA
    for mech, loci in catalog.items():
        if sym_norm in loci:
            return mech
        for locus in loci:
            if sym_norm.startswith(locus) and sym_norm != locus:
                return mech
    return ""


# ---------------------------------------------------------------------------
# Strict-vs-relaxed pass helpers (for callers building feasibility censuses
# or noise-class merges)
# ---------------------------------------------------------------------------


STRICT_MIC_TIERS: frozenset[str] = frozenset({"HIGH_R", "HIGH_S"})
RELAXED_MIC_TIERS: frozenset[str] = frozenset({"HIGH_R", "HIGH_S", "DECISIVE_R", "DECISIVE_S"})


def is_strict_r(tier: str) -> bool:
    return tier == "HIGH_R"


def is_strict_s(tier: str) -> bool:
    return tier == "HIGH_S"


def is_relaxed_r(tier: str) -> bool:
    return tier in ("HIGH_R", "DECISIVE_R")


def is_relaxed_s(tier: str) -> bool:
    return tier in ("HIGH_S", "DECISIVE_S")
