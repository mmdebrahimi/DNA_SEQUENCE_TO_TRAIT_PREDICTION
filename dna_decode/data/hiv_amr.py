"""HIV-1 NNRTI resistance determinant catalog — a SECOND viral target-site decoder.

After influenza A neuraminidase (the first viral cell), this is a second viral target within
the viral kingdom: *HIV-1* non-nucleoside reverse-transcriptase-inhibitor (NNRTI) resistance via
the established major NNRTI drug-resistance mutations (DRMs) on **reverse transcriptase (RT)**.

WHY THIS ONE IS DIFFERENT (the reason it exists): it is the FIRST target-site cell that can be
VALIDATED against a FREE, isolate-level, INDEPENDENT laboratory genotype-phenotype dataset — the
Stanford HIVDB genotype-phenotype dataset (PhenoSense fold-decreased susceptibility, Monogram). The
phenotype is a wet-lab IC50/fold-change INDEPENDENT of any genotype-interpretation rule, so scoring
this catalog against it is a genuine genotype->independent-phenotype test (it clears the project's
circular-label gate). Validation MUST use the PhenoSense fold-change, NEVER HIVDB's own Sierra/GRT-IS
interpretation (rule-vs-rule = circular); a naive-Sierra / Stanford-R-script baseline is reported
alongside, headlining the delta (the "validate the wrapper vs the underlying tool" discipline).

PROVENANCE (no fabrication, no build-on-unverified): the NNRTI major DR positions are taken VERBATIM
from the Stanford HIVDB genotype-phenotype dataset page's "NNRTI Major Drug Resistance Positions"
definition (the exact positions used to build the high-quality-filtered dataset; page last updated
2026-04-22): 100I, 101P, 103N, 106A/M, 181C/I/A, 188C/L/H, 190A/E/S/Q, 230L. The wild-type residue at
each position is the consensus-B amino acid (the reference the page links; settled). This independently
CONFIRMS the consensus-B wild-type at RT-188 is Tyrosine (Y188C/L/H) — NOT the "G188" a web summarizer
garbled. Cross-reference: Bennett et al. 2009, PLoS ONE 4(3):e4724 (WHO surveillance DRM list). Required
HIVDB citation: Rhee et al. 2003, Nucleic Acids Res 31:298-303.

NUMBERING: HIV-1 RT amino-acid positions, consensus-B numbering. Mutation form "<wt><pos><mut>"
(1-letter, e.g. K103N) — the same shorthand HIVDB uses (the dataset's CompMutList column).

SCOPE (v0 — honest): CLASS-LEVEL. v0 treats "any established major NNRTI DRM present" as a reduced-NNRTI-
susceptibility signal and maps EVERY NNRTI drug to the same RT major-DRM set, because the Stanford source
defines the positions at the NNRTI-CLASS level. This is deliberately a first cut: per-drug differential
resistance is real and large (e.g. K103N strongly affects efavirenz/nevirapine but SPARES etravirine/
rilpivirine; Y181C affects the second-generation NNRTIs), so v0 will OVER-CALL resistance for the drugs a
mutation spares. The v0 validation reports per-drug performance precisely (expected: strong on EFV/NVP,
weaker on ETR/RPV) — that honest gap is what motivates a v0.1 per-drug catalog (which needs the per-drug
DRM associations sourced separately). NRTI/PI/INI/CAI targets + minor/accessory NNRTI mutations are out of
v0 scope.
"""
from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# NNRTI major DRMs on RT (consensus-B numbering), VERBATIM from the Stanford HIVDB
# genotype-phenotype dataset page's "NNRTI Major Drug Resistance Positions" (2026-04-22),
# expanded to <wt><pos><mut> using the consensus-B wild-type residue at each position.
#   100I        -> L100I            (wt L)
#   101P        -> K101P            (wt K)
#   103N        -> K103N            (wt K)
#   106A/M      -> V106A, V106M     (wt V)
#   181C/I/A    -> Y181C, Y181I, Y181A   (wt Y)
#   188C/L/H    -> Y188C, Y188L, Y188H   (wt Y  <- confirms NOT "G188")
#   190A/E/S/Q  -> G190A, G190E, G190S, G190Q  (wt G)
#   230L        -> M230L            (wt M)
# ---------------------------------------------------------------------------
NNRTI_RT_MAJOR_DRMS: set[str] = {
    "L100I",
    "K101P",
    "K103N",
    "V106A", "V106M",
    "Y181C", "Y181I", "Y181A",
    "Y188C", "Y188L", "Y188H",
    "G190A", "G190E", "G190S", "G190Q",
    "M230L",
}

# Consensus-B wild-type residue at each catalogued RT position (documentation + a self-check).
_RT_WT = {100: "L", 101: "K", 103: "K", 106: "V", 181: "Y", 188: "Y", 190: "G", 230: "M"}

# The NNRTI drugs. v0 = CLASS-LEVEL: each maps to the SAME RT major-DRM set (see module docstring
# SCOPE note — per-drug differential resistance is v0.1). delavirdine (DLV) is the old first-gen NNRTI;
# doravirine (DOR) the newest. The dataset carries a per-drug fold-change column for whichever it covers.
HIV_NNRTI_DRUGS = ("nevirapine", "efavirenz", "etravirine", "rilpivirine", "doravirine", "delavirdine")

HIV_RESISTANCE_MUTATIONS: dict[str, dict[str, set[str]]] = {
    drug: {"RT": set(NNRTI_RT_MAJOR_DRMS)} for drug in HIV_NNRTI_DRUGS
}

HIV_DRUG_CLASS = {drug: "NNRTI" for drug in HIV_NNRTI_DRUGS}

# Mechanisms a major-NNRTI-DRM RT scan CANNOT see — a SUSCEPTIBLE (no-major-DRM) call cannot rule these
# out. This is the HIV analogue of the bacterial efflux / fungal aneuploidy / influenza different-class
# blind spot.
HIV_UNDETECTABLE_MECHANISMS = sorted({
    "minor_or_accessory_NNRTI_mutations",        # e.g. E138K/A/G/Q, V179D/F, K101E, V108I, H221Y, P225H
    "novel_or_uncatalogued_NNRTI_substitution",
    "per_drug_differential_resistance",          # v0 is class-level (K103N spares ETR/RPV) -> v0.1 caveat
    "NRTI_PI_INSTI_CAI_resistance",              # different target/drug classes entirely
    "non_subtypeB_specific_resistance_pathway",  # subtype-specific pathways outside the consensus-B catalog
    "mixture_or_minority_variant_below_sanger",
})

# Provenance carried into the caveat (immutable).
_SOURCE = ("Stanford HIVDB genotype-phenotype dataset 'NNRTI Major Drug Resistance Positions' "
           "(2026-04-22); consensus-B wild-types; xref Bennett 2009 PLoS ONE 4(3):e4724; "
           "cite Rhee 2003 Nucleic Acids Res 31:298-303")


def gene_for_drug(drug: str) -> str | None:
    """Single target gene per drug (RT for every NNRTI in v0). Routes the caller + picks the reference."""
    genes = HIV_RESISTANCE_MUTATIONS.get(drug.lower())
    return next(iter(genes), None) if genes else None


def supported_hiv_drugs() -> list[str]:
    return sorted(HIV_RESISTANCE_MUTATIONS)


def resistance_mutations_for(drug: str) -> dict[str, set[str]]:
    """gene -> major-DRM set for an HIV NNRTI drug. Raises KeyError on unknown drug."""
    d = drug.lower()
    if d not in HIV_RESISTANCE_MUTATIONS:
        raise KeyError(f"no HIV resistance catalog for drug {drug!r}; "
                       f"configured: {supported_hiv_drugs()}")
    return HIV_RESISTANCE_MUTATIONS[d]


def is_resistance_mutation(drug: str, gene: str, substitution: str) -> bool:
    """True iff <substitution> (e.g. 'K103N') in <gene> is a major NNRTI DRM for <drug>."""
    return substitution in HIV_RESISTANCE_MUTATIONS.get(drug.lower(), {}).get(gene, set())


@dataclass(frozen=True)
class HIVCall:
    """A deterministic HIV NNRTI R/S call (shape mirrors AntiviralCall / FungalCall / amr_rules)."""
    prediction: str            # R / S / INDETERMINATE
    drug: str
    determinants: list[str]    # e.g. ["RT:K103N"]
    undetectable_mechanisms: list[str]
    rule: str
    caveat: str


def call_from_observed_substitutions(drug: str, observed: dict[str, set[str]]) -> HIVCall:
    """Deterministic NNRTI R/S call from a virus's observed RT substitutions.

    `observed` = {gene: {substitutions present}} (e.g. parsed from the HIVDB dataset's CompMutList, or
    a BLAST caller). Rule (v0, class-level): R iff RT carries >=1 established major NNRTI DRM. An S call
    surfaces HIV_UNDETECTABLE_MECHANISMS — it means 'no major NNRTI DRM found', NOT 'definitely NNRTI-
    susceptible' (minor mutations, per-drug differences, and non-RT mechanisms are invisible to v0).
    """
    if drug.lower() not in HIV_RESISTANCE_MUTATIONS:
        return HIVCall("INDETERMINATE", drug, [], [], "hiv_nnrti_major_drm_v0",
                       f"no catalog for {drug!r}")
    hits = []
    for gene, subs in resistance_mutations_for(drug).items():
        for s in observed.get(gene, set()):
            if s in subs:
                hits.append(f"{gene}:{s}")
    pred = "R" if hits else "S"
    undetectable = HIV_UNDETECTABLE_MECHANISMS if pred == "S" else []
    caveat = (f"deterministic RT major-NNRTI-DRM call (NNRTI, consensus-B numbering; CLASS-LEVEL v0 — "
              f"per-drug differential resistance not yet modelled). Source: {_SOURCE}. Validate against "
              f"PhenoSense fold-change, NEVER HIVDB's own Sierra interpretation (circular). "
              + ("An S call cannot rule out minor/accessory NNRTI mutations, per-drug differences, or "
                 "non-RT / different-class resistance." if pred == "S" else ""))
    return HIVCall(pred, drug, sorted(hits), undetectable, "hiv_nnrti_major_drm_v0", caveat)


# ===========================================================================
# NRTI (v0 — POSITION-BASED; the second HIV drug class)
# ===========================================================================
# Stanford HIVDB genotype-phenotype dataset "NRTI Major Drug Resistance Positions" (2026-04-22):
# 41, 65, 70, 74, 75, 151, 184, 210, 215 — POSITIONS ONLY (unlike NNRTI, Stanford did NOT publish
# mutant-level NRTI majors on that page). So v0 NRTI is POSITION-BASED: any non-consensus mutation at one
# of these RT positions counts as a major NRTI DRM present. HONEST CAVEAT (vs the NNRTI mutant-level v0):
# this DELIBERATELY OVER-CALLS — T215 revertants (T215S/C/D/E/I/V) + V75 polymorphisms are non-consensus
# at a major position but are NOT resistant -> reduced specificity (esp. AZT/D4T). The validation
# quantifies the per-drug spec hit. v0.1 = a MUTANT-SPECIFIC catalog (data-derived from the OLS
# coefficients, or a sourced SDRM list). Consensus-B wild-types per the HIVDB-linked reference.
NRTI_MAJOR_POSITIONS = (41, 65, 70, 74, 75, 151, 184, 210, 215)
NRTI_RT_WT = {41: "M", 65: "K", 70: "K", 74: "L", 75: "V", 151: "Q", 184: "M", 210: "L", 215: "T"}
NRTI_DRUGS = ("lamivudine", "abacavir", "zidovudine", "stavudine", "didanosine", "tenofovir")

NRTI_UNDETECTABLE_MECHANISMS = sorted({
    "minor_or_accessory_NRTI_mutations",
    "per_drug_differential_resistance",          # position-based v0 over-calls (215 revertants etc.)
    "NNRTI_PI_INSTI_CAI_resistance",
    "non_subtypeB_specific_resistance_pathway",
    "mixture_or_minority_variant_below_sanger",
})

_NRTI_SOURCE = ("Stanford HIVDB genotype-phenotype dataset 'NRTI Major Drug Resistance Positions' "
                "(2026-04-22; positions only); consensus-B wild-types; cite Rhee 2003")


def supported_nrti_drugs() -> list[str]:
    return sorted(NRTI_DRUGS)


def is_nrti_major_position_mutation(substitution: str) -> bool:
    """True iff `substitution` ('<wt><pos><mut>', e.g. 'M184V') is at an NRTI major position."""
    digits = "".join(c for c in substitution[1:] if c.isdigit())
    try:
        return int(digits) in NRTI_MAJOR_POSITIONS
    except ValueError:
        return False


def call_nrti_from_observed(drug: str, observed: dict[str, set[str]]) -> HIVCall:
    """Deterministic POSITION-BASED NRTI R/S call from observed RT substitutions.

    Rule (v0): R iff RT carries >=1 non-consensus substitution at an NRTI major position (Stanford's
    filtered-dataset definition). DELIBERATELY over-calls (revertants/polymorphisms at a major position) —
    an honest coarse first cut whose specificity hit the validation quantifies. An S call surfaces
    NRTI_UNDETECTABLE_MECHANISMS (NOT 'definitely NRTI-susceptible')."""
    if drug.lower() not in NRTI_DRUGS:
        return HIVCall("INDETERMINATE", drug, [], [], "hiv_nrti_major_position_v0",
                       f"no NRTI catalog for {drug!r}")
    hits = sorted(f"RT:{s}" for s in observed.get("RT", set()) if is_nrti_major_position_mutation(s))
    pred = "R" if hits else "S"
    undetectable = NRTI_UNDETECTABLE_MECHANISMS if pred == "S" else []
    caveat = (f"deterministic RT major-NRTI-POSITION call (NRTI, consensus-B numbering; POSITION-BASED v0 "
              f"— over-calls 215 revertants / 75 polymorphisms; mutant-specific catalog is v0.1). Source: "
              f"{_NRTI_SOURCE}. Validate against PhenoSense fold-change, NEVER HIVDB's own Sierra "
              f"interpretation (circular). "
              + ("An S call cannot rule out minor NRTI mutations, per-drug differences, or non-RT / "
                 "different-class resistance." if pred == "S" else ""))
    return HIVCall(pred, drug, hits, undetectable, "hiv_nrti_major_position_v0", caveat)
