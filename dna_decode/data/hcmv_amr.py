"""Human cytomegalovirus (HCMV) antiviral-resistance catalog - the 4th validated-against-free-measured-label
viral cell (after HIV-1, SARS-CoV-2 Mpro, influenza NA), and the FIRST herpesvirus cell.

WHY THIS ONE (the free-independent-label breakthrough, extended to a new pathogen family): HCMV resistance is
the field's GOLD-STANDARD recombinant-phenotyping system - a resistance mutation is confirmed by transferring
it into a lab strain (marker transfer / BAC mutagenesis) and MEASURING the fold-increase in EC50 (ganciclovir/
cidofovir/foscarnet) or nM EC50 (letermovir). Those measured fold-changes are published in FREE open-access
compilations (Chou et al.), and they are wet-lab measurements INDEPENDENT of any genotype-interpretation rule
- so scoring a curated target-site catalog against them clears the project's circular-label gate, exactly as
the HIV (Stanford HIVDB PhenoSense) and SARS-CoV-2 (CoV-RDB) cells do. This mirrors those cells' architecture:
a curated target-site mutation catalog + a phenotyped-BENIGN class (the specificity anchor, like SARS-CoV-2
Mpro P132H).

TARGETS (v0 - GCV/CDV/FOS/letermovir):
  - UL97 (phosphotransferase) -> ganciclovir (GCV is a prodrug UL97 activates; a UL97 loss-of-activation
    mutation confers GCV resistance). ~5-10x EC50 canonical band.
  - UL54 (DNA polymerase) -> ganciclovir / cidofovir / foscarnet, with CROSS-RESISTANCE by pol domain
    (regions IV/V -> GCV+CDV; regions II/VI/delta-C -> FOS). Per-mutation drug set is catalogued, not assumed.
  - UL56 (terminase) -> letermovir (LMV). C325W/Y = high-grade (>8000x); R369/V236 = mid-grade.

PROVENANCE (no fabrication - every entry is tied to a fetched open-access source; unverifiable positions were
EXCLUDED, not invented):
  - PMC3262590 (Chou review) - UL54 fold-change Table 1 + the recombinant-phenotyped no-effect BENIGN list
    (Table 2) + the canonical UL97 GCV set. **Primary provenance for the catalog.**
  - PMC5483911 (Chou 2017 JCM) - UL97 codon 591-603 point + deletion fold-changes (dual cell-line HFF/ARPE).
  - AAC 10.1128/aac.00922-18 (Chou 2018) - the definitive UL56 letermovir recombinant EC50 table.
  - PMC9759347 - UL56 letermovir compilation + novel recombinants + its benign class.
NUMBERING: UL97/UL54/UL56 codon numbering referenced to lab strain AD169 (Chou-lab recombinant-phenotyping
convention). Mutation form "<wt><pos><mut>" (e.g. M460V), same shorthand as the HIV/SARS-CoV-2 cells;
in-frame deletions keep the source's form (e.g. "del591-594", "595del") as opaque catalog keys.

EXCLUDED as unverified (real literature positions, but no fold-change tied to a fetched table - not
fabricated): UL97 A594T / L595F / C607F / C607Y; UL56 E237G. Maribavir (UL27/UL97 F342Y...) is a separate
future cell (different drug; out of the GCV/CDV/FOS/LMV scope).

SCOPE (v0 - honest): MUTANT-LEVEL (not position-based) because HCMV genes carry benign polymorphisms - a
position rule would over-call. Wheel-only `--observed UL97:M460V,UL54:F412L` (genome-FASTA mode = v0.1, needs
committed UL97/UL54/UL56 CDS references + a BLAST caller, mirroring how the fungal/SARS-CoV-2 cells started).
Validation is IN-DISTRIBUTION against the SAME measured fold-changes the catalog is curated from (a knowledge
baseline, like the SARS-CoV-2 CoV-RDB cell) - an INDEPENDENT number needs held-out recombinant-phenotyping
studies (v0.1). An S call means "no catalogued resistance mutation in the drug's target gene(s)", NOT
"definitely susceptible".
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# --- drug -> target gene(s) routing (a drug is scored ONLY against its mechanism genes) -------------------
# valganciclovir is the oral prodrug of ganciclovir -> identical UL97/UL54 mechanism.
GENES_FOR_DRUG: dict[str, tuple[str, ...]] = {
    "ganciclovir": ("UL97", "UL54"),
    "valganciclovir": ("UL97", "UL54"),
    "cidofovir": ("UL54",),
    "foscarnet": ("UL54",),
    "letermovir": ("UL56",),
}

# --- resistance catalog: GENE -> { mutation : frozenset(drugs it confers resistance to) } ----------------
# Sourced VERBATIM (see module docstring). UL97 = all ganciclovir; UL54 = per-mutation drug set from the
# fold-change table (cross-resistance is catalogued, never assumed); UL56 = all letermovir.
_G = "ganciclovir"
_C = "cidofovir"
_F = "foscarnet"
_L = "letermovir"

# valganciclovir is the oral prodrug of ganciclovir -> IDENTICAL UL97/UL54 mechanism. Normalize it to
# ganciclovir for the catalog membership check (the catalog stores drug sets under 'ganciclovir'); the
# ORIGINAL drug name is preserved in the call/record.
_MECH_ALIAS: dict[str, str] = {"valganciclovir": "ganciclovir"}

UL97_RESISTANCE: dict[str, frozenset[str]] = {m: frozenset({_G}) for m in (
    # canonical GCV markers (PMC3262590 / PMC3773841)
    "M460V", "M460I", "H520Q", "C592G", "A594V", "L595S", "C603W",
    # codon 591-603 point + in-frame deletions with measured fold-change (PMC5483911 Table 1)
    "A591V", "del591-594", "595del", "595del9", "596del", "597del2", "597del3",
    "599del", "600del", "600del2", "601del", "601del2", "601del3",
)}

UL54_RESISTANCE: dict[str, frozenset[str]] = {
    # PMC3262590 Table 1 (fold-increase EC50; per-mutation drug set)
    "F412L": frozenset({_G, _C}), "F412S": frozenset({_G, _C}), "P488R": frozenset({_G, _C}),
    "K500N": frozenset({_G, _C}), "C539R": frozenset({_G, _C}), "L545W": frozenset({_G, _C}),
    "T552N": frozenset({_G, _F}), "Q578H": frozenset({_G, _C, _F}), "S585A": frozenset({_F}),
    "F595I": frozenset({_G}), "L802V": frozenset({_F}), "P829S": frozenset({_G}),
    "L862F": frozenset({_G}), "V946L": frozenset({_G}), "L957F": frozenset({_F}),
    # triple GCV+CDV+FOS cross-resistant (PMC3773841 text - qualitative multi-drug)
    "V812L": frozenset({_G, _C, _F}), "A834P": frozenset({_G, _C, _F}),
    "del981-982": frozenset({_G, _C, _F}),
}

UL56_RESISTANCE: dict[str, frozenset[str]] = {m: frozenset({_L}) for m in (
    # AAC 10.1128/aac.00922-18 Table 1 + PMC9759347
    "V231L", "V236M", "V236A", "L257F", "C325Y", "C325W", "C325F", "C325R",
    "A365S", "L328V", "C25F", "R369S", "R369G", "R369T", "R369K", "Q234R", "V363I",
)}

RESISTANCE_BY_GENE: dict[str, dict[str, frozenset[str]]] = {
    "UL97": UL97_RESISTANCE, "UL54": UL54_RESISTANCE, "UL56": UL56_RESISTANCE,
}

# --- BENIGN class (recombinant-phenotyped, no significant fold-change) - the specificity anchor -----------
# A benign polymorphism must return S, never R (like SARS-CoV-2 Mpro P132H). These are NOT in the resistance
# catalog, so the call naturally returns S; the sets are kept for the specificity tests + provenance.
UL97_BENIGN: frozenset[str] = frozenset({"K599E", "T601M"})                          # PMC5483911 Table 1
UL54_BENIGN: frozenset[str] = frozenset({                                            # PMC3262590 Table 2
    "S291P", "C304S", "E315D", "F357L", "V377A", "K415R", "P448S", "R512H", "C524R", "E530G",
    "V544A", "D576G", "A635T", "P648S", "S649P", "A692V", "V694I", "S695T", "A714V", "N757K",
    "I837V", "D879G", "S897P", "L926V", "A972V", "G993D",
})
UL56_BENIGN: frozenset[str] = frozenset({                                            # PMC9759347 Table 3
    "L243P", "S262C", "L328I", "H335Y", "E339G", "K350R", "N368I", "T399I",
})
BENIGN_BY_GENE: dict[str, frozenset[str]] = {"UL97": UL97_BENIGN, "UL54": UL54_BENIGN, "UL56": UL56_BENIGN}

HCMV_UNDETECTABLE_MECHANISMS: list[str] = [
    "minor_or_novel_substitution_below_catalog",
    "in_frame_deletion_not_in_v0_point_catalog",
    "non_target_gene_or_epistatic_pathway",
    "per_drug_differential_within_cross_resistance_domain",
    "maribavir_UL27_UL97_resistance",   # a different drug/target - separate future cell
]

_SOURCE = ("Chou et al. recombinant-phenotyping compilations - PMC3262590 (UL54 fold-change + benign) / "
           "PMC5483911 (UL97 591-603) / AAC 10.1128/aac.00922-18 + PMC9759347 (UL56 letermovir); "
           "AD169 codon numbering")

_POINT = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY*])$")


@dataclass(frozen=True)
class HCMVCall:
    """A deterministic HCMV R/S call (shape mirrors HIVCall / SARSCoV2Call - duck-typed by the CLI record
    builder `_target_site_record`)."""
    prediction: str            # R / S / INDETERMINATE
    drug: str
    determinants: list[str]    # e.g. ["UL97:M460V", "UL54:Q578H"]
    undetectable_mechanisms: list[str]
    rule: str
    caveat: str


def genes_for_hcmv_drug(drug: str) -> tuple[str, ...]:
    """Target gene(s) for an HCMV drug - routes the caller + (v0.1) the reference. () for an unknown drug."""
    return GENES_FOR_DRUG.get(drug.lower(), ())


def all_supported_hcmv_drugs() -> list[str]:
    """Every HCMV drug the decoder routes (v0 = GCV/valGCV/CDV/FOS/letermovir)."""
    return sorted(GENES_FOR_DRUG)


def is_hcmv_resistance_mutation(gene: str, substitution: str, drug: str | None = None) -> bool:
    """True iff `substitution` in `gene` is a catalogued resistance mutation (optionally for a specific drug)."""
    cat = RESISTANCE_BY_GENE.get(gene, {})
    if substitution not in cat:
        return False
    if drug is None:
        return True
    d = _MECH_ALIAS.get(drug.lower(), drug.lower())
    return d in cat[substitution]


def wt_consistency_ok() -> bool:
    """Self-integrity gate (v0 substitute for a committed-CDS translation check): every catalogued POINT
    mutation with the same (gene, position) must assert the SAME wild-type residue - a typo like a second WT
    at one position fails loudly. Deletion-form keys are opaque and skipped. `tests/test_hcmv_amr.py` pins it."""
    for gene, cat in RESISTANCE_BY_GENE.items():
        seen: dict[int, str] = {}
        for mut in list(cat) + list(BENIGN_BY_GENE.get(gene, frozenset())):
            m = _POINT.match(mut)
            if not m:
                continue
            wt, pos = m.group(1), int(m.group(2))
            if pos in seen and seen[pos] != wt:
                return False
            seen[pos] = wt
    return True


def call_hcmv_observed(drug: str, observed: dict[str, set[str]]) -> HCMVCall:
    """Deterministic HCMV R/S call from a virus's observed target-gene substitutions.

    `observed` = {gene: {substitutions}} keyed by 'UL97'/'UL54'/'UL56'. Rule (v0, MUTANT-LEVEL): R iff any of
    the DRUG'S target gene(s) carries >=1 catalogued resistance mutation that confers resistance TO THIS DRUG
    (cross-resistance is per-mutation catalogued, not assumed). An S call surfaces HCMV_UNDETECTABLE_MECHANISMS
    - it means 'no catalogued resistance mutation in the drug's target gene(s)', NOT 'definitely susceptible'."""
    genes = genes_for_hcmv_drug(drug)
    if not genes:
        return HCMVCall("INDETERMINATE", drug, [], [], "hcmv_target_site_v0",
                        f"no HCMV catalog for {drug!r} (supported: {', '.join(all_supported_hcmv_drugs())})")
    d = _MECH_ALIAS.get(drug.lower(), drug.lower())
    hits: list[str] = []
    for gene in genes:
        cat = RESISTANCE_BY_GENE[gene]
        for sub in sorted(observed.get(gene, set())):
            if sub in cat and d in cat[sub]:
                hits.append(f"{gene}:{sub}")
    pred = "R" if hits else "S"
    undetectable = HCMV_UNDETECTABLE_MECHANISMS if pred == "S" else []
    caveat = (f"deterministic HCMV target-site call ({'/'.join(genes)}, AD169 numbering; MUTANT-LEVEL v0). "
              f"Source: {_SOURCE}. IN-DISTRIBUTION vs the measured recombinant fold-change the catalog is "
              f"curated from (knowledge baseline; an independent number needs held-out phenotyping studies)."
              + (" An S call cannot rule out minor/novel substitutions, in-frame deletions, epistasis, "
                 "per-drug differences within a cross-resistance domain, or non-target-gene pathways."
                 if pred == "S" else ""))
    return HCMVCall(pred, drug, hits, undetectable, "hcmv_target_site_v0", caveat)
