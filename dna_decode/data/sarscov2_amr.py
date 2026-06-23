"""SARS-CoV-2 Mpro (3CLpro) antiviral resistance catalog — the next free-independent-label viral cell.

A SECOND coronavirus-grade target-site decoder (after influenza NA + HIV), and the project's next
free-independent-label cell: Mpro (3CL protease / nsp5) inhibitor resistance — **nirmatrelvir** (Paxlovid),
**ensitrelvir** (Xocova), **lufotrelvir** — via the established treatment/selection-emergent Mpro
substitutions.

WHY THIS ONE (the reason it exists): like the HIV cell, it can be VALIDATED against a FREE, variant-level,
INDEPENDENT laboratory genotype-phenotype resource — the **Stanford Coronavirus Antiviral & Resistance
Database (CoV-RDB)**, which curates lab-measured in-vitro susceptibility (fold-change EC50) + resistance
SELECTION results from peer-reviewed studies (the coronavirus analogue of Stanford HIVDB PhenoSense). The
phenotype is a wet-lab measurement INDEPENDENT of any genotype-interpretation rule, so scoring this catalog
against it clears the project's circular-label gate. Validation MUST use the measured fold-change, never a
genotype-derived prediction (rule-vs-rule = circular).

PROVENANCE (no fabrication): the Mpro major-substitution set is sourced VERBATIM from CoV-RDB's
`invitro_selection_results` Mpro-inhibitor selection records (Nirmatrelvir/Ensitrelvir/Lufotrelvir ×
_3CLpro), assembled at `data/raw/sarscov2/invitro_selection_antiviral.csv` (CoV-RDB
`hivdb/covid-drdb-payload`, MIT). The wild-type residue at every catalogued position is the SARS-CoV-2
Mpro reference (Wuhan-Hu-1 NC_045512.2 nsp5, 10055-10972, 306 aa) translation — committed at
`data/sarscov2_ref/SARSCoV2_Mpro_NC045512_cds.fna` and asserted residue-for-residue by
`tests/test_sarscov2_caller.py` (the same reference-integrity gate as the HIV/fungal cells; catalytic dyad
H41/C145 + the key nirmatrelvir position E166 verified).

NUMBERING: Mpro amino-acid positions 1-306 (Mpro = nsp5 protein numbering, N-terminus SGFRKM...). Mutation
form "<wt><pos><mut>" (e.g. E166V) — the same shorthand the HIV cell uses.

SCOPE (v0 — honest): CLASS-LEVEL + MUTANT-LEVEL. Every Mpro inhibitor maps to the SAME Mpro major-substitution
set (CoV-RDB defines selection at the per-drug level but the v0 cell is class-level, mirroring the HIV NNRTI
v0). MUTANT-LEVEL matching (not position-based) BECAUSE Mpro carries benign lineage polymorphisms at/near
resistance positions (e.g. P132H is fixed in Omicron and is NOT nirmatrelvir-resistant; K90R is a lineage
marker) — a position-based rule would over-call every Omicron genome. Two further honest caveats: (a) the
set is SELECTION-DERIVED, so some positions are weak/incidental resistance contributors — the fold-change
validation quantifies precision; (b) per-drug differential resistance (nirmatrelvir vs ensitrelvir) is real
and is v0.1. RdRp/remdesivir (a different target gene with a -1 ribosomal frameshift) is a separate v0.1 cell.
"""
from __future__ import annotations

from dataclasses import dataclass

# Wild-type residue at each catalogued Mpro position — the committed reference (Wuhan-Hu-1 NC_045512.2 nsp5)
# translation. `tests/test_sarscov2_caller.py` re-translates the reference and asserts it matches this map at
# every catalogued position (the integrity gate; a swapped/corrupt reference fails loudly).
MPRO_WT: dict[int, str] = {
    2: "G", 3: "F", 6: "M", 21: "T", 32: "L", 46: "S", 49: "M", 50: "L", 51: "N", 54: "Y",
    80: "H", 89: "L", 90: "K", 98: "T", 99: "P", 106: "I", 108: "P", 116: "A", 118: "Y", 126: "Y",
    127: "Q", 128: "C", 129: "A", 135: "T", 138: "G", 139: "S", 140: "F", 141: "L", 142: "N", 144: "S",
    147: "S", 155: "D", 160: "C", 166: "E", 167: "L", 169: "T", 171: "V", 172: "H", 173: "A", 186: "V",
    188: "R", 191: "A", 192: "Q", 193: "A", 194: "A", 197: "D", 200: "I", 202: "V", 203: "N", 204: "V",
    208: "L", 216: "D", 219: "F", 222: "R", 252: "P", 256: "Q", 294: "F", 295: "D", 296: "V", 298: "R",
    299: "Q", 301: "S", 304: "T", 305: "F",
}

# Mpro major resistance substitutions (MUTANT-LEVEL), VERBATIM from CoV-RDB Mpro-inhibitor selection records.
# Includes the canonical nirmatrelvir set (E166V/A, L50F, S144A, A173V/T, H172Y, L167F, T21I, ...).
MPRO_MAJOR_DRMS: frozenset[str] = frozenset({
    "G2V", "F3S", "M6I", "M6V", "T21I", "L32I", "S46F", "M49L", "L50F", "N51Y",
    "Y54C", "Y54S", "H80P", "L89F", "K90R", "T98I", "P99L", "I106S", "P108S", "A116T",
    "Y118H", "Y126F", "Q127R", "C128Y", "A129V", "T135I", "G138S", "S139A", "S139P", "F140L",
    "L141F", "N142G", "N142R", "S144A", "S147Y", "D155A", "D155G", "C160F", "C160Y", "E166A",
    "E166V", "L167F", "T169I", "V171S", "H172Q", "H172Y", "A173T", "A173V", "V186A", "R188G",
    "A191S", "A191T", "A191V", "Q192R", "A193P", "A194S", "D197Y", "I200T", "V202F", "N203D",
    "N203H", "N203K", "N203S", "V204F", "L208W", "D216G", "D216Y", "F219S", "R222Q", "P252L",
    "Q256L", "F294L", "F294V", "D295G", "D295N", "D295Y", "V296G", "R298G", "Q299K", "Q299P",
    "S301P", "T304A", "T304I", "F305L",
})

# v0 Mpro-inhibitor drugs (class-level — all share the same Mpro major-substitution set).
SARSCOV2_MPRO_DRUGS = ("nirmatrelvir", "ensitrelvir", "lufotrelvir")
SARSCOV2_DRUG_CLASS = {d: "Mpro_inhibitor" for d in SARSCOV2_MPRO_DRUGS}

# Mechanisms an Mpro major-substitution scan CANNOT see — an S (no-major-substitution) call cannot rule
# these out (the SARS-CoV-2 analogue of the HIV minor-DRM / bacterial efflux blind spot).
SARSCOV2_UNDETECTABLE_MECHANISMS = sorted({
    "minor_or_accessory_Mpro_substitutions",
    "novel_or_uncatalogued_Mpro_substitution",
    "per_drug_differential_resistance",          # v0 is class-level (nirmatrelvir vs ensitrelvir differ) -> v0.1
    "RdRp_remdesivir_or_other_target_resistance",  # different target gene / drug class (RdRp is v0.1)
    "non_Mpro_resistance_pathway",
    "mixture_or_minority_variant_below_consensus",
})

_SOURCE = ("Stanford CoV-RDB invitro_selection_results (Mpro-inhibitor in-vitro selection: "
           "Nirmatrelvir/Ensitrelvir/Lufotrelvir × 3CLpro; hivdb/covid-drdb-payload, MIT); "
           "Mpro WT = Wuhan-Hu-1 NC_045512.2 nsp5 reference translation")


@dataclass(frozen=True)
class SARSCoV2Call:
    """A deterministic SARS-CoV-2 Mpro R/S call (shape mirrors HIVCall / FungalCall — duck-typed by the CLI
    record builder `_target_site_record`)."""
    prediction: str            # R / S / INDETERMINATE
    drug: str
    determinants: list[str]    # e.g. ["Mpro:E166V"]
    undetectable_mechanisms: list[str]
    rule: str
    caveat: str


def gene_for_sarscov2_drug(drug: str) -> str | None:
    """Target gene for a SARS-CoV-2 drug (Mpro for every v0 Mpro inhibitor). Routes the caller + reference."""
    return "Mpro" if drug.lower() in SARSCOV2_DRUG_CLASS else None


def supported_sarscov2_drugs() -> list[str]:
    return sorted(SARSCOV2_MPRO_DRUGS)


def all_supported_sarscov2_drugs() -> list[str]:
    """Every SARS-CoV-2 drug the decoder routes (v0 = Mpro inhibitors)."""
    return sorted(SARSCOV2_MPRO_DRUGS)


def is_mpro_major_drm(substitution: str) -> bool:
    """True iff `substitution` ('<wt><pos><mut>', e.g. 'E166V') is a catalogued Mpro major resistance mutation."""
    return substitution in MPRO_MAJOR_DRMS


def call_sarscov2_observed(drug: str, observed: dict[str, set[str]]) -> SARSCoV2Call:
    """Deterministic Mpro R/S call from a virus's observed Mpro substitutions.

    `observed` = {gene: {substitutions}} keyed by 'Mpro'. Rule (v0, class-level, MUTANT-LEVEL): R iff Mpro
    carries >=1 catalogued major resistance substitution. An S call surfaces SARSCOV2_UNDETECTABLE_MECHANISMS
    — it means 'no major Mpro substitution found', NOT 'definitely susceptible' (minor/novel substitutions,
    per-drug differences, and RdRp/other-target mechanisms are invisible to v0)."""
    if drug.lower() not in SARSCOV2_DRUG_CLASS:
        return SARSCoV2Call("INDETERMINATE", drug, [], [], "sarscov2_mpro_major_drm_v0",
                            f"no SARS-CoV-2 catalog for {drug!r}")
    hits = sorted(f"Mpro:{s}" for s in observed.get("Mpro", set()) if s in MPRO_MAJOR_DRMS)
    pred = "R" if hits else "S"
    undetectable = SARSCOV2_UNDETECTABLE_MECHANISMS if pred == "S" else []
    caveat = (f"deterministic Mpro major-substitution call (Mpro inhibitor, Wuhan-Hu-1 numbering; CLASS-LEVEL "
              f"+ MUTANT-LEVEL v0 — selection-derived set, per-drug differential resistance not yet modelled). "
              f"Source: {_SOURCE}. Validate against CoV-RDB measured fold-change, NEVER a genotype-derived "
              f"prediction (circular). "
              + ("An S call cannot rule out minor/novel Mpro substitutions, per-drug differences, or "
                 "RdRp/other-target resistance." if pred == "S" else ""))
    return SARSCoV2Call(pred, drug, hits, undetectable, "sarscov2_mpro_major_drm_v0", caveat)
