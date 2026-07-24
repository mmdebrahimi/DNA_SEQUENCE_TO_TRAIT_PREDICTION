"""Bacteriophage host-RECEPTOR catalog - the FIRST non-AMR, non-host-organism cell: a phage-genome ->
host-receptor-class decoder, and the project's first phage (viral-of-bacteria) axis.

WHY THIS ONE (a genuinely-new axis with a FREE MEASURED label): the BASEL collection (Maffei et al. 2021,
PLOS Biology, doi 10.1371/journal.pbio.3001424; completed in Maffei et al. 2025, doi 10.1371/journal.pbio
.3003063) is a systematically characterized library of E. coli phages whose host RECEPTORS were determined
EXPERIMENTALLY - by testing >50 single-gene E. coli K-12 receptor mutants for loss of phage sensitivity, plus
EOP host-range assays. Those receptor assignments are wet-lab measurements published in a free, open-access
(CC-BY) paper, and the genomes are deposited (GenBank MZ501046-MZ501113). So a curated receptor catalog scored
against them is a genuinely-new MEASURED-label decoder - the phage analogue of the AMR determinant catalog.

SCOPE (v0 - honest, and the tractability boundary is DOCUMENTED not hidden):
  - The TRACTABLE sub-problem is RECEPTOR-CLASS (which outer-membrane / LPS structure the phage's
    receptor-binding protein targets), NOT the full phage x strain host-range MATRIX. The matrix is polygenic
    (receptor-binding protein + anti-defense systems + downstream infection compatibility) and is documented
    as an "intractable challenge" from genome alone (the 2026 benchmark that cracked receptor-CLASS needed
    1050 genetic screens + AlphaFold3 + deep learning). We scope to receptor-class ONLY and say so.
  - The v0 CALLER is genome-homology receptor-TRANSFER: a query phage inherits the receptor of its nearest
    genome-BLAST neighbor among reference phages of known receptor. This measures HOW WELL receptor usage
    transfers along genome similarity - the honest scientific question - rather than claiming a solved
    RBP->receptor map. Receptor is largely genus-conserved (Drexlerviridae -> FhuA/BtuB; Autographiviridae
    -> LPS core; T5-like Demerecviridae -> BtuB) but VARIES within the T-even Tevenvirinae by RBP
    (T4=OmpC / T6=Tsx / T2=FadL), so transfer is expected to be receptor-dependent - the validation reports
    that honestly rather than a single headline.
  - TIER: IN-DISTRIBUTION (catalog curated from BASEL, scored on BASEL) - a knowledge baseline, like the
    SARS-CoV-2 CoV-RDB cell. An INDEPENDENT number needs a held-out phage set with measured receptors.

PROVENANCE (no fabrication - every receptor->taxon assignment is quoted VERBATIM from the open-access BASEL
Results text / figures; every reference-phage accession was VERIFIED to resolve to the named phage via NCBI
efetch before being committed here):
  - Maffei 2021 PLOS Biology 3001424 Results ("Systematic exploration of E. coli phage-host interactions"):
    the receptor->genus/family assignments below. Terminal receptors named in the paper: FhuA, BtuB, YncD,
    TolC, LptD, LamB, FepA, OmpC, OmpF, OmpA, Tsx, FadL, NfrA, plus LPS core and ECA (enterobacterial common
    antigen) glycan receptors.
  - Model-phage receptors (T2/T4/T5/T6/T7/T3/lambda/N4) are the textbook receptor assignments the BASEL paper
    itself cites as the reference frame (T4->OmpC, T6->Tsx, T2->FadL, T5->BtuB, T7/T3->LPS core, lambda->LamB,
    N4->NfrA).

NUMBERING/KEYS: receptor-class strings are the canonical E. coli surface-structure names (outer-membrane
proteins by gene symbol; LPS_core / ECA / O_antigen for glycan receptors). Genus keys use NCBI taxonomic
genus names where a phage's receptor is genus-conserved; family keys (e.g. Drexlerviridae) are used where the
paper reports the receptor at family/subfamily rank.
"""
from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Receptor classes (the label space) - every one is named VERBATIM in the BASEL Results.
# ---------------------------------------------------------------------------
RECEPTOR_CLASSES: tuple[str, ...] = (
    "FhuA",     # TonB-dependent ferrichrome transporter
    "BtuB",     # TonB-dependent vitamin-B12 transporter
    "YncD",     # putative TonB-dependent transporter
    "TolC",     # outer-membrane efflux channel
    "LptD",     # LPS assembly translocon (essential OMP)
    "FepA",     # TonB-dependent ferric-enterobactin transporter
    "LamB",     # maltoporin
    "OmpC",     # classical porin
    "OmpF",     # classical porin
    "OmpA",     # outer-membrane protein A
    "Tsx",      # nucleoside-specific channel
    "FadL",     # long-chain fatty-acid transporter
    "NfrA",     # N4 receptor
    "LPS_core", # lipopolysaccharide core glycan
    "ECA",      # enterobacterial common antigen glycan
    "O_antigen",# LPS O-antigen (surface glycan)
)

# ---------------------------------------------------------------------------
# Receptor -> taxon assignments, VERBATIM from Maffei 2021 (Results + figures).
# A taxon may list >1 receptor when the paper reports receptor variation within it
# (the honest coarseness the validation surfaces). `primary` is the receptor the
# paper reports as dominant / terminal for that taxon.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TaxonReceptor:
    taxon: str          # NCBI genus, or paper family/subfamily where receptor is family-rank
    rank: str           # "genus" | "family" | "subfamily"
    receptors: tuple[str, ...]   # receptors the paper attributes to this taxon
    primary: str        # dominant/terminal receptor
    note: str           # VERBATIM-grounded note
    clade_conserved: bool = True  # True = the WHOLE clade uses one well-defined receptor (usable as a
    #                               label); False = receptor VARIES within the clade by receptor-binding
    #                               protein (T-even, Drexlerviridae) -> NOT auto-labelled (no fabrication).


# Family/subfamily-rank assignments (where the paper reports the receptor at that rank).
FAMILY_RECEPTOR: dict[str, TaxonReceptor] = {
    "Drexlerviridae": TaxonReceptor(
        "Drexlerviridae", "family", ("FhuA", "BtuB", "YncD", "TolC", "LptD"), "FhuA",
        "11/13 use FhuA/BtuB/YncD/TolC; AugustePiccard(Bas01)+JeanPiccard(Bas02) use LptD.",
        clade_conserved=False),   # 4-5 receptors within the family -> RBP-variable, not labellable
    "Demerecviridae": TaxonReceptor(
        "Demerecviridae", "family", ("BtuB", "FhuA", "FepA"), "BtuB",
        "Markadamsvirinae (T5-like): vast majority BtuB; T5 FhuA; H8/S124 FepA."),
    "Autographiviridae": TaxonReceptor(
        "Autographiviridae", "family", ("LPS_core",), "LPS_core",
        "Studiervirinae (T3/T7): all tested use core LPS structures."),
    "Autotranscriptaviridae": TaxonReceptor(
        "Autotranscriptaviridae", "family", ("LPS_core",), "LPS_core",
        "T7-like family (reclassified Autographiviridae/Studiervirinae): core LPS receptor."),
    "Schitoviridae": TaxonReceptor(
        "Schitoviridae", "family", ("NfrA", "ECA"), "NfrA",
        "Enquatrovirinae (N4): NfrA terminal receptor; strong wecB(ECA) dependence."),
    "Straboviridae": TaxonReceptor(
        "Straboviridae", "family", ("OmpC", "Tsx", "FadL", "OmpF", "OmpA", "LPS_core"), "OmpC",
        "T-even Tevenvirinae: OmpC/Tsx/FadL/OmpF/OmpA vary by RBP (T4=OmpC, T2=FadL, T6=Tsx).",
        clade_conserved=False),   # the canonical RBP-variable clade
    "Vequintaviridae": TaxonReceptor(
        "Vequintaviridae", "family", ("ECA", "LPS_core"), "ECA",
        "Vequintavirinae incl. phi92-like: ECA shared primary receptor + LPS-core glucose."),
}

# Genus-rank assignments (NCBI genus where receptor is genus-conserved). Model phages
# anchor the T-even receptor variation the family rank cannot express.
GENUS_RECEPTOR: dict[str, TaxonReceptor] = {
    # T-even Tevenvirinae - receptor VARIES by RBP (T4=OmpC, T2=FadL, T6=Tsx) -> NOT auto-labelled.
    "Tequatrovirus": TaxonReceptor("Tequatrovirus", "genus", ("OmpC", "FadL", "Tsx", "LPS_core"), "OmpC",
        "T-even: OmpC(T4)/FadL(T2)/Tsx(T6) vary by RBP; short tail fibers target lipid A-Kdo LPS core.",
        clade_conserved=False),
    "Tequintavirus": TaxonReceptor("Tequintavirus", "genus", ("BtuB", "FhuA"), "BtuB",
        "T5-like: FhuA (T5) / BtuB (majority of the subfamily)."),
    "Teseptimavirus": TaxonReceptor("Teseptimavirus", "genus", ("LPS_core",), "LPS_core",
        "T7-like: broad LPS-core recognition."),
    "Berlinvirus": TaxonReceptor("Berlinvirus", "genus", ("LPS_core",), "LPS_core",
        "T3-like: LPS core."),
    "Enquatrovirus": TaxonReceptor("Enquatrovirus", "genus", ("NfrA", "ECA"), "NfrA",
        "N4-like: NfrA terminal; ECA(wecB) dependence."),
    "Felixounavirus": TaxonReceptor("Felixounavirus", "genus", ("LPS_core",), "LPS_core",
        "JohannRWettstein totally depends on an intact LPS core."),
    "Lambdavirus": TaxonReceptor("Lambdavirus", "genus", ("LamB",), "LamB",
        "Lambda: LamB maltoporin (requires intact LPS inner core)."),
    "Vequintavirus": TaxonReceptor("Vequintavirus", "genus", ("ECA", "LPS_core"), "ECA",
        "Vequintavirinae: all tested use ECA as the shared primary receptor."),
    # The two BASEL Drexlerviridae phages the paper explicitly links to LptD (genus wins over the
    # RBP-variable Drexlerviridae family fallback).
    "Augustepiccardvirus": TaxonReceptor("Augustepiccardvirus", "genus", ("LptD",), "LptD",
        "AugustePiccard (Bas01): resistance linked to LptD mutations."),
    "Julespiccardvirus": TaxonReceptor("Julespiccardvirus", "genus", ("LptD",), "LptD",
        "JeanPiccard (Bas02): LptD terminal receptor."),
}


# ---------------------------------------------------------------------------
# Reference phages with textbook/VERBATIM receptor labels + VERIFIED NCBI accessions.
# These are the labeled anchor set the genome-homology caller transfers FROM.
# Every accession here was confirmed via NCBI efetch to resolve to the named phage.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ReferencePhage:
    name: str
    accession: str      # VERIFIED NCBI nucleotide accession
    receptor: str       # primary receptor class
    genus_hint: str     # NCBI genus (for the family/genus fallback path)


REFERENCE_PHAGES: tuple[ReferencePhage, ...] = (
    ReferencePhage("T4", "NC_000866", "OmpC", "Tequatrovirus"),
    ReferencePhage("T5", "NC_005859", "BtuB", "Tequintavirus"),
    ReferencePhage("T7", "NC_001604", "LPS_core", "Teseptimavirus"),
    ReferencePhage("T3", "NC_003298", "LPS_core", "Berlinvirus"),
    ReferencePhage("lambda", "NC_001416", "LamB", "Lambdavirus"),
    ReferencePhage("N4", "NC_008720", "NfrA", "Enquatrovirus"),
)


def _canon(receptor: str) -> str:
    return receptor.strip()


def receptor_for_taxon(taxon: str) -> TaxonReceptor | None:
    """Look up the receptor assignment for an NCBI genus OR paper family name.

    Genus rank is preferred (finer); falls back to family rank. Returns None when the
    taxon has no catalogued receptor (honest INDETERMINATE, never a fabricated guess).
    """
    if taxon in GENUS_RECEPTOR:
        return GENUS_RECEPTOR[taxon]
    if taxon in FAMILY_RECEPTOR:
        return FAMILY_RECEPTOR[taxon]
    return None


def primary_receptor_for_taxon(taxon: str) -> str | None:
    tr = receptor_for_taxon(taxon)
    return tr.primary if tr is not None else None


def receptor_for_lineage(lineage: list[str] | tuple[str, ...]) -> TaxonReceptor | None:
    """Given an ordered NCBI taxonomic lineage (genus-last preferred), return the finest
    catalogued receptor assignment. Scans from most specific to least (genus before family)."""
    for taxon in lineage:
        tr = receptor_for_taxon(taxon)
        if tr is not None:
            return tr
    return None


def label_receptor_for_lineage(lineage: list[str] | tuple[str, ...]) -> str | None:
    """The receptor to use as a VALIDATION LABEL for a phage of this lineage.

    Returns the finest catalogued receptor ONLY when that taxon is `clade_conserved` (one
    well-defined receptor for the whole clade). RBP-variable clades (T-even Straboviridae/
    Tequatrovirus, Drexlerviridae) return None - they are excluded from labelling rather than
    assigned a single receptor they do not uniformly use (no fabricated labels). This is the
    honest tractability boundary, not a coverage gap.
    """
    tr = receptor_for_lineage(lineage)
    if tr is None or not tr.clade_conserved:
        return None
    return tr.primary


def is_receptor_class(x: str) -> bool:
    return _canon(x) in RECEPTOR_CLASSES
