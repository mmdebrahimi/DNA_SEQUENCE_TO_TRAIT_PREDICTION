"""Human eye-colour deterministic decoder — the first OFF-PATHOGEN cell (2026-06-28).

The wide-net scan (`wiki/wide_net_gene_phenotype_source_scan_2026-06-28.md`) flagship: the deterministic
gene->trait pattern (curated rule + free measured labels), generalized from AMR to a human visible trait.
This is the RULE side (the OpenSNP measured-label side is `scripts/eye_colour_opensnp_validate.py`).

v0 = a SINGLE-LOCUS rule on **rs12913832** (HERC2 intron-86 enhancer regulating OCA2) — the single strongest
blue/brown eye-colour predictor (alone explains ~74% of blue-vs-brown variance in Europeans). This is the
honest "major-effect v0" (like the position-based v0 in the HIV arc). The full **IrisPlex 6-SNP** model
(rs12913832 + rs1800407 OCA2 + rs12896399 SLC24A4 + rs16891982 SLC45A2 + rs1393350 TYR + rs12203592 IRF4)
is the v0.1 upgrade — but its multinomial-logit COEFFICIENTS must be SOURCED from Walsh 2011 (FSI:Genetics),
NOT fabricated; v0.1 is deliberately deferred until those coefficients are in hand.

STRAND DISCIPLINE (load-bearing — sourced 2026-06-28, caught a memory-inversion): in dbSNP forward orientation
GG->blue, AA/AG->brown; but 23andMe/many DTC chips report rs12913832 on the COMPLEMENTARY strand as C/T
(CC->blue, TT/CT->brown). rs12913832 is an A/G SNP = NON-palindromic, so a strand-AGNOSTIC allele-class map is
unambiguous: the BLUE allele is {G (dbSNP), C (complement)}, the BROWN allele is {A (dbSNP), T (complement)}.
We classify by allele CLASS, never by raw letter, so the call is correct whichever strand the raw file uses.

HONEST BOUNDS: ~3% of GG/CC (blue-allele-homozygous) Europeans have brown eyes; the rule is European-calibrated
(ancestry-confounded outside it) — so validation MUST be within-ancestry (the within-lineage de-confounding
discipline). Intermediate/green eyes are a genuine third class the single locus only weakly resolves (the
heterozygote bin); v0 reports `intermediate` for heterozygotes rather than forcing blue/brown.
"""
from __future__ import annotations

EYE_COLOUR_LOCUS = "rs12913832"
# strand-agnostic allele classes (A/G non-palindromic -> no ambiguity). G(dbSNP)=C(complement)=blue allele.
_BLUE_ALLELES = frozenset({"G", "C"})
_BROWN_ALLELES = frozenset({"A", "T"})

# The IrisPlex 6-SNP panel (rsIDs are standard; v0.1 needs the Walsh-2011 coefficients SOURCED, not invented).
IRISPLEX_SNPS = ("rs12913832", "rs1800407", "rs12896399", "rs16891982", "rs1393350", "rs12203592")


def _alleles(genotype: str) -> list[str]:
    """Normalize a raw genotype string ('GG', 'A/G', 'CT', '--', 'GA') -> list of single-letter alleles."""
    if not genotype:
        return []
    return [c.upper() for c in genotype if c.upper() in ("A", "C", "G", "T")]


def call_eye_colour(genotype: str) -> dict:
    """Single-locus rs12913832 eye-colour call (strand-agnostic). genotype = the raw 2-allele call.

    blue-allele {G,C} homozygote -> blue; brown-allele {A,T} homozygote -> brown; heterozygote ->
    intermediate; missing/no-call -> INDETERMINATE. Returns the prediction + the matched allele classes.
    """
    al = _alleles(genotype)
    if len(al) != 2:
        return {"locus": EYE_COLOUR_LOCUS, "genotype": genotype, "prediction": "INDETERMINATE",
                "reason": "not a 2-allele call", "n_blue_alleles": None, "n_brown_alleles": None}
    n_blue = sum(1 for a in al if a in _BLUE_ALLELES)
    n_brown = sum(1 for a in al if a in _BROWN_ALLELES)
    if n_blue + n_brown != 2:                       # an allele outside {A,C,G,T-as-class} -> can't classify
        pred = "INDETERMINATE"
    elif n_blue == 2:
        pred = "blue"
    elif n_brown == 2:
        pred = "brown"
    else:
        pred = "intermediate"                       # heterozygote: single locus only weakly resolves
    return {"locus": EYE_COLOUR_LOCUS, "genotype": genotype, "prediction": pred,
            "n_blue_alleles": n_blue, "n_brown_alleles": n_brown,
            "rule": "rs12913832 single-locus v0 (strand-agnostic; blue={G,C}, brown={A,T})"}
