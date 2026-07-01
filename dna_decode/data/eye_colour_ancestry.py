"""rs12913832 ancestry-confound quantification (the D:-free first cut of M2).

The eye-colour v0/v0.1 cells carry a standing caveat: "rs12913832 is European-calibrated → ancestry-
confounded." This module turns that hand-wave into committed numbers. The blue allele (G) frequency across
1000 Genomes phase-3 SUPER-populations (sourced from Ensembl REST, variation/human/rs12913832?pops=1) shows
the blue allele is strongly EUROPEAN-CONCENTRATED — so in a Europe-majority cohort (OpenSNP), the SNP
genotype is partly a tag for European ancestry, which itself correlates with light eye colour.

HONEST SCOPE (load-bearing): this BOUNDS the confound's STRUCTURAL basis (the SNP is ancestry-informative);
it does NOT prove the v0 0.993 is inflated, because rs12913832 is ALSO a known CAUSAL variant (HERC2
regulatory element controlling OCA2 expression). The actual disentangler is the WITHIN-European re-score
(does the SNP predict eye colour AMONG Europeans only) — that needs the per-user OpenSNP genotypes + an
ancestry axis and is D:-gated (deferred). This module is the confound's magnitude, not its resolution.
"""
from __future__ import annotations

import re

BLUE_ALLELE = "G"   # rs12913832 G = blue (matches the validated v0 strand-agnostic rule: blue={G,C})

# Blue-allele (G) frequency by 1000G phase-3 super-population. SOURCED: Ensembl REST, fetched 2026-06-30.
RS12913832_BLUE_FREQ_1000G: dict[str, float] = {
    "EUR": 0.6362,
    "AMR": 0.2017,
    "SAS": 0.0706,
    "AFR": 0.0280,
    "EAS": 0.0020,
    "ALL": 0.1773,
}


# self-reported-ancestry keyword classifier (for the within-European re-score, M2-full).
# Conservative: a stratum is assigned ONLY when the text has that continent's keywords and NOT the other's;
# 'mixed'/ambiguous/Ashkenazi(distinct cluster)/bare-'american' -> MIXED_UNKNOWN (excluded from the clean test).
_EUROPEAN_KW = (
    "european", "caucasian", "caucasion", "british", "irish", "ireland", "scottish", "scotland",
    "english", "england", "welsh", "wales", "scandinav", "nordic", "finnish", "finland", "german",
    "french", "italian", "spanish", "iberian", "portuguese", "russian", "polish", "czech", "dutch",
    "norwegian", "swedish", "danish", "greek", "hungarian", "austrian", "romanian", "ukrain", "slav",
)
_NON_EUROPEAN_KW = (
    "african", "hispanic", "latino", "mexican", "asian", "chinese", "japanese", "korean",
    "south asian", "east asian", "native american", "amerindian", "middle eastern", "arab",
    "persian", "filipino", "vietnamese", "thai", "caribbean", "black", "indian", "berber", "hmong",
)
_AMBIGUOUS_KW = ("mixed", "ashkenazi", "jewish", "adopted")


def classify_self_reported_ancestry(*texts: str) -> str:
    """Classify free-text self-reported ancestry -> EUROPEAN / NON_EUROPEAN / MIXED_UNKNOWN.

    Conservative by design: EUROPEAN only if a European keyword is present AND no non-European keyword;
    symmetric for NON_EUROPEAN; anything mixed/ambiguous/both/neither -> MIXED_UNKNOWN (excluded from the
    clean within-ancestry test). This is a keyword proxy for SELF-REPORTED ancestry, not genetic inference."""
    blob = " ".join(t for t in texts if t).lower()
    if not blob.strip():
        return "MIXED_UNKNOWN"
    # word-boundary match: prevents 'asian' matching inside 'caucasian' (a real OpenSNP collision).
    def _hit(kws):
        return any(re.search(r"\b" + re.escape(k) + r"\b", blob) for k in kws)
    if _hit(_AMBIGUOUS_KW):
        return "MIXED_UNKNOWN"
    has_eur = _hit(_EUROPEAN_KW)
    has_non = _hit(_NON_EUROPEAN_KW)
    if has_eur and not has_non:
        return "EUROPEAN"
    if has_non and not has_eur:
        return "NON_EUROPEAN"
    return "MIXED_UNKNOWN"


def confound_summary() -> dict:
    """Quantify how ancestry-informative rs12913832 is (the confound's structural magnitude)."""
    superpops = {k: v for k, v in RS12913832_BLUE_FREQ_1000G.items() if k != "ALL"}
    eur = superpops["EUR"]
    others = {k: v for k, v in superpops.items() if k != "EUR"}
    min_pop = min(others, key=others.__getitem__)
    return {
        "blue_allele": BLUE_ALLELE,
        "eur_blue_freq": eur,
        "max_nonEUR_blue_freq": max(others.values()),
        "min_blue_freq_pop": min_pop,
        "min_blue_freq": others[min_pop],
        "eur_over_min_ratio": round(eur / others[min_pop], 1) if others[min_pop] else None,
        "is_ancestry_informative": eur >= 0.5 and max(others.values()) <= 0.25,
        "interpretation": (
            "blue allele European-concentrated -> rs12913832 is ancestry-informative; in a Europe-majority "
            "cohort it partly tags ancestry. NOT a proof of inflated accuracy (the SNP is also causal at "
            "HERC2/OCA2). Disentangler = within-European re-score (D:-gated, deferred)."
        ),
        "source": "Ensembl REST 1000G phase_3 super-pop frequencies, fetched 2026-06-30",
    }
