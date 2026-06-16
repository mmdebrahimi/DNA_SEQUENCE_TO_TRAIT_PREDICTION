"""NON-FROZEN experimental drug-rule overlay — scorer-local rules NOT in the deployed surface.

The reproducibility freeze (2026-06-13) froze the 6-drug deployed decoder: `amr_rules.DRUG_RULE`
(`threshold + subclass_any`, count/OR semantics) + `mic_tiers.DRUG_BREAKPOINTS`/`supported_drugs`.
This module is the EXPERIMENTAL overlay for candidate drugs whose rule shape the frozen engine cannot
express and which are validated ONLY on the external-validation arm (never the frozen report card /
shipped_decoder_surface). Adding a rule here touches NO frozen file.

First entry: trimethoprim-sulfamethoxazole (TMP-SMX / co-trimoxazole). Its determinant rule is
`(>=1 acquired sul gene) AND (>=1 acquired dfr gene)` — AND-across-two-gene-families, a shape the
frozen `DRUG_RULE` (count of determinants matching ANY subclass token) cannot represent. Empirically
required (Sci234, measured TRISUL MIC): sul+dfr 69/70 R; sul-only 0/40 R; dfr-only 1/6 R; neither
1/117 R -> an OR rule would mis-call the sul-only/dfr-only strata; AND separates cleanly.

The MIC tiering reuses the FROZEN, pure `mic_tiers.classify_tier` with breakpoints passed in (no frozen
edit): co-trimoxazole Enterobacterales breakpoints are expressed as the TRIMETHOPRIM COMPONENT.
"""
from __future__ import annotations

import re

from dna_decode.data.mic_tiers import classify_tier

# Branding — every artifact a scorer emits from this overlay carries these so a cell is never confused
# with a frozen deployed-decoder claim.
RULE_STATUS = "EXPERIMENTAL_SCORED"
RULE_SCOPE = "scorer_local"
DRUG = "trimethoprim-sulfamethoxazole"

# Co-trimoxazole Enterobacterales breakpoints, TRIMETHOPRIM-COMPONENT units.
# Source to authoritatively pin before any deployed claim: EUCAST v14.0 (Enterobacterales TMP-SMX
# S <= 2 / R > 4) + CLSI M100 2024 (S <= 2/38, R >= 4/76). EUCAST==CLSI here so classify_tier never
# returns AMBIGUOUS on this drug. [unverified-in-plan — see plan Risk Flags]
COTRIMOXAZOLE_BREAKPOINTS: dict[str, float] = {
    "clsi_r": 4.0, "clsi_s": 2.0, "eucast_r": 4.0, "eucast_s": 2.0,
}

# Explicit acquired-family regexes (NOT loose startswith): match the curated acquired sulfonamide
# (sul1-4) + trimethoprim (dfrA/dfrB + numbered allele) determinants, after stripping any allele/suffix
# decoration. Deliberately EXCLUDE regulators / look-alikes: `sulR`, `sulP`, a bare `dfrA` without an
# allele number, `folP`/`folA` (target point-mutation genes — a different, AMRFinder-invisible mechanism).
SUL_RE = re.compile(r"^sul[1-4]$", re.IGNORECASE)
DFR_RE = re.compile(r"^dfr[AB]\d+[a-z]?$", re.IGNORECASE)


def _normalize(symbol: str) -> str:
    """Strip common allele/suffix decoration so a family regex can match (e.g. 'sul2_1' -> 'sul2')."""
    s = (symbol or "").strip()
    # drop a trailing _<digits> allele-copy suffix some callers append (sul2_1, dfrA17_2)
    s = re.sub(r"_\d+$", "", s)
    return s


def is_sul(symbol: str) -> bool:
    return bool(SUL_RE.match(_normalize(symbol)))


def is_dfr(symbol: str) -> bool:
    return bool(DFR_RE.match(_normalize(symbol)))


def tmp_smx_call(gene_symbols: list[str]) -> dict:
    """Scorer-local TMP-SMX rule: R iff (>=1 sul) AND (>=1 dfr), else S.

    Returns prediction + the matched determinants + the exact rule text (for the artifact)."""
    matched_sul = sorted({g for g in gene_symbols if is_sul(g)})
    matched_dfr = sorted({g for g in gene_symbols if is_dfr(g)})
    prediction = "R" if (matched_sul and matched_dfr) else "S"
    return {
        "prediction": prediction,
        "matched_sul": matched_sul,
        "matched_dfr": matched_dfr,
        "rule_text": "R iff (>=1 acquired sul[1-4]) AND (>=1 acquired dfrA/B<allele>); else S "
                     "(folP/folA target point-mutation TMP-R is an acquired-gene blind-spot)",
    }


def cotrimoxazole_tier(mic_tokens: list[float], distinct_calls: set[str] | None = None) -> str:
    """MIC-tier for co-trimoxazole via the FROZEN classify_tier + the overlay breakpoints (no frozen edit)."""
    return classify_tier(list(mic_tokens), set(distinct_calls or set()), COTRIMOXAZOLE_BREAKPOINTS)
