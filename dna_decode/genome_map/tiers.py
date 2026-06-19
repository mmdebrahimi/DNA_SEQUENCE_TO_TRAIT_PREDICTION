"""Function/confidence tier classifier (Step 3).

Assigns every Bakta feature its single highest-confidence NON-determinant tier
from the v1 4-tier schema. The top tier — `determinant-phenotype` — is NOT
assigned here: it is applied by the overlay/assembler (Step 4/5) ONLY when a
HIGH-confidence AMRFinder determinant join hits the feature. This module
produces the tier for everything else:

    curated-molecular-function  — a named gene_symbol OR a specific product
    homology-only-hypothesis    — low-confidence wording (putative/by similarity/DUF…)
    unknown                     — hypothetical / empty product, no gene call

Returns the (tier, classification_reason) pair; the REASON (the matched signal)
is retained in the map so Step 6's G1 audit computes from the map alone
(brainstorm catch C4). Vocabulary is seeded from the Step-2 manifest
(`tier_vocab`).
"""
from __future__ import annotations

from dna_decode.genome_map import (
    TIER_CURATED_FUNCTION,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.tier_vocab import (
    LOW_CONFIDENCE_PATTERNS,
    UNKNOWN_PATTERNS,
    matches_any,
)


def classify_feature_tier(product: str, gene_symbol: str = "") -> tuple[str, str]:
    """Classify one feature into a NON-determinant tier + return the reason.

    Precedence (high -> low), honesty-first:
      1. A firm `gene_symbol` (a named ortholog call) -> curated-molecular-function.
         A named gene is a confident assignment; the product wording is kept as
         secondary evidence by the caller.
      2. A `hypothetical protein`/`unknown function` product (or empty + no gene)
         -> unknown. Genuinely uninformative; never relabelled to a lower tier.
      3. Low-confidence wording (`putative`/`by similarity`/`DUF…`) -> the
         homology-only-hypothesis tier (this is the G1 type-(a) DEMOTION).
      4. Any remaining non-empty (specific) product -> curated-molecular-function.
      5. Otherwise -> unknown.

    Returns (tier, classification_reason). The reason names the matched signal so
    the demotion / call is auditable from the map.
    """
    product = (product or "").strip()
    gene_symbol = (gene_symbol or "").strip()

    # 1. A firm gene-symbol call is curated regardless of product wording.
    if gene_symbol:
        return TIER_CURATED_FUNCTION, f"named gene symbol '{gene_symbol}'"

    # 2. Explicitly unknown wording (or no product at all) -> unknown.
    unk = matches_any(product, UNKNOWN_PATTERNS)
    if unk is not None:
        return TIER_UNKNOWN, f"unknown-wording product ('{unk}')"
    if not product:
        return TIER_UNKNOWN, "empty product, no gene symbol"

    # 3. Low-confidence wording -> homology-only-hypothesis (the honest demotion).
    low = matches_any(product, LOW_CONFIDENCE_PATTERNS)
    if low is not None:
        return TIER_HOMOLOGY_HYPOTHESIS, f"low-confidence wording ('{low}')"

    # 4. A specific product with no weak wording is a confident functional call.
    return TIER_CURATED_FUNCTION, "specific product wording"
