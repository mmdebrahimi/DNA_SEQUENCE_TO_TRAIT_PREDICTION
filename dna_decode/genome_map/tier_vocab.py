"""Tier-boundary wording vocabulary (Step 3).

The low-confidence + unknown product-wording patterns that separate
`curated-molecular-function` from `homology-only-hypothesis` from `unknown`.
These are SEEDED + validated against the Step-2 tool-surface manifest's REAL
db-light product vocabulary (not guessed) — see the wiki manifest's
`bakta.product_vocabulary_sample` / `low_confidence_product_examples`.

Kept as a separate module so the patterns are auditable + the Step-2 manifest
can be diffed against them when the spike runs.
"""
from __future__ import annotations

# Wording that signals a WEAK, homology-only inference (the
# `homology-only-hypothesis` tier). Lower-cased substring match on the product.
# Validated against Bakta db-light wording (manifest-seeded).
LOW_CONFIDENCE_PATTERNS: tuple[str, ...] = (
    "putative",
    "probable",
    "by similarity",
    "domain-containing protein",
    "domain-containing",
    "uncharacterized",
    "predicted",
    "duf",                # DUFxxxx = "domain of unknown function" family
    "-like protein",
    "similar to",
)

# Wording (or emptiness) that signals NO functional call at all (the `unknown`
# tier). The product is genuinely uninformative.
UNKNOWN_PATTERNS: tuple[str, ...] = (
    "hypothetical protein",
    "unknown function",
    "uncharacterized protein",
)


def matches_any(product: str, patterns: tuple[str, ...]) -> str | None:
    """Return the first pattern that is a substring of the lower-cased product, else None."""
    p = (product or "").lower()
    for pat in patterns:
        if pat in p:
            return pat
    return None
