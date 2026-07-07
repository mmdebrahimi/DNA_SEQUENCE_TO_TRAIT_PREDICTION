"""CYP2D6 hybrid IDENTITY (v0.4) — a three-level abstaining classifier over the PSV D6-fraction profile.

Phase B of the hybrid arc: turns the read-level PSV evidence (per-region CYP2D6/CYP2D7 D6-fraction, from
`scripts/cyp2d6_psv_evidence.py`) into a hybrid-ALLELE call. Validated at full N on the structural panel
(`wiki/cyp2d6_psv_phaseb_falsifier.json`): non-hybrid specificity 1.0 (26/26 — a normal/dup/deletion is NEVER
mis-called a hybrid), *68 4/4, *36 6/8 (the 2 misses are subtle exon-9 gene-conversions).

THREE LEVELS (brainstorm R2 — never force a hybrid-positive into a specific allele):
  * resolved identity: `*68` (5' CYP2D6 / 3' CYP2D7) | `*13` (opposite; 5' CYP2D7 / 3' CYP2D6) |
    `*36` (exon-9-tip CYP2D7 conversion)
  * `hybrid_present_identity_unresolved` — a hybrid signal too weak/ambiguous to name (e.g. a subtle *36)
  * `evidence_not_callable` — insufficient callable PSVs (the callability gate)

HONEST TIER: high-SPECIFICITY (a resolved call is trustworthy — spec 1.0 in validation). Sensitivity is
per-allele: *68 clean, *36 partial (subtle conversions -> unresolved), *13 single-sample-validated
(UNPOWERED — n=1). NEVER a clinical tool. The D6-fraction is already copy-normalized (a fraction, not an
absolute depth), so a pure duplication's flat-but-elevated profile does not trip the directional signal.
Reference tool: Cyrius (117-PSV, WGS).
"""
from __future__ import annotations

import statistics

# Validated thresholds (full-N panel: spec 1.0). The D6-fraction directional shift + the exon-9-tip dip.
DIRECTIONAL_THRESHOLD = 0.15     # |5'-mean − 3'-mean| >= this -> a directional fusion (*68 / *13)
EXON9_DIP_THRESHOLD = 0.15       # (gene-median − downstream_exon9) >= this -> a *36 exon-9 conversion
MIN_CALLABLE_PSV = 50            # callability gate: below this -> evidence_not_callable

_FIVE_PRIME = ("exon1", "upstream_exon1", "intron1")           # 5' end of CYP2D6 (high genomic coord)
_THREE_PRIME = ("downstream_exon9", "intron6", "exon6")        # 3' end (low genomic coord)


def _mean(region_d6_fraction: dict, keys) -> float | None:
    v = [region_d6_fraction[k] for k in keys if region_d6_fraction.get(k) is not None]
    return statistics.mean(v) if v else None


def classify_hybrid_identity(region_d6_fraction: dict, n_callable: int,
                             directional: float = DIRECTIONAL_THRESHOLD,
                             exon9_dip: float = EXON9_DIP_THRESHOLD,
                             min_callable: int = MIN_CALLABLE_PSV) -> dict:
    """Region D6-fraction profile -> three-level hybrid-identity call. `region_d6_fraction` maps region ->
    median D6-fraction (from the PSV evidence table); `n_callable` is the count of callable PSVs."""
    if n_callable < min_callable:
        return {"call": "evidence_not_callable", "confidence": "n/a",
                "note": f"only {n_callable} callable PSVs (< {min_callable}); cannot resolve hybrid identity",
                "features": {}}
    five = _mean(region_d6_fraction, _FIVE_PRIME)
    three = _mean(region_d6_fraction, _THREE_PRIME)
    fp_mp = round(five - three, 3) if (five is not None and three is not None) else None
    present = [v for v in region_d6_fraction.values() if v is not None]
    gene_med = statistics.median(present) if present else None
    ex9 = region_d6_fraction.get("downstream_exon9")
    ex9d = round(gene_med - ex9, 3) if (gene_med is not None and ex9 is not None) else None
    feats = {"five_prime_minus_three_prime": fp_mp, "exon9_tip_dip": ex9d, "n_callable": n_callable}

    if fp_mp is not None and fp_mp >= directional:
        call, star = "*68", "CYP2D6-CYP2D7 hybrid; 5' CYP2D6 / 3' CYP2D7 (intron-1 breakpoint)"
    elif fp_mp is not None and fp_mp <= -directional:
        call, star = "*13", "CYP2D7-CYP2D6 hybrid; 5' CYP2D7 / 3' CYP2D6 (opposite direction) — UNPOWERED (n=1)"
    elif ex9d is not None and ex9d >= exon9_dip:
        call, star = "*36", "CYP2D6-CYP2D7 hybrid; exon-9 CYP2D7 conversion (subtle; some conversions missed)"
    else:
        return {"call": "hybrid_present_identity_unresolved", "confidence": "low",
                "note": "a hybrid signal is present (from the CYP2D7-depth detector) but the PSV profile is "
                        "too weak/ambiguous to name the allele (e.g. a subtle *36 exon-9 conversion)",
                "features": feats}
    return {"call": call, "resolved": True, "confidence": "high", "description": star,
            "note": ("read-level PSV D6-fraction profile -> hybrid identity (high-specificity; spec 1.0 in "
                     "validation). Reference tool: Cyrius. NOT a clinical tool."),
            "features": feats}
