"""Staphylococcus aureus ciprofloxacin — curated organism-specific rule (Tier-3/4 cell 2, NON-frozen).

Second AMR-Portal Tier-3/4 curated cell (plan `plans/AMR_Portal_Tier34_Curation_Plan_2026-07-01.md`),
mechanism class A (target-site point mutation). S. aureus fluoroquinolone resistance is driven by QRDR
substitutions in BOTH topoisomerase targets (unlike N. gonorrhoeae where gyrA is primary + parC accessory):
grlA/parC (Ser80, Glu84) is the FIRST-step target (parC S80F alone typically already lifts cipro above the
breakpoint), and gyrA (Ser84, Ser85, Glu88) is the SECOND-step / high-level target. So EITHER a gyrA or a
grlA/parC QRDR mutation confers resistance. (AMRFinderPlus reports the S. aureus grlA gene under the `parC`
symbol, e.g. `parC_S80F`.)

NON-FROZEN + scorer-local (same posture as `neisseria_amr` / `tb_*`): NOT in the frozen `amr_rules.py` /
`calibrated_amr_rules.json` deployed surface. Scored on the EBI AMR Portal via
`scripts/staph_cipro_amr_portal_validate.py`; endorsed only if it clears the spec>=0.85 falsifier.
"""
from __future__ import annotations

import re

DRUG = "ciprofloxacin"
ORGANISM = "Staphylococcus aureus"

# QRDR resistance codons — both genes are PRIMARY for S. aureus FQ resistance.
QRDR_CODONS = {"gyrA": frozenset({84, 85, 88}), "parC": frozenset({80, 84})}
_POINT_RE = re.compile(r"^(gyrA|parC)_([A-Z])(\d+)([A-Z*])$")


def call_sa_ciprofloxacin(symbols: list[str]) -> dict:
    """Predict cipro R/S for S. aureus. R iff >=1 QRDR point mutation in gyrA (84/85/88) OR parC/grlA (80/84)."""
    hits = {"gyrA": [], "parC": []}
    for s in symbols:
        m = _POINT_RE.match((s or "").strip())
        if not m:
            continue
        gene, pos = m.group(1), int(m.group(3))
        if pos in QRDR_CODONS[gene]:
            hits[gene].append(s.strip())
    resistant = bool(hits["gyrA"] or hits["parC"])
    return {
        "prediction": "R" if resistant else "S",
        "matched_gyrA_qrdr": hits["gyrA"], "matched_parC_qrdr": hits["parC"],
        "rule": "gyrA (Ser84/Ser85/Glu88) OR grlA/parC (Ser80/Glu84) QRDR point mutation -> R",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


# --- rifampicin (rpoB RRDR) -------------------------------------------------------------------------------
DRUG_RIF = "rifampicin"
_RPOB_RE = re.compile(r"^rpoB_[A-Z]\d+[A-Z*]$")


def call_sa_rifampicin(symbols: list[str]) -> dict:
    """Predict rifampicin R/S for S. aureus from its determinant symbols.

    R iff AMRFinderPlus reports ANY rpoB RRDR point substitution (e.g. rpoB_H481Y, rpoB_L466S, rpoB_A477V).
    Unlike a fixed-codon QRDR list, rif resistance spans many RRDR codons (~464-529) — but AMRFinder's rpoB
    POINT calls are ALL curated resistance mutations (it does not emit benign rpoB polymorphisms as AMR),
    so presence of any is the correct deterministic determinant. Frozen surface untouched (scorer-local)."""
    matched = [s.strip() for s in symbols if _RPOB_RE.match((s or "").strip())]
    return {
        "prediction": "R" if matched else "S",
        "matched_rpoB": matched,
        "rule": "any AMRFinder rpoB RRDR resistance point mutation -> R",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }
