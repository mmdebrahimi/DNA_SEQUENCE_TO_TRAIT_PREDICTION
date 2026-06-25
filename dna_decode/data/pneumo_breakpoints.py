"""S. pneumoniae β-lactam MIC breakpoints (NON-FROZEN) — the breakpoint-discipline foundation.

The pneumococcus AMR go/no-go found β-lactam R/S is **breakpoint-AMBIGUOUS**: the SAME MIC gives different
R/S under meningitis vs non-meningitis (IV) vs oral breakpoints (penicillin acc 0.969 meningitis → 0.881
non-meningitis on the same genomes). So a pneumo β-lactam cell MUST carry the breakpoint CONTEXT explicitly;
it can never collapse to one R/S. This module is that context set (CLSI M100 S. pneumoniae), kept OUT of the
frozen `mic_tiers.DRUG_BREAKPOINTS` (which is E. coli) so the frozen surface is untouched.

Each entry: drug -> {context: {"S": <=S, "R": >=R}} in µg/mL. `classify(drug, context, mic)` -> R/I/S.
"""
from __future__ import annotations

# CLSI M100 S. pneumoniae β-lactam breakpoints (µg/mL). meningitis/nonmeningitis = parenteral (IV);
# oral = oral penicillin V. I (intermediate) is the band strictly between S and R.
PNEUMO_BETALACTAM_BREAKPOINTS: dict[str, dict[str, dict[str, float]]] = {
    "penicillin": {
        "meningitis":     {"S": 0.06, "R": 0.12},
        "non_meningitis": {"S": 2.0,  "R": 8.0},    # I = 4
        "oral":           {"S": 0.06, "R": 2.0},    # penicillin V (oral)
    },
    "ceftriaxone": {
        "meningitis":     {"S": 0.5, "R": 2.0},
        "non_meningitis": {"S": 1.0, "R": 4.0},
    },
    "cefotaxime": {
        "meningitis":     {"S": 0.5, "R": 2.0},
        "non_meningitis": {"S": 1.0, "R": 4.0},
    },
    "meropenem": {
        "meningitis":     {"S": 0.25, "R": 1.0},
        "non_meningitis": {"S": 0.25, "R": 1.0},
    },
    "cefuroxime": {  # parenteral
        "non_meningitis": {"S": 1.0, "R": 2.0},
    },
}

CONTEXTS = ("meningitis", "non_meningitis", "oral")


def breakpoints_for(drug: str, context: str) -> dict[str, float] | None:
    return PNEUMO_BETALACTAM_BREAKPOINTS.get((drug or "").strip().lower(), {}).get(context)


def classify(drug: str, context: str, mic: float | None) -> str | None:
    """R / I / S for `mic` under the (drug, context) breakpoint, or None if mic/context unknown."""
    bp = breakpoints_for(drug, context)
    if bp is None or mic is None:
        return None
    return "S" if mic <= bp["S"] else "R" if mic >= bp["R"] else "I"
