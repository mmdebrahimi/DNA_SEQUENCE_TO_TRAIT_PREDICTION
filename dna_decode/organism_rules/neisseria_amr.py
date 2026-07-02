"""Neisseria gonorrhoeae ciprofloxacin — curated organism-specific determinant rule (Tier-3 cell, NON-frozen).

First AMR-Portal Tier-3/4 curated cell (plan `plans/AMR_Portal_Tier34_Curation_Plan_2026-07-01.md`).
Mechanism class A (target-site point mutation) — genuinely decodable. NG cipro resistance is driven by
**gyrA QRDR** substitutions at Ser91 and Asp95 (the canonical first-/second-step fluoroquinolone-resistance
mutations; Belland 1994, WHO GC surveillance): gyrA S91F is by far the dominant determinant, with S91Y and
D95N/A/G also resistance-conferring. parC (Ser87/Ser88/Glu91) mutations are ACCESSORY — they raise MIC on a
gyrA background but rarely confer resistance alone — so the primary binary rule is gyrA-QRDR-driven.

NON-FROZEN + scorer-local (like `organism_rules/tb_*` + `experimental_drug_rules`): this is NOT in the frozen
`amr_rules.py` / `calibrated_amr_rules.json` deployed surface (those stay byte-pinned). Scored on the EBI AMR
Portal via `scripts/neisseria_cipro_amr_portal_validate.py`; endorsed only if it clears the spec>=0.85
falsifier on the powered provenance-disjoint set.

Input = the isolate's AMRFinderPlus point-mutation symbols (the Portal reports them in `amr_element_symbol`,
e.g. `gyrA_S91F`, with `element_subtype=POINT`).
"""
from __future__ import annotations

import re

DRUG = "ciprofloxacin"
ORGANISM = "Neisseria gonorrhoeae"

# QRDR resistance codons (any substitution at these positions confers/contributes to FQ resistance).
GYRA_QRDR_CODONS = frozenset({91, 95})          # Ser91, Asp95 — primary
PARC_QRDR_CODONS = frozenset({87, 88, 91})      # Ser87, Ser88, Glu91 — accessory (recorded, not required)
_POINT_RE = re.compile(r"^(gyrA|parC)_([A-Z])(\d+)([A-Z*])$")


def _qrdr_hits(symbols: list[str]) -> dict[str, list[str]]:
    """Parse AMRFinder point-mutation symbols -> {'gyrA': [...], 'parC': [...]} of QRDR-codon substitutions."""
    hits = {"gyrA": [], "parC": []}
    for s in symbols:
        m = _POINT_RE.match((s or "").strip())
        if not m:
            continue
        gene, _wt, pos, _sub = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        codons = GYRA_QRDR_CODONS if gene == "gyrA" else PARC_QRDR_CODONS
        if pos in codons:
            hits[gene].append(s.strip())
    return hits


DRUG_TET = "tetracycline"


def call_ng_tetracycline(symbols: list[str]) -> dict:
    """Predict tetracycline R/S for N. gonorrhoeae.

    Two determinants: acquired **tet(M)** (plasmid ribosomal-protection -> HIGH-level R) and the **rpsJ V57M**
    ribosomal S10 point mutation (chromosomal -> LOW-level R, often at/above the CLSI R>=2 breakpoint). Either
    -> R. (mtrR/porB efflux variants contribute to multidrug tolerance but are not the primary tet
    determinant.) Non-frozen / scorer-local; endorsed only if it clears the spec>=0.85 falsifier on the AMR
    Portal (the rpsJ-only isolates are the spec risk if V57M sits below the breakpoint)."""
    tetm = [s.strip() for s in symbols if s and s.strip().startswith("tet(M)")]
    rpsj = [s.strip() for s in symbols if (s or "").strip() == "rpsJ_V57M"]
    # tet(M)-ONLY: rpsJ V57M is common + LOW-level (frequently below the CLSI R>=2 breakpoint), so including
    # it collapses specificity (empirically 0.35 on the AMR Portal). Demoted to accessory context.
    return {
        "prediction": "R" if tetm else "S",
        "matched_tetM": tetm, "accessory_rpsJ_V57M": rpsj,
        "rule": "tet(M) (high-level ribosomal protection) -> R (rpsJ V57M accessory-only: low-level, over-calls)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


def call_ng_ciprofloxacin(symbols: list[str]) -> dict:
    """Predict cipro R/S for N. gonorrhoeae from its determinant symbols.

    R iff >=1 gyrA QRDR (Ser91/Asp95) substitution is present (the deterministic primary rule). parC QRDR
    hits are surfaced as accessory context but do NOT change the binary call on their own."""
    hits = _qrdr_hits(symbols)
    resistant = bool(hits["gyrA"])
    return {
        "prediction": "R" if resistant else "S",
        "matched_gyrA_qrdr": hits["gyrA"],
        "accessory_parC_qrdr": hits["parC"],
        "rule": "gyrA QRDR (Ser91/Asp95) point mutation -> R (parC accessory-only)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }
