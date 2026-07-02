"""Enterococcus faecium AMR — curated organism-specific rules (Tier-3/4 Phase C, NON-frozen).

Gram-positive; the E. coli DRUG_RULE default MUST NOT be sprayed here (intrinsic-gene over-call trap). Three
decodable cells, each with its own determinant biology:

- **ciprofloxacin** (class A): gyrA QRDR (Ser83, Glu87) OR parC QRDR (Ser80, Glu84) point mutation -> R.
- **tetracycline** (class B): any acquired tet ribosomal-protection/efflux gene — tet(M)/tet(L)/tet(S)/tet(O) -> R.
- **gentamicin** (class B, INTRINSIC-EXCLUSION critical): HIGH-level gentamicin resistance requires the
  bifunctional **aph(2'')** enzyme (aac(6')-Ie-aph(2'')-Ia). The chromosomal **aac(6')-Ii** is INTRINSIC to
  E. faecium (low-level, kanamycin — NOT gentamicin) and MUST be excluded, else spec collapses. Rule: any
  symbol carrying `aph(2'')` -> R; aac(6')-Ii alone -> S.

NON-FROZEN / scorer-local; endorsed only if spec>=0.85 on the AMR Portal.
"""
from __future__ import annotations

import re

ORGANISM = "Enterococcus faecium"
_QRDR_RE = re.compile(r"^(gyrA|parC)_([A-Z])(\d+)([A-Z*])$")
_QRDR_CODONS = {"gyrA": frozenset({83, 87}), "parC": frozenset({80, 84})}
_TET_RE = re.compile(r"^tet\(")


def call_efm_ciprofloxacin(symbols: list[str]) -> dict:
    hits = [s.strip() for s in symbols
            if (m := _QRDR_RE.match((s or "").strip())) and int(m.group(3)) in _QRDR_CODONS[m.group(1)]]
    return {"prediction": "R" if hits else "S", "matched_qrdr": hits,
            "rule": "gyrA (Ser83/Glu87) OR parC (Ser80/Glu84) QRDR point mutation -> R",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}


def call_efm_tetracycline(symbols: list[str]) -> dict:
    hits = [s.strip() for s in symbols if _TET_RE.match((s or "").strip())]
    return {"prediction": "R" if hits else "S", "matched_tet": hits,
            "rule": "any acquired tet gene tet(M)/tet(L)/tet(S)/tet(O) -> R",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}


def call_efm_gentamicin(symbols: list[str]) -> dict:
    """HIGH-level gent = bifunctional aph(2''). Intrinsic aac(6')-Ii (kanamycin, low-level) EXCLUDED."""
    hits = [s.strip() for s in symbols if "aph(2'')" in (s or "")]
    return {"prediction": "R" if hits else "S", "matched_aph2": hits,
            "rule": "aph(2'') bifunctional (high-level) -> R; intrinsic aac(6')-Ii excluded",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}
