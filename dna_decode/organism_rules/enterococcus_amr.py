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

# Glycopeptide van clusters. The RESISTANCE DETERMINANT is the D-Ala-D-X ligase — the BARE `van<letter>`
# symbol (vanA, vanB, vanD, vanM, ...); the accessory genes carry a cluster suffix (vanH-A / vanX-A /
# vanR-A / vanS-A / vanY-A / vanZ-A). Match the ligase (no suffix).
_VAN_LIGASE_RE = re.compile(r"^van([A-N])$")
# Ligase cluster-type -> teicoplanin phenotype. vanA/vanD/vanM = high-level, teico-R; vanB/vanC/vanG/vanN =
# teico-SUSCEPTIBLE (the classic vanA-vs-vanB glycopeptide split). vanE is teico-S; vanL rare, treat teico-S.
_TEICO_R_LIGASES = frozenset({"A", "D", "M"})
# Regulator genes of the two-component induction system (any cluster suffix). Their ABSENCE alongside
# present structural van genes = a vancomycin-VARIABLE enterococcus (VVE): operon present but non-expressed
# -> phenotypically susceptible. Disclosed (not hard-coded into the call) — the genotype is resistance-capable.
_VAN_REGULATOR_RE = re.compile(r"^van[RS]-")
_VAN_STRUCTURAL_RE = re.compile(r"^van[HXY]-")


def _van_ligases(symbols: list[str]) -> list[str]:
    return [s.strip() for s in symbols if _VAN_LIGASE_RE.match((s or "").strip())]


def _van_expression_note(symbols: list[str]) -> str | None:
    """Flag potential non-expression (VVE): structural van genes present but regulator (vanR/vanS) absent."""
    syms = [(s or "").strip() for s in symbols]
    has_struct = any(_VAN_STRUCTURAL_RE.match(s) for s in syms)
    has_reg = any(_VAN_REGULATOR_RE.match(s) for s in syms)
    if has_struct and not has_reg:
        return ("POSSIBLE_NON_EXPRESSION_VVE: structural van genes present but vanR/vanS regulator absent -> "
                "operon may be non-induced (vancomycin-variable enterococcus); genotype is resistance-CAPABLE")
    return None


def call_efm_vancomycin(symbols: list[str]) -> dict:
    """Any acquired van ligase (vanA/B/D/M/...) -> R (standard determinant convention). VVE disclosed."""
    ligases = _van_ligases(symbols)
    return {"prediction": "R" if ligases else "S", "matched_van_ligases": ligases,
            "expression_note": _van_expression_note(symbols),
            "rule": "any acquired van ligase (vanA/vanB/vanD/vanM/...) -> R; VVE non-expression disclosed",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}


def call_efm_teicoplanin(symbols: list[str]) -> dict:
    """vanA/vanD/vanM ligase -> teico-R; vanB/vanC/vanG/vanN only -> teico-S (vanA-vs-vanB split)."""
    ligases = _van_ligases(symbols)
    teico_r = [g for g in ligases if _VAN_LIGASE_RE.match(g).group(1) in _TEICO_R_LIGASES]
    return {"prediction": "R" if teico_r else "S", "matched_teico_r_ligases": teico_r,
            "all_van_ligases": ligases, "expression_note": _van_expression_note(symbols),
            "rule": "vanA/vanD/vanM -> teicoplanin-R; vanB/vanC/vanG/vanN alone -> teicoplanin-S",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}


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
