"""Campylobacter jejuni/coli tet + gent — curated rules (Tier-3 completion, NON-frozen).

Campylobacter already has the endorsed `Campylobacter|ciprofloxacin` registry rule (qrdr_point). These add:
- **tetracycline** (class B): acquired **tet(O)** ribosomal-protection (+ mosaic tet(O/M/O)) -> R. This is
  DISTINCT from the E. coli efflux tetA set (why the generic default rule was guardrail-blocked here).
- **gentamicin** (class B): TRUE gentamicin enzymes only — aph(2''), aac(3), aac(6')-aph bifunctional. The
  Campylobacter AMINOGLYCOSIDE determinants in the AMR Portal are dominated by **aad9** (spectinomycin/
  streptomycin) + **spw**, which do NOT confer gentamicin resistance and are EXCLUDED (intrinsic-exclusion).
  If no true gent marker is present the rule calls S — expected to under-detect (honest INDETERMINATE).

NON-FROZEN / scorer-local; endorsed only if spec>=0.85 on the AMR Portal.
"""
from __future__ import annotations

DRUG_TET = "tetracycline"
DRUG_GENT = "gentamicin"
_GENT_MARKERS = ("aph(2'')", "aac(3)", "aac(6')-aph", "aac(6')-Ie")


def call_cj_tetracycline(symbols: list[str]) -> dict:
    hits = [s.strip() for s in symbols if (s or "").strip().startswith("tet(O")]
    return {"prediction": "R" if hits else "S", "matched_tetO": hits,
            "rule": "tet(O)-family ribosomal protection -> R",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}


def call_cj_gentamicin(symbols: list[str]) -> dict:
    hits = [s.strip() for s in symbols if any(m in (s or "") for m in _GENT_MARKERS)]
    return {"prediction": "R" if hits else "S", "matched_gent": hits,
            "rule": "true gentamicin enzyme aph(2'')/aac(3)/aac(6')-aph -> R (aad9/spw non-gent excluded)",
            "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local"}
