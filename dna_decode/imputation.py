"""Fail-closed LD-imputation pre-processor for the deterministic decoder (Phase 4 — productionizes V2).

Impute an UNCALLABLE determinant genotype from a linked tag SNP using a FROZEN, data-derived LD map, so the
deterministic rule can call otherwise-ABSTAINED samples. Two invariants make this safe to put in front of the
frozen decoder:

  1. **Fail-closed.** Impute ONLY when the tag maps to its majority target at purity >= `min_purity`
     (default 0.90); otherwise return ABSTAIN. A wrong impute is worse than an honest ABSTAIN.
  2. **Provenance-tagged.** An imputed call is NEVER confused with a directly-typed one — every result
     carries `provenance` ("direct" / "imputed:<tag>=<gt>@<purity>" / "abstain:<reason>").

The map is a COMMITTED, data-derived artifact (`data/imputation/*.json`, frozen by
`scripts/impute_determinant_abstain.py --dump-map`), NOT fabricated. Validated 2026-07-04
(`wiki/impute_abstain_abo_result_2026-07-04.md`: ABO O-deletion imputed at 98.9%). This module NEVER touches
the frozen AMR surface — it is a pure input-completion layer in front of the deterministic call.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

_UNCALLABLE = {"", "--", "-", "NA", "N/A", None}
DEFAULT_MIN_PURITY = 0.90                      # fail-closed threshold; below this -> ABSTAIN

_REPO = Path(__file__).resolve().parent.parent
# registry: determinant target rsid -> committed frozen LD-map path
IMPUTATION_MAPS: dict[str, Path] = {
    "rs8176719": _REPO / "data" / "imputation" / "abo_rs8176719_from_rs657152.json",   # ABO O-status
}


@dataclass(frozen=True)
class Imputation:
    genotype: str | None
    provenance: str
    confidence: float | None


@dataclass(frozen=True)
class LdImputer:
    target: str
    tag: str
    table: dict                                # {tag_gt: {"majority": gt, "purity": float, "n": int}}
    min_purity: float = DEFAULT_MIN_PURITY

    @classmethod
    def from_json(cls, path: str | Path, min_purity: float = DEFAULT_MIN_PURITY) -> "LdImputer":
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(target=d["target"], tag=d["tag"], table=d["map"], min_purity=min_purity)

    @classmethod
    def for_target(cls, target_rsid: str, min_purity: float = DEFAULT_MIN_PURITY) -> "LdImputer | None":
        p = IMPUTATION_MAPS.get(target_rsid)
        return cls.from_json(p, min_purity) if (p and p.exists()) else None

    def impute(self, target_gt, tag_gt) -> Imputation:
        """Fail-closed: pass through a called target; else impute from the tag only at purity>=min_purity."""
        if target_gt not in _UNCALLABLE and str(target_gt).strip() not in _UNCALLABLE:
            return Imputation(target_gt, "direct", 1.0)               # already typed -> no imputation
        if tag_gt in _UNCALLABLE or str(tag_gt).strip() in _UNCALLABLE:
            return Imputation(None, "abstain:no-tag", None)
        entry = self.table.get(str(tag_gt).strip())
        if entry is None:
            return Imputation(None, f"abstain:tag-genotype-unseen({tag_gt})", None)
        pur = float(entry["purity"])
        if pur < self.min_purity:
            return Imputation(None, f"abstain:low-purity({pur}<{self.min_purity})", pur)
        return Imputation(entry["majority"], f"imputed:{self.tag}={tag_gt}@{pur}", pur)


def call_with_imputation(call_fn: Callable[[str], str], target_gt, tag_gt, imputer: LdImputer,
                         abstain_value: str = "INDETERMINATE") -> dict:
    """Impute-then-call: run the FROZEN deterministic `call_fn` on the (possibly imputed) genotype.
    Returns {call, provenance, confidence}. On ABSTAIN, `call` is `abstain_value` (the decoder still abstains,
    honestly, when no confident tag) — never a guessed call."""
    imp = imputer.impute(target_gt, tag_gt)
    if imp.genotype is None:
        return {"call": abstain_value, "provenance": imp.provenance, "confidence": imp.confidence}
    return {"call": call_fn(imp.genotype), "provenance": imp.provenance, "confidence": imp.confidence}
