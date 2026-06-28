"""Standing report card for the AMR Portal INDEPENDENT validation (the durable, visible surface).

Read-only roll-up of `wiki/amr_portal_independent_scores.json` (produced by
`scripts/amr_portal_score_independent.py`) into `wiki/amr_portal_independent_report_card.{md,json}`. Exit 0
always — a REPORT, not a gate. NAMESPACE-SEPARATE from the frozen NCBI-PD report card, the HIV card, and the
external-cohort card (the shared-key silent-overwrite lesson): a distinct file, distinct schema.

Per-cell tier (honest, NO aggregate headline):
  SCORED_INDEPENDENT  — powered (>=10 R and >=10 S) + sens/spec present.
  UNDERPOWERED        — a class < 10 (or 0) — reported with whatever CI exists, never hidden.
Rule routing is shown per cell (DRUG_RULE default for E. coli/Shigella; the OPT-IN calibrated registry for
Klebsiella/Salmonella cipro — now independently validated by this very card).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
SCORES = REPO / "wiki" / "amr_portal_independent_scores.json"
OUT_MD = REPO / "wiki" / "amr_portal_independent_report_card.md"
OUT_JSON = REPO / "wiki" / "amr_portal_independent_report_card.json"

def _calibrated_keys() -> set[str]:
    """The set of `<registry_organism>|<drug>` keys whose verdict is CALIBRATED in the frozen rules.
    Data-driven (NOT a hardcoded cell list, which drifts as cells are added — 2026-06-28: Campylobacter
    cipro was mislabeled drug_rule_default by the old hand-list)."""
    import json as _json
    rules = _json.loads((REPO / "dna_decode" / "data" / "calibrated_amr_rules.json").read_text(encoding="utf-8"))
    rules = rules.get("rules", rules)
    return {k for k, v in rules.items() if isinstance(v, dict) and v.get("verdict") == "CALIBRATED"}


def _is_calibrated(org: str, drug: str, calibrated_keys: set[str]) -> bool:
    """A cell is calibrated-routed iff its registry organism|drug is a CALIBRATED key. The scorer's
    RULE_ORGANISM maps the AMR-Portal organism name -> the call_resistance registry organism."""
    from scripts.amr_portal_score_independent import RULE_ORGANISM
    reg = RULE_ORGANISM.get(org)
    return bool(reg) and f"{reg}|{drug}" in calibrated_keys


def _fmt(x):
    return f"{x:.3f}" if isinstance(x, (int, float)) else "—"


def _ci(c):
    return f"[{c[0]:.3f}, {c[1]:.3f}]" if isinstance(c, (list, tuple)) and len(c) == 2 else "—"


def build(scores: dict) -> dict:
    rows = []
    calibrated_keys = _calibrated_keys()
    for key, r in scores.items():
        org, drug = r["organism"], r["drug"]
        powered = bool(r.get("powered"))
        tier = "SCORED_INDEPENDENT" if powered else "UNDERPOWERED"
        routing = "calibrated_registry" if _is_calibrated(org, drug, calibrated_keys) else "drug_rule_default"
        rows.append({"organism": org, "drug": drug, "tier": tier, "routing": routing,
                     "n_R": r["n_R"], "n_S": r["n_S"], "sens": r["sens"], "spec": r["spec"],
                     "accuracy": r["accuracy"], "sens_ci95": r.get("sens_wilson95"),
                     "spec_ci95": r.get("spec_wilson95"), "n_indeterminate": r.get("n_indeterminate", 0)})
    rows.sort(key=lambda x: (x["organism"], x["drug"]))
    return {"schema": "amr-portal-independent-report-card-v1",
            "status_field": "PROVENANCE_DISJOINT_INDEPENDENT_ACCESSION_LEVEL",
            "n_cells": len(rows),
            "n_scored_independent": sum(r["tier"] == "SCORED_INDEPENDENT" for r in rows),
            "n_underpowered": sum(r["tier"] == "UNDERPOWERED" for r in rows),
            "cells": rows}


def render_md(card: dict) -> str:
    L = ["# AMR Portal — INDEPENDENT validation report card",
         "",
         f"Standing roll-up of the frozen deterministic decoder scored on the EBI AMR Portal (CABBAGE) "
         f"provenance-disjoint, measured-AST isolates. **{card['n_scored_independent']} SCORED_INDEPENDENT** "
         f"+ {card['n_underpowered']} UNDERPOWERED of {card['n_cells']} cells. NO aggregate headline (per-cell "
         f"truth only).",
         "",
         "## Honesty rails",
         "- **Independent at the ACCESSION level (upper bound).** Disjoint vs CRyPTIC + our tuning cohorts by "
         "BioSample/ERS/GCA; BioSample cross-archive resolution would only TIGHTEN it.",
         "- **Genotype = the AMR Portal's own AMRFinderPlus run** (different operator → more independent; "
         "AMRFinder-version a named caveat). **Phenotype = wet-lab MIC/disk** (non-circular).",
         "- **Rule applied UNCHANGED**; the frozen surface (`amr_rules.py` + `calibrated_amr_rules.json`) is "
         "byte-unchanged. NAMESPACE-SEPARATE from the NCBI-PD / HIV / external-cohort cards.",
         "- **Calibrated registry now independently validated:** `Salmonella|ciprofloxacin` (broad) + "
         "`Klebsiella|ciprofloxacin` (qrdr_point + oqxAB-exclusion) + `Campylobacter|ciprofloxacin` "
         "(qrdr_point; added 2026-06-28, C. jejuni acc 0.981 / C. coli 0.995) — OPT-IN configs whose "
         "provenance asked for an independent cohort; this card IS that cohort. (Promoting them to DEFAULT "
         "mutates the sha-pinned frozen file → a deliberate ratify-first freeze-amendment, NOT done here.)",
         "",
         "## Cells (per-organism, measured AST, provenance-disjoint)",
         "| Organism | Drug | Tier | Routing | nR / nS | sens (95% CI) | spec (95% CI) | acc |",
         "|---|---|---|---|---|---|---|---|"]
    for r in card["cells"]:
        L.append(f"| {r['organism']} | {r['drug']} | {r['tier']} | {r['routing']} | {r['n_R']}/{r['n_S']} | "
                 f"{_fmt(r['sens'])} {_ci(r['sens_ci95'])} | {_fmt(r['spec'])} {_ci(r['spec_ci95'])} | "
                 f"{_fmt(r['accuracy'])} |")
    L += ["",
          "## Provenance",
          "Scores `wiki/amr_portal_independent_scores.json` (`scripts/amr_portal_score_independent.py`, frozen "
          "rule). Feasibility `wiki/amr_portal_feasibility_result_2026-06-23.md`; validation memo "
          "`wiki/amr_portal_independent_validation_2026-06-23.md`. Rebuild: "
          "`uv run python scripts/build_amr_portal_report_card.py`."]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    if not SCORES.exists():
        print(f"ERROR: scores not found at {SCORES} (run scripts/amr_portal_score_independent.py first)",
              file=sys.stderr)
        return 2
    card = build(json.loads(SCORES.read_text(encoding="utf-8")))
    OUT_JSON.write_text(json.dumps(card, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render_md(card), encoding="utf-8")
    print(f"[amr-portal-card] {card['n_scored_independent']} SCORED_INDEPENDENT / {card['n_underpowered']} "
          f"UNDERPOWERED of {card['n_cells']} cells -> {OUT_MD.name}, {OUT_JSON.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
