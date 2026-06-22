"""HIV viral-decoder trust surface — a read-only roll-up of the HIV validation artifacts (Wave B Rec 3).

DESIGN DECISION (best-judgment override of "fold HIV into the bacterial report card"): the bacterial
report card (`build_validation_report_card.py`) is the NCBI-PD PROVENANCE-DISJOINT external stress test;
the HIV cell is validated IN-DISTRIBUTION against the Stanford HIVDB PhenoSense dataset. Those are DIFFERENT
rigour modalities, and the project's namespace-separation discipline (the Oxford external-validation arm)
says: do NOT conflate them in one grid. So HIV gets its OWN first-class trust surface here, explicitly
labelled as the in-distribution HIVDB-PhenoSense modality — NOT the bacterial provenance-disjoint one.

Read-only: exit 0 always (a report, not a gate). Rolls up whatever HIV validation JSONs exist in wiki/.
"""
from __future__ import annotations

import glob
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WIKI = REPO / "wiki"


def _latest(prefix: str) -> dict | None:
    files = sorted(glob.glob(str(WIKI / f"{prefix}*.json")))
    if not files:
        return None
    try:
        return json.loads(Path(files[-1]).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def build() -> dict:
    nnrti_val = _latest("hiv_nnrti_v0_validation_")
    nnrti_base = _latest("hiv_nnrti_baseline_vs_ols_")
    nrti_val = _latest("hiv_nrti_v0_validation_")
    nrti_mut = _latest("hiv_nrti_mutant_catalog_v0_1_")
    nrti_sub = _latest("hiv_nrti_within_subtype_")

    cells = []
    # NNRTI: AUC (call separates fold) + catalog-vs-OLS balacc delta
    if nnrti_val:
        base = (nnrti_base or {}).get("per_drug", {})
        for drug, m in nnrti_val.get("per_drug", {}).items():
            b = base.get(drug, {})
            cells.append({
                "drug_class": "NNRTI", "drug": drug, "catalog": "mutant-specific (Stanford majors)",
                "n": m.get("n_isolates_with_fold"),
                "auc_call_separates_fold": m.get("auc_call_separates_fold"),
                "catalog_balacc": (b.get("catalog") or {}).get("balanced_accuracy"),
                "ols_baseline_balacc": (b.get("ols_baseline") or {}).get("balanced_accuracy"),
                "delta_ols_minus_catalog": b.get("delta_ols_minus_catalog_balacc"),
                "subtype_transfer": None, "v0_1_mutant_gain": None,
            })
    # NRTI: catalog vs OLS + mutant-specific v0.1 gain + subtype transfer
    if nrti_val:
        mut = (nrti_mut or {}).get("per_drug", {})
        sub = (nrti_sub or {}).get("per_drug", {})
        for drug, m in nrti_val.get("per_drug", {}).items():
            if "catalog" not in m:
                continue
            mm = mut.get(drug, {})
            ss = sub.get(drug, {})
            nonb = (ss.get("non-B") or {})
            cells.append({
                "drug_class": "NRTI", "drug": drug, "catalog": "position-based v0 (+ mutant-specific v0.1)",
                "n": m.get("n_isolates"),
                "auc_call_separates_fold": (m.get("catalog") or {}).get("auc_call_separates_fold"),
                "catalog_balacc": (m.get("catalog") or {}).get("balanced_accuracy"),
                "ols_baseline_balacc": (m.get("ols_baseline") or {}).get("balanced_accuracy"),
                "delta_ols_minus_catalog": m.get("delta_ols_minus_catalog_balacc"),
                "v0_1_mutant_gain": mm.get("balacc_gain_v0_1_minus_v0"),
                "subtype_transfer": (f"non-B balacc {nonb.get('balanced_accuracy')} (n={nonb.get('n')})"
                                     if nonb.get("balanced_accuracy") is not None else "under-powered"),
            })
    return {
        "artifact": "hiv_decoder_report_card", "schema": "hiv-report-card-v0",
        "modality": "IN-DISTRIBUTION validation vs Stanford HIVDB PhenoSense fold-change (independent "
                    "wet-lab IC50). DISTINCT from the bacterial NCBI-PD provenance-disjoint report card — "
                    "NOT conflated (a different, more external rigour modality).",
        "label_independence": "PhenoSense fold-change is NOT HIVDB's own Sierra interpretation (non-circular)",
        "underlying_tool_baseline": "Stanford DRMcv.R OLS, reimplemented (sklearn)",
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use",
        "n_cells": len(cells),
        "honest_caveats": [
            "in-distribution (HIVDB), NOT provenance-disjoint -> a lower external-rigour bar than the bacterial card",
            "NNRTI = mutant-specific (excellent on 1st-gen EFV/NVP); NRTI v0 = position-based (over-calls, "
            "fixed by the deconfounded mutant-specific v0.1 for 5/6 drugs; ddI keeps position-based)",
            "non-B subtype transfer is under-powered (data ~96% subtype B)",
        ],
        "cells": cells,
    }


def render_md(rc: dict, generated: str) -> str:
    lines = [f"# HIV viral-decoder report card ({generated})", ""]
    lines.append(f"**Modality:** {rc['modality']}")
    lines.append(f"**Label independence:** {rc['label_independence']}.")
    lines.append(f"**Underlying-tool baseline:** {rc['underlying_tool_baseline']}.")
    lines.append("")
    lines.append("| Class | Drug | n | AUC (call sep. fold) | catalog balacc | OLS balacc | Δ(OLS−cat) | v0.1 gain | subtype transfer |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for c in rc["cells"]:
        lines.append(f"| {c['drug_class']} | {c['drug']} | {c['n']} | {c['auc_call_separates_fold']} | "
                     f"{c['catalog_balacc']} | {c['ols_baseline_balacc']} | {c['delta_ols_minus_catalog']} | "
                     f"{c.get('v0_1_mutant_gain') if c.get('v0_1_mutant_gain') is not None else '-'} | "
                     f"{c.get('subtype_transfer') or '-'} |")
    lines.append("")
    lines.append("## Honest caveats")
    for cv in rc["honest_caveats"]:
        lines.append(f"- {cv}")
    lines.append("")
    lines.append(f"Citation: {rc['citation']}.")
    return "\n".join(lines)


def main(argv=None) -> int:
    today = _date.today().isoformat()
    rc = build()
    (WIKI / "hiv_decoder_report_card.json").write_text(json.dumps(rc, indent=2), encoding="utf-8")
    (WIKI / "hiv_decoder_report_card.md").write_text(render_md(rc, today), encoding="utf-8")
    print(render_md(rc, today))
    print(f"\n[wrote {WIKI / 'hiv_decoder_report_card.md'} + .json | {rc['n_cells']} cells]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
