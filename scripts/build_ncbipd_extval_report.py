"""Roll-up report card for the NCBI-PD no-compute external-validation track (NON-FROZEN organism cells).

Reads the per-organism `wiki/*_ncbipd_extval_*.json` artifacts + renders one honest roll-up:
`wiki/ncbipd_extval_report_card.{md,json}`. Every cell keeps its own tier (raw + lineage-collapsed
sens/spec, verdict); NO aggregate headline (per project discipline).

NAMESPACE NOTE: this is a SEPARATE report from the frozen `external_validation_report_card` (the Oxford
arm) — these are NON-FROZEN organism cells (gono/campy/staph/pneumo) validated on NCBI-PD's own AMRFinder
calls, a different track. Emitting into the frozen card would risk the shared-key overwrite trap.

  uv run python scripts/build_ncbipd_extval_report.py
"""
from __future__ import annotations

import glob
import json
import re
from datetime import date as _date
from pathlib import Path

# organism-artifact glob (the generic scorer writes wiki/<organism>_ncbipd_extval_<date>.json)
ARTIFACTS = sorted(glob.glob("wiki/*_ncbipd_extval_*.json"))


def _fmt(x):
    return "—" if x is None else (f"{x:.3f}" if isinstance(x, (int, float)) else str(x))


def _load_blindness() -> dict:
    """(organism, drug) -> invisible_fraction from the determinant-blindness atlas, if built."""
    p = Path("wiki/determinant_blindness_atlas.json")
    if not p.exists():
        return {}
    atlas = json.loads(p.read_text(encoding="utf-8"))
    return {(c["organism"], c["drug"]): c["invisible_fraction"] for c in atlas.get("cells", [])}


def main() -> int:
    blindness = _load_blindness()
    rows = []
    for path in ARTIFACTS:
        art = json.loads(Path(path).read_text(encoding="utf-8"))
        org = art.get("organism", "?")
        for drug, r in (art.get("results") or {}).items():
            b = r.get("binary") or {}
            lin = r.get("lineage_collapsed") or {}
            rows.append({
                "organism": org, "drug": drug, "n": b.get("n_scored"),
                "n_R": r.get("n_R"), "n_S": r.get("n_S"),
                "raw_sens": b.get("sens"), "raw_spec": b.get("spec"),
                "lin_sens": lin.get("sens"), "lin_spec": lin.get("spec"),
                "lin_discordant": lin.get("n_discordant"),
                "invisible_fraction": blindness.get((org, drug)),
                "verdict": r.get("headline"),
            })
    rows.sort(key=lambda x: (x["organism"], x["drug"]))
    endorsed = [r for r in rows if r["verdict"] == "SCORED_ENDORSED"]

    lines = [
        "# NCBI-PD no-compute external-validation report card",
        "",
        f"_Generated {_date.today().isoformat()} from {len(ARTIFACTS)} organism artifacts. Roll-up of the "
        "reusable NCBI-PD external-validation substrate (`scripts/score_ncbipd_extval.py`): a pure metadata "
        "join + score of NON-FROZEN organism cells against NCBI Pathogen Detection's OWN published "
        "AMRFinderPlus calls, with no-compute lineage-collapse via NCBI-PD SNP clusters._",
        "",
        "**Honest tier (no aggregate headline):** provenance-disjoint (frozen accessions excluded) but NOT "
        "methodology-independent (same AMRFinderPlus + same cell). RAW sens/spec is clonality-inflated; the "
        "**lineage-collapsed** columns (one vote per NCBI-PD SNP cluster) are the honest number. A cell "
        "`SCORED_ENDORSED` only if powered (>=5/class), non-degenerate, spec >= 0.85, sens >= 0.5.",
        "",
        f"**{len(endorsed)} of {len(rows)} cells SCORED_ENDORSED** across "
        f"{len({r['organism'] for r in rows})} organisms.",
        "",
        "| organism | drug | n (R/S) | raw sens/spec | **lineage sens/spec** | disc. | blind. | verdict |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        mark = "✅ " if r["verdict"] == "SCORED_ENDORSED" else ""
        lines.append(
            f"| {r['organism']} | {r['drug']} | {r['n']} ({r['n_R']}R/{r['n_S']}S) | "
            f"{_fmt(r['raw_sens'])}/{_fmt(r['raw_spec'])} | "
            f"**{_fmt(r['lin_sens'])}/{_fmt(r['lin_spec'])}** | {_fmt(r['lin_discordant'])} | "
            f"{_fmt(r['invisible_fraction'])} | {mark}{r['verdict']} |")
    lines += [
        "",
        "## Notes",
        "- **`blind.` = determinant-invisible fraction** (`wiki/determinant_blindness_atlas.md`): of the "
        "measured-R isolates, the fraction the cell calls non-R because no rule-firing catalog determinant is "
        "present. A cell can be SCORED_ENDORSED (spec holds) yet highly blind (e.g. gono tetracycline spec "
        "1.0 but 0.68 invisible — high-level TRNG only). It is a DESCRIPTIVE honesty column, not a metric the "
        "endorsement gate uses.",
        "- **Lineage-collapse is no-compute** — NCBI-PD publishes per-isolate SNP clusters "
        "(`<PDG>.reference_target.all_isolates.tsv` → `PDS_acc`), collapsed via "
        "`clonality.cluster_weighted_confusion` (no Mash/Docker). Every endorsed cell HOLDS at the lineage "
        "level → the rules decode mechanism, not clonal structure.",
        "- **DEGENERATE guard**: a cell predicting all-one-class (e.g. gono azithromycin all-S) is never "
        "endorsed even at spec/sens 1.0.",
        "- **NON-FROZEN cells**; the frozen decoder surface is byte-unchanged throughout (`verify_lock` OK).",
        "- Per-cell detail: the `wiki/*_ncbipd_extval_*.md` result docs.",
    ]
    Path("wiki/ncbipd_extval_report_card.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    Path("wiki/ncbipd_extval_report_card.json").write_text(
        json.dumps({"_schema": "ncbipd-extval-report-card-v1", "date": _date.today().isoformat(),
                    "n_endorsed": len(endorsed), "n_cells": len(rows), "cells": rows}, indent=2),
        encoding="utf-8")
    print(f"{len(endorsed)}/{len(rows)} cells SCORED_ENDORSED across "
          f"{len({r['organism'] for r in rows})} organisms")
    print("wrote wiki/ncbipd_extval_report_card.{md,json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
