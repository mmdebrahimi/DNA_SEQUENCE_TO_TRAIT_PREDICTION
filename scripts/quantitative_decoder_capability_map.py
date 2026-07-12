"""Quantitative-decoder capability map — a standing rollup of every CALIBRATED cell (2026-07-12).

The deterministic decoder emits R/S calls; the QUANTITATIVE layer (Family B) emits coverage-valid prediction
INTERVALS. This rolls up every calibrated cell into one map so "what can the decoder put a calibrated number
on, and how tight" is answerable at a glance:

  - HIV drug-resistance fold-change (Stanford PhenoSense; split-conformal on log10 fold) — 4 classes x drugs.
  - CRyPTIC TB BMD-MIC (wet-lab; censoring-aware split-conformal on log2 MIC) — the drug panel.

Read-only roll-up (exit 0 always; a REPORT, not a gate). Reads the latest calibration artifacts under wiki/;
prefers the TB PANEL artifact if present, else the RIF+INH base. Frozen decoder surface untouched.
"""
from __future__ import annotations

import argparse
import glob
import json
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
INFORMATIVE_R2 = 0.05   # cell is INFORMATIVELY calibrated iff its model beats the marginal label (R2 > 5%)


def _informative(r2):
    return bool(r2 is not None and r2 > INFORMATIVE_R2)


def _latest(*patterns):
    files = []
    for p in patterns:
        files += glob.glob(str(REPO / "wiki" / p))
    return sorted(files)[-1] if files else None


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8")) if path else None


def rows_from_hiv(d):
    """HIV per_class -> flat rows. cover key = cover_90, fold = fold_factor_90, calibrated_90."""
    rows = []
    if not d:
        return rows
    for cls, block in d.get("per_class", {}).items():
        for drug, m in block.get("per_drug", {}).items():
            if not m.get("powered"):
                rows.append({"pathogen": "HIV-1", "unit": f"{cls}/{drug}", "label": "PhenoSense fold-change",
                             "n": m.get("n", 0), "r2": None, "cover_90": None, "fold_factor": None,
                             "calibrated": None, "status": "under-powered"})
                continue
            rows.append({"pathogen": "HIV-1", "unit": f"{cls}/{drug}", "label": "PhenoSense fold-change",
                         "n": m["n"], "r2": m.get("r2_oof"), "cover_90": m.get("cover_90"),
                         "fold_factor": m.get("fold_factor_90"), "calibrated": bool(m.get("calibrated_90")),
                         "informative": _informative(m.get("r2_oof")),
                         "status": "calibrated" if m.get("calibrated_90") else "powered-not-calibrated"})
    return rows


def rows_from_tb(d):
    """TB per_drug -> flat rows. cover key = cover_resolved_90, fold = interval_fold_factor."""
    rows = []
    if not d:
        return rows
    for drug, m in d.get("per_drug", {}).items():
        if not m.get("powered"):
            rows.append({"pathogen": "M. tuberculosis", "unit": drug, "label": "CRyPTIC BMD-MIC",
                         "n": m.get("n", 0), "r2": None, "cover_90": None, "fold_factor": None,
                         "calibrated": None, "status": f"under-powered ({m.get('note','')})".strip()})
            continue
        rows.append({"pathogen": "M. tuberculosis", "unit": drug, "label": "CRyPTIC BMD-MIC",
                     "n": m["n"], "r2": m.get("r2_oof_resolved"), "cover_90": m.get("cover_resolved_90"),
                     "fold_factor": m.get("interval_fold_factor"), "calibrated": bool(m.get("calibrated")),
                     "informative": (m["informative"] if "informative" in m
                                     else _informative(m.get("r2_oof_resolved"))),
                     "status": "calibrated" if m.get("calibrated") else "powered-not-calibrated"})
    return rows


def build(hiv_path=None, tb_path=None):
    hiv_path = hiv_path or _latest("hiv_quantitative_calibration_*.json")
    # prefer the panel artifact if it exists
    tb_path = tb_path or _latest("tb_mic_calibration_panel_*.json") or _latest("tb_mic_calibration_*.json")
    hiv, tb = _load(hiv_path), _load(tb_path)
    rows = rows_from_hiv(hiv) + rows_from_tb(tb)
    calibrated = [r for r in rows if r["status"] == "calibrated"]
    informative = [r for r in calibrated if r.get("informative")]
    powered = [r for r in rows if r["calibrated"] is not None]
    covers = [r["cover_90"] for r in calibrated if r["cover_90"] is not None]
    return {
        "artifact": "quantitative_decoder_capability_map", "schema": "qdcm-v1",
        "sources": {"hiv": Path(hiv_path).name if hiv_path else None,
                    "tb": Path(tb_path).name if tb_path else None},
        "summary": {
            "n_cells_total": len(rows), "n_powered": len(powered), "n_calibrated": len(calibrated),
            "n_informative": len(informative),
            "n_coverage_valid_only": len(calibrated) - len(informative),
            "n_under_powered": sum(1 for r in rows if r["calibrated"] is None),
            "coverage_range": [round(min(covers), 4), round(max(covers), 4)] if covers else None,
            "pathogens": sorted({r["pathogen"] for r in rows}),
        },
        "honest_caveats": [
            "CALIBRATED != INFORMATIVE: conformal coverage holds even for a useless model (the interval widens to "
            "the marginal label spread). `informative` (R2>0.05) flags the cells whose determinant model actually "
            "beats the mean — only those intervals carry genotype signal. The coverage-valid-only cells (some TB "
            "second-line drugs with few determinants + rare resistance) hit 0.90 coverage trivially.",
            "Every cell is IN-DISTRIBUTION vs its own knowledge base (HIV: Stanford catalog features; TB: WHO "
            "catalogue determinants) — the interval is a coverage guarantee, NOT an independent-validation claim.",
            "Split-conformal gives MARGINAL (population-level) coverage, not per-genotype. TB MIC coverage is on "
            "the RESOLVED (uncensored) subset; censored isolates are scored only for consistency.",
            "This map rolls up the QUANTITATIVE (interval) layer only — the deterministic R/S decoder + its "
            "provenance-disjoint report card are separate surfaces.",
        ],
        "rows": rows,
    }


def render_md(res, generated):
    s = res["summary"]
    cov = res["summary"]["coverage_range"]
    L = [f"# Quantitative-decoder capability map ({generated})", "",
         f"**{s['n_informative']}/{s['n_powered']} powered cells INFORMATIVELY calibrated** "
         f"({s['n_calibrated']} coverage-valid, of which {s['n_coverage_valid_only']} are coverage-valid-only "
         f"— model no better than the marginal) at 90% target coverage across {', '.join(s['pathogens'])}"
         + (f"; calibrated-cell coverage spans {cov[0]}–{cov[1]}." if cov else "."),
         f"({s['n_under_powered']} cells under-powered; {s['n_cells_total']} cells total.)", "",
         f"Sources: HIV `{res['sources']['hiv']}` · TB `{res['sources']['tb']}`.", "",
         "The deterministic decoder emits R/S; this is the QUANTITATIVE layer — coverage-valid prediction "
         "intervals. `cover_90` = held-out coverage (target 0.90); `fold_factor` = the prediction interval "
         "expressed as ×/÷ (fold-change for HIV, MIC dilutions for TB). `informative` = the model beats the "
         "marginal label (R²>0.05) — only these intervals carry genotype signal.", "",
         "| pathogen | class/drug | label | n | R² | cover_90 | interval | calibrated | informative |",
         "|---|---|---|---|---|---|---|---|---|"]
    order = {"calibrated": 0, "powered-not-calibrated": 1}
    for r in sorted(res["rows"], key=lambda r: (r["pathogen"], order.get(r["status"], 2),
                                                not r.get("informative"), -(r["n"] or 0))):
        if r["calibrated"] is None:
            L.append(f"| {r['pathogen']} | {r['unit']} | {r['label']} | {r['n']} | — | — | — | {r['status']} | — |")
            continue
        L.append(f"| {r['pathogen']} | {r['unit']} | {r['label']} | {r['n']} | {r['r2']} | "
                 f"**{r['cover_90']}** | ×/÷{r['fold_factor']} | {'YES' if r['calibrated'] else 'no'} | "
                 f"{'YES' if r.get('informative') else 'no (≈marginal)'} |")
    L += ["", "## Honest caveats"] + [f"- {c}" for c in res["honest_caveats"]]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hiv", default=None)
    ap.add_argument("--tb", default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    today = _date.today().isoformat()
    res = build(a.hiv, a.tb)
    out = a.out or (REPO / "wiki" / f"quantitative_decoder_capability_map_{today}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
