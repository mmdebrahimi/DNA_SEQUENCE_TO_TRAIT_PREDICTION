"""Does a DETERMINISTIC conservation score compete with learned models in the molecular-phenotype regime?

The BLOSUM benchmark (scripts/dms_variant_effect_benchmark.py) showed a deterministic SUBSTITUTION-SEVERITY
floor of ~0.21 median Spearman, ~2x below AlphaMissense (~0.42). But substitution matrices are position-BLIND.
The honest deterministic competitor is a POSITION-SPECIFIC conservation score (an independent-sites model over a
per-protein MSA) -- the "Site-Independent" baseline.

Rather than hand-roll an MSA pipeline (risk: a non-canonical baseline / parameter-search against the learned
models), this reads ProteinGym's AUTHORITATIVE, published per-assay zero-shot Spearman table (v1.x), which
already ran the Site-Independent deterministic baseline AND the learned models (ESM-1v / EVE / GEMME) on
IDENTICAL variant rows. This module extracts + audits that comparison by selection type (Activity=function /
Binding / Stability=abundance / ...) on the human assays, and anchors it to this project's own BLOSUM floor +
AlphaMissense numbers. Grounded in the canonical benchmark; non-circular (ProteinGym is not tuned to this
question). Inputs (free, committed by ProteinGym) cached on D:.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

DETERMINISTIC = "Site-Independent"
LEARNED = ["ESM-1v (ensemble)", "EVE (ensemble)", "GEMME"]
SELECTION_ORDER = ["Activity", "Binding", "Stability", "Expression", "OrganismalFitness"]


def load(reference_csv: Path, spearman_csv: Path) -> pd.DataFrame:
    ref = pd.read_csv(reference_csv)
    sp = pd.read_csv(spearman_csv).rename(columns={"DMS ID": "DMS_id"})
    selc = "coarse_selection_type" if "coarse_selection_type" in ref.columns else "selection_type"
    keep = ["DMS_id", "UniProt_ID", selc, "taxon"]
    m = sp.merge(ref[keep], on="DMS_id", how="left").rename(columns={selc: "selection"})
    return m


def _median(g: pd.DataFrame, col: str) -> float | None:
    v = pd.to_numeric(g[col], errors="coerce").dropna()
    return round(float(v.median()), 4) if len(v) else None


def run(reference_csv: Path, spearman_csv: Path, blosum_json: Path | None = None,
        am_json: Path | None = None, human_only: bool = True) -> dict:
    m = load(reference_csv, spearman_csv)
    if human_only:
        m = m[m["taxon"].astype(str).str.contains("Human", case=False, na=False)]
    rows = []
    for st in SELECTION_ORDER:
        g = m[m["selection"] == st]
        if len(g) < 3:
            continue
        rec = {"selection": st, "n": len(g), "site_independent": _median(g, DETERMINISTIC)}
        for lm in LEARNED:
            rec[lm] = _median(g, lm)
        # deterministic competes on this modality if it >= the WORST learned model here (beats a learned model)
        learned_vals = [rec[lm] for lm in LEARNED if rec[lm] is not None]
        rec["beats_a_learned_model"] = bool(rec["site_independent"] is not None and learned_vals
                                            and rec["site_independent"] >= min(learned_vals))
        rec["gap_to_best_learned"] = (round(max(learned_vals) - rec["site_independent"], 4)
                                      if rec["site_independent"] is not None and learned_vals else None)
        rows.append(rec)
    overall = {"n": len(m), "site_independent": _median(m, DETERMINISTIC),
               **{lm: _median(m, lm) for lm in LEARNED}}
    # project anchors (this repo's own numbers, for context)
    anchors = {}
    for name, path in (("blosum_floor", blosum_json), ("alphamissense", am_json)):
        if path and Path(path).exists():
            d = json.loads(Path(path).read_text(encoding="utf-8"))
            anchors[name] = d.get("overall_median_spearman") or d.get("am_median_spearman_paired") or \
                d.get("am_median_spearman")
    act = next((r for r in rows if r["selection"] == "Activity"), None)
    verdict = "NO_FUNCTION_ROW"
    if act:
        si = act["site_independent"]
        any_beat = any(r["beats_a_learned_model"] for r in rows)
        if si is not None and si >= 0.45:
            verdict = "DETERMINISTIC_CONSERVATION_MATCHES_TARGET"
        elif si is not None and si >= 0.35 and any_beat:
            verdict = "DETERMINISTIC_CONSERVATION_LARGELY_COMPETES"
        else:
            verdict = "DETERMINISTIC_CONSERVATION_TRAILS"
    return {"predictor_deterministic": DETERMINISTIC, "learned_reference": LEARNED,
            "n_human_assays": len(m), "by_selection": rows, "overall_human": overall,
            "project_anchors": anchors,
            "function_site_independent": act["site_independent"] if act else None,
            "function_best_learned": (max(v for lm in LEARNED if (v := act.get(lm)) is not None)
                                      if act else None),
            "verdict": verdict}


def main(argv=None) -> int:
    import argparse
    DG = Path("D:/dna_decode_cache/proteingym")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, default=DG / "pg_reference.csv")
    ap.add_argument("--spearman", type=Path, default=DG / "pg_spearman_dms.csv")
    ap.add_argument("--blosum", type=Path, default=REPO / "wiki" / "dms_benchmark_big_scores.json")
    ap.add_argument("--am", type=Path, default=REPO / "wiki" / "dms_am_big_scores.json")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "dms_conservation_benchmark_scores.json")
    a = ap.parse_args(argv)
    res = run(a.reference, a.spearman, a.blosum, a.am)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"human assays={res['n_human_assays']} | verdict={res['verdict']}")
    print(f"{'selection':16} {'n':>3} {'Site-Indep':>10} " + " ".join(f"{lm.split()[0]:>8}" for lm in LEARNED)
          + " beats?")
    for r in res["by_selection"]:
        print(f"{r['selection']:16} {r['n']:>3} {r['site_independent']:>10} "
              + " ".join(f"{r[lm]:>8}" for lm in LEARNED) + f"  {r['beats_a_learned_model']}")
    print(f"\nfunction: Site-Independent={res['function_site_independent']} vs best-learned="
          f"{res['function_best_learned']} | anchors={res['project_anchors']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
