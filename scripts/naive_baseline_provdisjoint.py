"""Paired naive-AMRFinder baseline vs the FROZEN call_resistance rule across the 10
provenance-disjoint SCORED cells (the PRIMARY trust surface).

Generalizes scripts/oxford_naive_baseline.py from the single Oxford external cohort to the
report card's 10 NCBI-PD provenance-disjoint cells. The committed rail: a curated policy
layer over a curated-DB tool must BEAT naive use of that tool -- and the report card's
headline numbers are the frozen rule alone, so this closes the wrapper-vs-tool gap on the
surface that matters most. Reuses each cell's ALREADY-CACHED per-isolate AMRFinder main.tsv
(data/raw/<slug>/amrfinder_runs/<acc>/main.tsv) -- NO Docker, NO re-download.

- naive  = R iff ANY main.tsv row whose Class is in mic_tiers.amrfinder_classes_for(drug)
           (broad class match; no subclass/point/threshold refinement; no abstain).
- frozen = the shipped call_resistance(main.tsv, drug, organism=registry_organism) prediction.
- metric = balanced accuracy (sens+spec)/2 (naive games one axis: call everything R).

M4 RECONCILIATION GUARD: the frozen confusion matrix recomputed here from the cache MUST
equal the committed artifact's metrics{tp,fp,tn,fn} before any delta is trusted; a mismatch
(cache incomplete / org-drug drift) -> the cell is reported RECONCILE_MISMATCH, no delta.
No frozen-surface change.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.mic_tiers import amrfinder_classes_for
from dna_decode.eval.amr_rules import call_resistance
from scripts.independent_cohort_validate import _conf

# (committed provdisjoint artifact, cached-runs slug dir) for each of the 10 SCORED cells.
CELLS = [
    ("wiki/provenance_disjoint_validation_campylobacter_cipro_2026-06-10.json",
     "data/raw/campylobacter_provdisjoint_ciprofloxacin"),
    ("wiki/provenance_disjoint_validation_escherichia_coli_shigella_cipro_2026-06-10.json",
     "data/raw/escherichia_coli_shigella_provdisjoint_ciprofloxacin"),
    ("wiki/provenance_disjoint_validation_escherichia_coli_shigella_ceftr_2026-06-12.json",
     "data/raw/escherichia_coli_shigella_provdisjoint_ceftriaxone"),
    ("wiki/provenance_disjoint_validation_escherichia_coli_shigella_genta_2026-06-12.json",
     "data/raw/escherichia_coli_shigella_provdisjoint_gentamicin"),
    ("wiki/provenance_disjoint_validation_escherichia_coli_shigella_tetra_2026-06-12.json",
     "data/raw/escherichia_coli_shigella_provdisjoint_tetracycline"),
    ("wiki/provenance_disjoint_validation_klebsiella_cipro_2026-06-10.json",
     "data/raw/klebsiella_provdisjoint_ciprofloxacin"),
    ("wiki/provenance_disjoint_validation_klebsiella_ceftr_2026-06-10.json",
     "data/raw/klebsiella_provdisjoint_ceftriaxone"),
    ("wiki/provenance_disjoint_validation_klebsiella_genta_2026-06-10.json",
     "data/raw/klebsiella_provdisjoint_gentamicin"),
    ("wiki/provenance_disjoint_validation_klebsiella_merop_2026-06-10.json",
     "data/raw/klebsiella_provdisjoint_meropenem"),
    ("wiki/provenance_disjoint_validation_klebsiella_tetra_2026-06-10.json",
     "data/raw/klebsiella_provdisjoint_tetracycline"),
]


def read_selected(slug: Path) -> list[tuple[str, str]]:
    """selected.tsv rows: <accession>\\t<R|S> (no header)."""
    out = []
    for ln in (slug / "selected.tsv").read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split("\t")
        if len(parts) >= 2 and parts[1].strip() in ("R", "S"):
            out.append((parts[0].strip(), parts[1].strip()))
    return out


def naive_predict(main_tsv: Path, drug: str) -> str:
    """R iff ANY determinant row's Class is in the drug's AMR-class set (no refinement)."""
    wanted = amrfinder_classes_for(drug)
    with open(main_tsv, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ci = header.index("Class")
        except ValueError:
            return "S"
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if ci < len(parts) and parts[ci].strip().upper() in wanted:
                return "R"
    return "S"


def _cm(conf: dict) -> tuple:
    return (conf["tp"], conf["fp"], conf["tn"], conf["fn"])


def main() -> int:
    cells_out = {}
    for artifact_path, slug_str in CELLS:
        art = json.loads(Path(artifact_path).read_text(encoding="utf-8"))
        drug = art["drug"]
        reg_org = art.get("registry_organism") or art["organism"]
        slug = Path(slug_str)
        runs = slug / "amrfinder_runs"
        key = f"{art['organism']}:{drug}"

        sel = read_selected(slug)
        frozen_pairs, naive_pairs, missing = [], [], 0
        for acc, label in sel:
            mt = runs / acc / "main.tsv"
            if not mt.exists():
                missing += 1
                continue
            y = 1 if label == "R" else 0
            fz = call_resistance(mt, drug, organism=reg_org)["prediction"]
            nv = naive_predict(mt, drug)
            frozen_pairs.append((fz, y))
            naive_pairs.append((nv, y))

        fz_conf, nv_conf = _conf(frozen_pairs), _conf(naive_pairs)
        committed = art.get("metrics", {})
        committed_cm = (committed.get("tp"), committed.get("fp"),
                        committed.get("tn"), committed.get("fn"))
        reconciled = _cm(fz_conf) == committed_cm

        entry = {"drug": drug, "registry_organism": reg_org, "n_selected": len(sel),
                 "n_cached_missing": missing, "frozen": fz_conf, "naive": nv_conf,
                 "committed_cm": list(committed_cm), "recomputed_cm": list(_cm(fz_conf)),
                 "reconciled": reconciled}
        if reconciled and fz_conf["sens"] is not None and nv_conf["sens"] is not None:
            fz_balacc = round((fz_conf["sens"] + fz_conf["spec"]) / 2, 4)
            nv_balacc = round((nv_conf["sens"] + nv_conf["spec"]) / 2, 4)
            d = round(fz_balacc - nv_balacc, 4)
            entry.update(frozen_balacc=fz_balacc, naive_balacc=nv_balacc, delta_balacc=d,
                         value_add_verdict=("CURATED_LAYER_ADDS_VALUE" if d >= 0.03
                                            else "NAIVE_BEATS_CURATED" if d <= -0.03
                                            else "NAIVE_TIES_CURATED"))
            tag = entry["value_add_verdict"]
            print(f"{key:38s} n={fz_conf['n_scored']:3d} frozen_balacc={fz_balacc:.3f} "
                  f"naive_balacc={nv_balacc:.3f} delta={d:+.3f}  {tag}")
        else:
            entry["value_add_verdict"] = "RECONCILE_MISMATCH"
            print(f"{key:38s} RECONCILE_MISMATCH committed={committed_cm} recomputed={_cm(fz_conf)} "
                  f"(missing {missing}/{len(sel)} cached)")
        cells_out[key] = entry

    scored = [c for c in cells_out.values() if c["value_add_verdict"] != "RECONCILE_MISMATCH"]
    adds = [c for c in scored if c["value_add_verdict"] == "CURATED_LAYER_ADDS_VALUE"]
    out = {
        "_schema": "provdisjoint-naive-comparator-v1",
        "date": _date.today().isoformat(),
        "baseline_definition": ("naive = R iff ANY AMRFinder determinant Class in "
                                "mic_tiers.amrfinder_classes_for(drug); frozen = call_resistance rule. "
                                "Metric = balanced accuracy. M4 reconciliation guard: frozen CM must "
                                "match the committed provdisjoint artifact before any delta is trusted."),
        "rail": "validate-wrapper-vs-underlying-tool, on the primary 10-cell provenance-disjoint surface.",
        "n_cells_reconciled": len(scored),
        "n_curated_adds_value": len(adds),
        "cells": cells_out,
        "frozen_surface_changed": False,
    }
    outp = Path(f"wiki/provdisjoint_naive_comparator_{_date.today().isoformat()}.json")
    outp.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n{len(scored)}/{len(cells_out)} cells reconciled; "
          f"{len(adds)}/{len(scored)} CURATED_LAYER_ADDS_VALUE -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
