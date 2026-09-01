"""Validate the DEPLOYED gentamicin v2 rule on the leakage-gated set. Offline; exit 0 always.

WHY THIS EXISTS AS A SCRIPT. The +0.369 figure that motivated v2 was computed inline and recorded in
prose. A lock cannot rest on prose: `wiki/prospective_lock_manifest_*.json` pins the decoder by sha256,
so the number that justifies the pinned decoder has to be re-derivable against that same decoder. This
scores the REAL `call_resistance` -- not a re-implementation of it -- so what is validated is what ships.

WHAT IT MEASURES. Every accession in the committed leakage-gated census sidecar that carries a
gentamicin R/S call and has cached AMRFinder output, scored with the deployed rule and with a v1
REFERENCE re-implementation for the before/after delta. Per cell, never pooled across organisms.

THE SPECIFICITY CAVEAT IS STRUCTURAL AND IS PRINTED EVERY RUN. No S-labelled `rmt` carrier exists in
any dataset checked (these isolates, the 150 local labelled ones, or 63 publicly-labelled carriers), so
"specificity unchanged" is an ABSENCE, not a bound. A rule cannot be shown to over-call on a population
containing none of the thing it newly counts.

Run: uv run python scripts/gentamicin_v2_validate.py
"""
from __future__ import annotations

import json
import re
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
WIKI = ROOT / "wiki"

RMT_GAP = re.compile(r"^(rmt[A-H]\d*|npmA\d*)$", re.I)
REGISTRY_ORGANISM = {"Escherichia_coli_Shigella": "Escherichia_coli_Shigella",
                     "Klebsiella": "Klebsiella_pneumoniae"}


def _conf(pairs) -> dict:
    tp = sum(1 for lab, pred in pairs if lab == "R" and pred)
    fn = sum(1 for lab, pred in pairs if lab == "R" and not pred)
    tn = sum(1 for lab, pred in pairs if lab == "S" and not pred)
    fp = sum(1 for lab, pred in pairs if lab == "S" and pred)
    n = tp + fn + tn + fp
    return {"n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "acc": round((tp + tn) / n, 3) if n else None,
            "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def main() -> int:
    from gentamicin_rmt_candidate import amrfinder_index, frozen_call, read_rows

    from dna_decode.eval.amr_rules import call_resistance, rule_for

    cfg = rule_for("gentamicin")
    if not cfg.get("symbol_rescue"):
        print("the deployed gentamicin rule carries NO symbol_rescue -- this validates v2; v1 is live.")
        return 0

    census = WIKI / "unscored_genome_label_census.json"
    if not census.exists():
        print("leakage-gated census sidecar absent -- run scripts/unscored_genome_label_census.py")
        return 0
    labels = json.loads(census.read_text(encoding="utf-8")).get("labels") or {}
    idx = amrfinder_index()

    by_cell: dict[tuple, dict] = {}
    for acc, rec in sorted(labels.items()):
        lab = (rec.get("calls") or {}).get("gentamicin", "").upper()
        if lab not in ("R", "S"):
            continue
        main_tsv = idx.get(acc)
        if main_tsv is None or not Path(main_tsv).exists():
            continue
        group = rec.get("group", "")
        cell = by_cell.setdefault((group, "gentamicin"),
                                  {"v2": [], "v1": [], "rmt_R": 0, "rmt_S": 0, "rescued": []})
        rows = read_rows(Path(main_tsv))
        carries_rmt = any(RMT_GAP.match((r.get("Element symbol") or "").strip()) for r in rows)
        # DEPLOYED rule (the real surface) vs a v1 REFERENCE re-implementation for the delta.
        pred_v2 = call_resistance(Path(main_tsv), "gentamicin",
                                  organism=REGISTRY_ORGANISM.get(group))["prediction"] == "R"
        pred_v1 = frozen_call(rows)
        cell["v2"].append((lab, pred_v2))
        cell["v1"].append((lab, pred_v1))
        if carries_rmt:
            cell["rmt_R" if lab == "R" else "rmt_S"] += 1
        if pred_v2 != pred_v1:
            cell["rescued"].append({"accession": acc, "label": lab,
                                    "v1": "R" if pred_v1 else "S", "v2": "R" if pred_v2 else "S"})

    cells = []
    for (org, drug), c in sorted(by_cell.items()):
        v1, v2 = _conf(c["v1"]), _conf(c["v2"])
        cells.append({"organism": org, "drug": drug, "v1": v1, "v2": v2,
                      "delta_sens": (None if v1["sens"] is None or v2["sens"] is None
                                     else round(v2["sens"] - v1["sens"], 3)),
                      "delta_spec": (None if v1["spec"] is None or v2["spec"] is None
                                     else round(v2["spec"] - v1["spec"], 3)),
                      "rmt_carriers_R": c["rmt_R"], "rmt_carriers_S": c["rmt_S"],
                      "n_calls_changed": len(c["rescued"]), "changed": c["rescued"][:40]})

    out = {"schema": "gentamicin-v2-validation-v1", "generated": date.today().isoformat(),
           "deployed_rule": {"symbol_rescue": cfg.get("symbol_rescue"),
                             "subclass_any": sorted(cfg.get("subclass_any") or []),
                             "threshold": cfg.get("threshold")},
           "cohort": "wiki/unscored_genome_label_census.json (leakage-gated via cohort_manifest)",
           "scored_against": "the DEPLOYED dna_decode.eval.amr_rules.call_resistance, not a re-implementation",
           "cells": cells,
           "specificity_caveat": (
               "Zero S-labelled rmt carriers exist in this cohort or in any other dataset checked, so an "
               "unchanged specificity is an ABSENCE of counter-examples, not a bound on over-calling. "
               "The over-call risk of the rescue is UNTESTED, not measured to be zero.")}

    print(f"\ngentamicin v2 (rescue {cfg['symbol_rescue']}) on the leakage-gated set\n")
    for c in cells:
        v1, v2 = c["v1"], c["v2"]
        print(f"  {c['organism']} x {c['drug']}  n={v2['n']} ({v2['tp']+v2['fn']}R/{v2['tn']+v2['fp']}S)")
        print(f"     v1  acc {v1['acc']}  sens {v1['sens']}  spec {v1['spec']}")
        print(f"     v2  acc {v2['acc']}  sens {v2['sens']}  spec {v2['spec']}   "
              f"[dsens {c['delta_sens']:+} dspec {c['delta_spec']:+}]" if c["delta_sens"] is not None
              else f"     v2  acc {v2['acc']}")
        print(f"     rmt carriers: {c['rmt_carriers_R']}R / {c['rmt_carriers_S']}S "
              f"| calls changed: {c['n_calls_changed']}")
    print(f"\n  SPECIFICITY CAVEAT: {out['specificity_caveat']}")

    dest = WIKI / f"gentamicin_v2_validation_{out['generated']}.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote wiki/{dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
