"""Systematic RULE-GAP screen: which determinants do the deployed cells MISS, and are they real gaps?

WHY
The 2026-08-24 prospective cohort exposed one gap by hand: E. coli x gentamicin missed 24 of 28 resistant
isolates that carry a 16S rRNA methyltransferase (`rmt*`/`armA`) which the frozen `Subclass=GENTAMICIN`
rule cannot see. That was found by eyeballing 28 false negatives on N=37.

`wiki/determinant_blindness_atlas.json` says the same shape may be widespread: across all 12 NCBI-PD
external-validation cells, EVERY invisible resistant isolate is `rule_limited` (a determinant IS present,
the rule just does not count it) and ZERO are `truly_invisible`. But the atlas is an AGGREGATE fraction --
it never names WHICH determinant. This screen names them, at up to N=110 instead of 37.

THE STATISTIC (the one that made the rmt finding convincing, not the raw count)
A token frequent among missed-R isolates proves nothing on its own -- passengers are frequent everywhere.
What separated `rmt` was that it was frequent in the missed-R set and ABSENT from the susceptible set. So
each candidate is scored on BOTH:
    miss_rate = carriers among missed-R / missed-R
    s_rate    = carriers among measured-S / measured-S
and ranked by the gap `miss_rate - s_rate`. A token common in S is a passenger, not a mechanism.

HONEST SCOPE -- this is HYPOTHESIS-GENERATING, not a verdict.
  * DESCRIPTIVE: it reports co-occurrence on one cohort. It does not establish that the determinant CAUSES
    the resistance, and a high gap score can still be a clonal artifact (a determinant that travels with a
    resistant lineage for unrelated reasons). The atlas's own `/innovate` note records that a burden-based
    RESCUE died exactly this way -- pooled signal collapsing inside SNP clusters.
  * So nothing here may enter `dna_decode/amr/uncounted.py::_MEASURED_GAPS` on this output alone. That
    table's bar is a MEASURED gap with a citable artifact; a screen hit is a candidate for that work.
  * A high miss rate is EXPECTED and correct for known determinant-invisible mechanisms (gonococcal
    azithromycin is mtr-efflux driven) -- there the honest answer stays "disclose the blindness".

PURE + OFFLINE: reads the committed cohort TSVs and calls the deployed caller functions on symbol lists.
No Docker, no network, no AMRFinder re-run. Frozen surface is READ-only.

Run: uv run python scripts/rule_gap_screen.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

MIN_MISSED = 3          # below this the rates are not interpretable; the cell is reported as underpowered
MIN_MISS_RATE = 0.50    # carried by at least half the missed-R isolates
MAX_S_RATE = 0.20       # and rare among the measured-S isolates (the passenger filter)


def load_cell(cohort_dir: Path):
    labels = {r["biosample"]: r for r in csv.DictReader(
        open(cohort_dir / "cohort.tsv", encoding="utf-8"), delimiter="\t")}
    dets = {r["biosample"]: [s for s in (r["determinants"] or "").split(";") if s]
            for r in csv.DictReader(open(cohort_dir / "determinants.tsv", encoding="utf-8"), delimiter="\t")}
    return labels, dets


def screen_drug(labels: dict, dets: dict, drug: str, call_fn) -> dict:
    """Partition the cell by (label, prediction) and rank candidate gap determinants. PURE."""
    missed_r, called_r, true_s = [], [], []
    for bs, row in labels.items():
        lab = row.get(drug, "")
        if lab not in ("R", "S"):
            continue
        d = dets.get(bs, [])
        try:
            pred = str(call_fn(d)["prediction"]).upper()
        except Exception:      # noqa: BLE001 -- a caller that refuses this isolate is not a gap signal
            continue
        if lab == "R":
            (called_r if pred == "R" else missed_r).append(d)
        else:
            true_s.append(d)

    n_missed, n_s = len(missed_r), len(true_s)
    out = {"drug": drug, "n_R": len(missed_r) + len(called_r), "n_missed_R": n_missed,
           "n_called_R": len(called_r), "n_S": n_s, "candidates": [], "status": "screened"}
    if n_missed < MIN_MISSED:
        out["status"] = "underpowered" if n_missed else "no_missed_R"
        return out

    miss_ct, s_ct = Counter(), Counter()
    for d in missed_r:
        miss_ct.update(set(d))
    for d in true_s:
        s_ct.update(set(d))

    # a token the rule already counts cannot be a gap: it is present in isolates the rule CALLED
    counted = set()
    for d in called_r:
        counted.update(d)

    cands = []
    for tok, n in miss_ct.items():
        miss_rate = n / n_missed
        s_rate = (s_ct[tok] / n_s) if n_s else 0.0
        if miss_rate >= MIN_MISS_RATE and s_rate <= MAX_S_RATE:
            cands.append({"determinant": tok, "n_in_missed_R": n, "miss_rate": round(miss_rate, 3),
                          "s_rate": round(s_rate, 3), "gap": round(miss_rate - s_rate, 3),
                          "also_in_called_R": tok in counted})
    out["candidates"] = sorted(cands, key=lambda c: -c["gap"])
    return out


def main() -> int:
    from scripts.score_ncbipd_extval import ORGANISMS

    report = {"_schema": "rule-gap-screen-v1", "date": _date.today().isoformat(),
              "thresholds": {"min_missed_R": MIN_MISSED, "min_miss_rate": MIN_MISS_RATE,
                             "max_s_rate": MAX_S_RATE},
              "honest_scope": ("DESCRIPTIVE co-occurrence on ONE cohort; hypothesis-generating only. A hit "
                               "may still be a clonal passenger. Nothing here may enter _MEASURED_GAPS "
                               "without its own validation."),
              "cells": []}

    for key, (name, cohort_dir, drugs) in sorted(ORGANISMS.items()):
        cdir = ROOT / cohort_dir
        if not (cdir / "cohort.tsv").exists():
            print(f"{name}: cohort absent ({cohort_dir}) -- skipped")
            continue
        labels, dets = load_cell(cdir)
        for drug, fn in sorted(drugs.items()):
            cell = screen_drug(labels, dets, drug, fn)
            cell["organism"] = name
            report["cells"].append(cell)
            top = cell["candidates"][:3]
            desc = ", ".join(f"{c['determinant']} (gap {c['gap']})" for c in top) or "-"
            print(f"{name:26s} {drug:14s} R={cell['n_R']:4d} missed={cell['n_missed_R']:4d} "
                  f"[{cell['status']}]  {desc}")

    out = ROOT / "wiki" / f"rule_gap_screen_{report['date']}.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    n_hits = sum(1 for c in report["cells"] if c["candidates"])
    print(f"\n{n_hits}/{len(report['cells'])} cells have >=1 candidate -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
