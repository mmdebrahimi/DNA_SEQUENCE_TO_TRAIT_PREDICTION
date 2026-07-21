"""Determinant-blindness atlas — quantify + SHIP the decoder's invisible fraction as an honesty surface.

The `/innovate` run (2026-07-21) killed the tempting move (rescue the false-negative ceiling by scoring the
determinants the DRUG_RULE filters): it dies for expression-driven resistance (azithro mtr, burden gap -0.6)
and is a pure CLONAL confound for cumulative-chromosomal resistance (tet pooled burden +6.2 collapses to -1.0
within SNP clusters). The ONE surviving move is the reframe: the blindness is not a bug to hide, it is an
OUTPUT to disclose. This atlas computes, per (organism, drug), the **determinant-invisible fraction** — of
the measured-R isolates, the fraction the deployed cell calls non-R because no rule-firing catalog determinant
is present — and splits it into TRULY-INVISIBLE (the genome carries zero determinant token at all) vs
RULE-LIMITED (a determinant is present but not one the rule counts). It is DESCRIPTIVE, not predictive, so it
is immune to the clonal confound that killed the burden rule.

Reuses the reusable NCBI-PD external-validation substrate (`score_ncbipd_extval.ORGANISMS`) verbatim: same
cohort/determinants TSVs, same per-isolate call path. No compute wall (pure metadata read).

  uv run python scripts/build_blindness_atlas.py
  -> wiki/determinant_blindness_atlas.{md,json}  (sorted most-blind first; NO aggregate headline)
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.score_ncbipd_extval import ORGANISMS  # noqa: E402


def blindness_row(triples: list[tuple[str, str, bool]]) -> dict:
    """PURE core. triples = [(label 'R'/'S', prediction, has_any_determinant), ...] for one (organism,drug).

    Determinant-invisible = a measured-R isolate the deployed cell calls non-R (no rule-firing determinant).
    Split: truly-invisible (genome carries ZERO determinant token) vs rule-limited (a determinant is present
    but not one the rule counts). invisible_fraction is over measured-R; None when n_R == 0 (unscorable)."""
    r_isolates = [t for t in triples if t[0] == "R"]
    n_r = len(r_isolates)
    invisible = [t for t in r_isolates if str(t[1]).upper() != "R"]
    n_inv = len(invisible)
    n_truly = sum(1 for t in invisible if not t[2])
    n_rule_limited = n_inv - n_truly
    return {
        "n_R": n_r,
        "n_invisible": n_inv,
        "invisible_fraction": (round(n_inv / n_r, 3) if n_r else None),
        "n_truly_invisible": n_truly,        # measured-R with zero determinant token at all
        "n_rule_limited": n_rule_limited,    # measured-R with a determinant present but not rule-counted
    }


def _cell_triples(cohort_dir: str, call_fn, drug: str) -> list[tuple[str, str, bool]]:
    labels = {r["biosample"]: r for r in csv.DictReader(open(f"{cohort_dir}/cohort.tsv"), delimiter="\t")}
    dets = {r["biosample"]: [s for s in (r["determinants"] or "").split(";") if s]
            for r in csv.DictReader(open(f"{cohort_dir}/determinants.tsv"), delimiter="\t")}
    triples = []
    for bs, row in labels.items():
        rs = row.get(drug, "")
        if rs not in ("R", "S"):
            continue
        d = dets.get(bs, [])
        pred = call_fn(d)["prediction"]
        if str(pred).upper() in ("R", "S"):
            triples.append((rs, pred, bool(d)))
    return triples


def build_atlas() -> list[dict]:
    rows = []
    for org, (name, cohort_dir, drugs) in ORGANISMS.items():
        if not Path(f"{cohort_dir}/cohort.tsv").exists():
            continue
        for drug, fn in drugs.items():
            triples = _cell_triples(cohort_dir, fn, drug)
            b = blindness_row(triples)
            if b["n_R"] == 0:
                continue  # unscorable for blindness (no measured-R)
            rows.append({"organism": name, "drug": drug, **b})
    # most-blind first; a fully-invisible cell (fraction 1.0) is the loudest honesty flag
    rows.sort(key=lambda x: (-(x["invisible_fraction"] or 0), x["organism"], x["drug"]))
    return rows


def main() -> int:
    rows = build_atlas()
    lines = [
        "# Determinant-blindness atlas",
        "",
        f"_Generated {_date.today().isoformat()} from the NCBI-PD external-validation cohorts "
        f"(`scripts/score_ncbipd_extval.ORGANISMS`). Per (organism, drug): of the measured-R isolates, the "
        f"fraction the deployed cell calls **non-R** because no rule-firing catalog determinant is present._",
        "",
        "**Why this exists.** The `/innovate` 2026-07-21 run KILLED the tempting move to *rescue* this "
        "false-negative ceiling by scoring the filtered determinants — it dies for expression-driven "
        "resistance (azithromycin mtr-efflux, burden gap −0.6) and is a pure **clonal confound** for "
        "cumulative-chromosomal resistance (tetracycline pooled burden +6.2 collapses to −1.0 within SNP "
        "clusters). The surviving move is to **disclose** the blindness, not hide it. This table is "
        "DESCRIPTIVE (not a predictor) → immune to that clonal confound.",
        "",
        "**Read it as:** a high `invisible fraction` means the determinant catalog structurally cannot see "
        "much of this cell's resistance — a known mechanism gap (efflux / regulatory / porin-loss / "
        "multi-locus-cumulative), not a rule bug. `truly-invisible` = the genome carries **zero** determinant "
        "token; `rule-limited` = a determinant is present but not one the rule counts. **No aggregate "
        "headline** — each cell stands alone.",
        "",
        "| organism | drug | n R | **invisible fraction** | invisible (truly / rule-limited) |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['organism']} | {r['drug']} | {r['n_R']} | **{r['invisible_fraction']}** | "
            f"{r['n_invisible']} ({r['n_truly_invisible']} / {r['n_rule_limited']}) |")
    lines += [
        "",
        "## Notes",
        "- **Descriptive honesty surface, not a predictor** — it reports where the catalog is blind; it does "
        "NOT attempt to call those isolates (the `/innovate` burden-rescue that would have is a closed "
        "negative: clonal / no-signal).",
        "- A high invisible fraction is expected + honest for known determinant-invisible mechanisms: gono "
        "azithromycin (mtr-efflux, 23S-independent — 100% invisible on this cohort), gono tetracycline "
        "(chromosomal cumulative), Klebsiella meropenem (porin-loss).",
        "- Reuses the NCBI-PD substrate (provenance-disjoint, NOT methodology-independent — same AMRFinderPlus "
        "+ same cell). NON-FROZEN cells; the frozen decoder surface is byte-unchanged.",
    ]
    Path("wiki/determinant_blindness_atlas.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    Path("wiki/determinant_blindness_atlas.json").write_text(
        json.dumps({"_schema": "determinant-blindness-atlas-v1", "date": _date.today().isoformat(),
                    "n_cells": len(rows), "cells": rows}, indent=2), encoding="utf-8")
    for r in rows:
        print(f"  {r['organism'][:22]:22s} {r['drug']:14s} nR={r['n_R']:3d} "
              f"invisible={r['invisible_fraction']} ({r['n_truly_invisible']}/{r['n_rule_limited']})")
    print(f"\n{len(rows)} cells -> wiki/determinant_blindness_atlas.{{md,json}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
