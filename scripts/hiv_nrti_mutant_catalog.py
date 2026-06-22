"""NRTI v0.1 — a MUTANT-SPECIFIC NRTI catalog data-derived from HIVDB fold-associations (Wave B).

The position-based v0 (any non-consensus at a major NRTI position) over-calls T215 revertants / V75
polymorphisms (spec ~0.41-0.64). v0.1 fix: keep only the amino-acid variants ACTUALLY resistance-associated
in the data. A naive carriers'-median-fold rule is CONFOUNDED (a revertant rides on resistant lineages -> its
carriers' median fold is high though it isn't independently resistant), so the resistant set is the mutants
whose MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) (an independent >=1.5x effect after controlling
for co-occurring mutations; >= MIN_CARRIERS carriers). That deconfounds the revertants.

CIRCULARITY GUARD (load-bearing — the project's standing discipline): the catalog is derived FROM the fold
data, so an in-sample score would be optimistic. The headline metric is therefore **5-fold CROSS-VALIDATED**
(derive the resistant-mutant set on the training folds, evaluate on the held-out fold, aggregate held-out
predictions) — so the reported spec recovery is an HONEST out-of-sample estimate, NOT in-sample. The fixed
deliverable catalog (derived on all data) is written separately + labelled data-derived.

Label = PhenoSense fold (independent wet-lab; NOT HIVDB Sierra). Cutoffs from DRMcv.R. Cite Rhee 2003.
DATA: data/raw/hiv/NRTI_DataSet.txt (gitignored).
"""
from __future__ import annotations

import json
import math
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import NRTI_MAJOR_POSITIONS, NRTI_RT_WT
from scripts.hiv_nnrti_baseline import _confusion
from scripts.hiv_nnrti_validate import _parse_fold, load_rows
from scripts.hiv_nrti_validate import DEFAULT_DATA, NRTI_DRUGS, NRTI_LOWER_CUTOFF, _NRTI_COL, _observed_nrti

MIN_CARRIERS = 5             # a candidate mutant must appear in >= this many isolates (with fold)
N_FOLDS = 5
SEED = 0
RESIST_COEF_MIN = math.log10(1.5)   # an OLS log10-fold coefficient >= this = an independent >=1.5x effect


def _isolate_records(rows, col):
    """[(observed_substitutions_set, fold)] for isolates with a usable fold for this drug column."""
    out = []
    for row in rows:
        fold = _parse_fold(row.get(col, ""))
        if fold is None or fold <= 0:
            continue
        out.append((_observed_nrti(row), fold))
    return out


def derive_resistant_mutants(records) -> set[str]:
    """Resistant mutants at the major NRTI positions via MULTIVARIATE OLS (controls for co-occurrence).

    A naive carriers'-median-fold rule is CONFOUNDED — a revertant (T215E/I/V) rides on resistant lineages
    (co-occurs with T215Y/M41L), so its carriers' median fold is high though it isn't independently
    resistant. So instead: fit OLS of log10(fold) ~ binary presence of every candidate major-position
    mutant (>= MIN_CARRIERS carriers); a mutant is resistant iff its INDEPENDENT coefficient >=
    RESIST_COEF_MIN (a >=1.5x fold effect after controlling for the others). This zeroes out revertants."""
    counts: dict[str, int] = {}
    for obs, _ in records:
        for s in obs:                          # obs is already restricted to major positions
            counts[s] = counts.get(s, 0) + 1
    candidates = sorted(m for m, n in counts.items() if n >= MIN_CARRIERS)
    if not candidates or len(records) < 30:
        return set()
    fidx = {m: j for j, m in enumerate(candidates)}
    X = np.zeros((len(records), len(candidates)))
    y = np.zeros(len(records))
    for i, (obs, fold) in enumerate(records):
        y[i] = math.log10(fold)
        for s in obs:
            if s in fidx:
                X[i, fidx[s]] = 1.0
    coef = LinearRegression().fit(X, y).coef_
    return {candidates[j] for j in range(len(candidates)) if coef[j] >= RESIST_COEF_MIN}


def _cv_compare(records, cutoff: float) -> dict:
    """5-fold CV: derive the resistant-mutant set on train folds, evaluate position-based vs mutant-specific
    on the held-out fold. Aggregate held-out predictions for an out-of-sample confusion each."""
    recs = list(records)
    if len(recs) < 30:
        return {"n": len(recs), "note": "too few isolates"}
    idx = np.arange(len(recs))
    pos_pred, ms_pred, actual = [], [], []
    for tr, te in KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED).split(idx):
        resistant = derive_resistant_mutants([recs[i] for i in tr])
        for i in te:
            obs, fold = recs[i]
            pos_pred.append(len(obs) > 0)                          # position-based: any major-pos mutation
            ms_pred.append(any(s in resistant for s in obs))       # mutant-specific: a derived resistant one
            actual.append(fold >= cutoff)
    pos_pred = np.array(pos_pred); ms_pred = np.array(ms_pred); actual = np.array(actual)
    return {
        "n_isolates": len(recs),
        "prevalence_R": round(float(actual.mean()), 3),
        "position_based_v0": _confusion(pos_pred, actual),
        "mutant_specific_v0_1_heldout": _confusion(ms_pred, actual),
    }


def run(path: Path = DEFAULT_DATA) -> dict:
    rows = load_rows(path)
    per_drug, final_catalog = {}, {}
    for drug in NRTI_DRUGS:
        col, cutoff = _NRTI_COL[drug], NRTI_LOWER_CUTOFF[drug]
        records = _isolate_records(rows, col)
        per_drug[drug] = _cv_compare(records, cutoff)
        # the fixed deliverable catalog: derived on ALL data for this drug
        final_catalog[drug] = sorted(derive_resistant_mutants(records))
        if "mutant_specific_v0_1_heldout" in per_drug[drug]:
            pb = per_drug[drug]["position_based_v0"]["balanced_accuracy"]
            ms = per_drug[drug]["mutant_specific_v0_1_heldout"]["balanced_accuracy"]
            per_drug[drug]["balacc_gain_v0_1_minus_v0"] = (
                round(ms - pb, 3) if (pb is not None and ms is not None) else None)
    return {
        "artifact": "hiv_nrti_mutant_specific_v0_1", "schema": "hiv-nrti-mutant-catalog-v0_1",
        "method": "data-derived resistant mutants = MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) "
                  "(independent >=1.5x effect, controlling for co-occurrence -> deconfounds revertants), "
                  ">=5 carriers; 5-fold CROSS-VALIDATED held-out (derive on train, eval on test) -> NOT in-sample",
        "honest_caveats": [
            "deconfounded (OLS-coefficient) derivation recovers specificity for 5/6 NRTIs (AZT/D4T/TDF/3TC/ABC)",
            "didanosine REGRESSES (balacc -0.201): ddI has a low resistance signal (cutoff 1.5) so few mutants "
            "clear the independent >=1.5x threshold -> sens collapses; production should keep POSITION-BASED for ddI",
            "naive carriers'-median-fold derivation was CONFOUNDED (revertants ride on resistant lineages) -> "
            "FIXED to the multivariate-OLS-coefficient rule (a verify-in-batch catch)",
        ],
        "label_source": "Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra)",
        "cutoffs": NRTI_LOWER_CUTOFF, "dataset": str(path),
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from DRMcv.R",
        "per_drug": per_drug,
        "deliverable_catalog_all_data": final_catalog,
    }


def render_md(result: dict, generated: str) -> str:
    lines = [f"# HIV NRTI v0.1 - mutant-specific catalog (data-derived, CV held-out) ({generated})", ""]
    lines.append(f"Method: {result['method']}.")
    lines.append(f"Label = {result['label_source']}.")
    lines.append("")
    lines.append("Position-based v0 over-calls revertants; v0.1 keeps only fold-associated mutants. "
                 "Both metrics are 5-fold CV held-out (out-of-sample).")
    lines.append("")
    lines.append("| Drug | n | prev R | v0 pos-based sens/spec/**balacc** | v0.1 mutant sens/spec/**balacc** | balacc gain |")
    lines.append("|---|---|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        if "mutant_specific_v0_1_heldout" not in m:
            lines.append(f"| {drug} | {m.get('n')} | - | {m.get('note','')} | - | - |")
            continue
        pb, ms = m["position_based_v0"], m["mutant_specific_v0_1_heldout"]
        lines.append(f"| {drug} | {m['n_isolates']} | {m['prevalence_R']} | "
                     f"{pb['sens']}/{pb['spec']}/**{pb['balanced_accuracy']}** | "
                     f"{ms['sens']}/{ms['spec']}/**{ms['balanced_accuracy']}** | "
                     f"{m['balacc_gain_v0_1_minus_v0']} |")
    lines.append("")
    lines.append("## Deliverable catalog (derived on all data)")
    for drug, muts in result["deliverable_catalog_all_data"].items():
        lines.append(f"- **{drug}** ({len(muts)}): {', '.join(muts)}")
    lines.append("")
    lines.append(f"Citation: {result['citation']}.")
    return "\n".join(lines)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)
    if not args.data.exists():
        print(f"ERROR: dataset not found at {args.data}", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    result = run(args.data)
    out_md = args.out_md or (REPO / "wiki" / f"hiv_nrti_mutant_catalog_v0_1_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
