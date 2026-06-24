"""HIV PI v0.1 — a MUTANT-SPECIFIC protease-inhibitor catalog data-derived from HIVDB fold-associations.

Mirrors the SHIPPED NRTI v0.1 arc (`scripts/hiv_nrti_mutant_catalog.py`) for the 8 protease inhibitors.
The position-based v0 (`PI_CLASS`, any non-consensus at a major protease position 30/32/33/46/47/48/50/54/
76/82/84/88/90) over-calls benign polymorphisms + accessory mutations that ride on resistant lineages
(the `hiv_pi_baseline_vs_ols` memo bounds the headroom: OLS beats the position catalog by +0.061..+0.131
balanced accuracy across all 8 PIs). v0.1 fix: keep only the amino-acid variants ACTUALLY resistance-
associated in the data. A naive carriers'-median-fold rule is CONFOUNDED (an accessory mutation co-occurs
with a major one -> its carriers' median fold is high though it isn't independently resistant), so the
resistant set is the mutants whose MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) (an independent
>=1.5x effect after controlling for co-occurring mutations; >= MIN_CARRIERS carriers). That deconfounds the
accessory/polymorphic positions.

DELTA-HONEST CUTOFF (load-bearing — the honest difference from the NRTI builder): there is NO per-drug
clinical PI cutoff sourced in this repo (the NRTI builder had `NRTI_LOWER_CUTOFF` from DRMcv.R; PI does
not). So BOTH the v0 position-based call AND the v0.1 mutant-specific call are scored at the SAME UNIFORM
illustrative fold>=3 boundary (`ILLUSTRATIVE_FOLD_CUTOFF`, the one `hiv_targetsite_validate` /
`hiv_targetsite_baseline` already use + flag). The headline is therefore the v0->v0.1 balanced-accuracy
GAIN at that fixed cutoff (a fair refinement signal), NOT an absolute-calibration claim. Sourcing per-drug
PI clinical cutoffs would upgrade this to absolute calibration (a v0.2 item; not fabricated here).

CIRCULARITY GUARD (the project's standing discipline): the catalog is derived FROM the fold data, so an
in-sample score would be optimistic. The headline metric is therefore 5-fold CROSS-VALIDATED (derive the
resistant-mutant set on the training folds, evaluate on the held-out fold, aggregate held-out predictions)
-> an HONEST out-of-sample estimate, NOT in-sample. The fixed deliverable catalog (derived on all data) is
written separately + labelled data-derived.

Label = PhenoSense fold (independent wet-lab; NOT HIVDB Sierra). INSTI deferred (its OLS-vs-catalog deltas
are thin/mixed: EVG -0.026, BIC -0.011, only CAB +0.314 at n=64 unstable; see hiv_insti_baseline_vs_ols).
DATA: data/raw/hiv/PI_DataSet.txt (gitignored; cite Rhee 2003).
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

from scripts.hiv_nnrti_baseline import _confusion  # noqa: E402
from scripts.hiv_nnrti_validate import _parse_fold, load_rows  # noqa: E402
from scripts.hiv_targetsite_validate import (  # noqa: E402
    ILLUSTRATIVE_FOLD_CUTOFF, _CLASS_SPECS, _observed_mutations,
)

MIN_CARRIERS = 5             # a candidate mutant must appear in >= this many isolates (with fold)
N_FOLDS = 5
SEED = 0
RESIST_COEF_MIN = math.log10(1.5)   # an OLS log10-fold coefficient >= this = an independent >=1.5x effect
DEFAULT_DATA = REPO / "data" / "raw" / "hiv" / "PI_DataSet.txt"


def _isolate_records(rows, col, cls):
    """[(observed_major_position_substitutions, fold)] for isolates with a usable fold for this drug column."""
    out = []
    for row in rows:
        fold = _parse_fold(row.get(col, ""))
        if fold is None or fold <= 0:
            continue
        out.append((_observed_mutations(row, cls), fold))
    return out


def derive_resistant_mutants(records) -> set[str]:
    """Resistant mutants at the major PI positions via MULTIVARIATE OLS (controls for co-occurrence).

    A naive carriers'-median-fold rule is CONFOUNDED — an accessory mutation co-occurs with a major one,
    so its carriers' median fold is high though it isn't independently resistant. So instead: fit OLS of
    log10(fold) ~ binary presence of every candidate major-position mutant (>= MIN_CARRIERS carriers); a
    mutant is resistant iff its INDEPENDENT coefficient >= RESIST_COEF_MIN (a >=1.5x fold effect after
    controlling for the others). This zeroes out the accessory/polymorphic riders."""
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
            pos_pred.append(len(obs) > 0)                          # position-based v0: any major-pos mutation
            ms_pred.append(any(s in resistant for s in obs))       # mutant-specific v0.1: a derived resistant one
            actual.append(fold >= cutoff)
    pos_pred = np.array(pos_pred); ms_pred = np.array(ms_pred); actual = np.array(actual)
    return {
        "n_isolates": len(recs),
        "prevalence_R": round(float(actual.mean()), 3),
        "position_based_v0": _confusion(pos_pred, actual),
        "mutant_specific_v0_1_heldout": _confusion(ms_pred, actual),
    }


def run(path: Path = DEFAULT_DATA) -> dict:
    cls, _, drug_cols = _CLASS_SPECS["PI"]
    cutoff = ILLUSTRATIVE_FOLD_CUTOFF
    rows = load_rows(path)
    per_drug, final_catalog = {}, {}
    for drug, col in drug_cols.items():
        records = _isolate_records(rows, col, cls)
        per_drug[drug] = _cv_compare(records, cutoff)
        # the fixed deliverable catalog: derived on ALL data for this drug
        final_catalog[drug] = sorted(derive_resistant_mutants(records))
        if "mutant_specific_v0_1_heldout" in per_drug[drug]:
            pb = per_drug[drug]["position_based_v0"]["balanced_accuracy"]
            ms = per_drug[drug]["mutant_specific_v0_1_heldout"]["balanced_accuracy"]
            per_drug[drug]["balacc_gain_v0_1_minus_v0"] = (
                round(ms - pb, 3) if (pb is not None and ms is not None) else None)
    gains = [m["balacc_gain_v0_1_minus_v0"] for m in per_drug.values()
             if m.get("balacc_gain_v0_1_minus_v0") is not None]
    n_improved = sum(1 for g in gains if g >= 0)
    return {
        "artifact": "hiv_pi_mutant_specific_v0_1", "schema": "hiv-pi-mutant-catalog-v0_1",
        "method": "data-derived resistant mutants = MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) "
                  "(independent >=1.5x effect, controlling for co-occurrence -> deconfounds accessory/"
                  "polymorphic riders), >=5 carriers; 5-fold CROSS-VALIDATED held-out (derive on train, "
                  "eval on test) -> NOT in-sample",
        "cutoff_note": "DELTA-HONEST: both v0 (position-based) and v0.1 (mutant-specific) scored at the SAME "
                       "UNIFORM illustrative fold>=3 (no per-drug PI clinical cutoff sourced in-repo); the "
                       "headline is the v0->v0.1 balanced-accuracy GAIN at that fixed cutoff, NOT an "
                       "absolute-calibration claim. Per-drug PI clinical cutoffs would upgrade to absolute "
                       "calibration (a v0.2 item; not fabricated here).",
        "n_drugs": len(drug_cols), "n_drugs_improved_or_held": n_improved,
        "mean_balacc_gain": round(float(np.mean(gains)), 3) if gains else None,
        "honest_caveats": [
            "DELTA-HONEST at the uniform fold>=3 cutoff (no per-drug PI clinical breakpoint sourced); the "
            "GAIN is the signal, not the absolute balacc",
            "deconfounded (multivariate-OLS-coefficient) derivation, mirroring the shipped NRTI v0.1 arc; "
            "5-fold CV held-out -> out-of-sample, not optimistic in-sample",
            "in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)",
            "INSTI deferred: its OLS-vs-catalog deltas are thin/mixed (EVG -0.026, BIC -0.011, only CAB "
            "+0.314 at n=64 unstable); PI is where the headroom is uniform across all 8 drugs",
        ],
        "label_source": "Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra)",
        "illustrative_lower_cutoff_fold": cutoff, "dataset": str(path),
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R",
        "per_drug": per_drug,
        "deliverable_catalog_all_data": final_catalog,
    }


def render_md(result: dict, generated: str) -> str:
    lines = [f"# HIV PI v0.1 - mutant-specific catalog (data-derived, CV held-out) ({generated})", ""]
    lines.append(f"Method: {result['method']}.")
    lines.append("")
    lines.append(f"**Cutoff (delta-honest):** {result['cutoff_note']}")
    lines.append("")
    lines.append(f"Label = {result['label_source']}.")
    lines.append("")
    lines.append(f"Position-based v0 over-calls accessory/polymorphic riders at the 13 major protease "
                 f"positions; v0.1 keeps only fold-associated mutants. Both metrics are 5-fold CV held-out "
                 f"(out-of-sample). **{result['n_drugs_improved_or_held']}/{result['n_drugs']} PI drugs "
                 f"improve-or-hold; mean balacc gain {result['mean_balacc_gain']}.**")
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
        lines.append(f"- **{drug}** ({len(muts)}): {', '.join(muts) if muts else '(none cleared the threshold)'}")
    lines.append("")
    lines.append("## Honest caveats")
    lines += [f"- {c}" for c in result["honest_caveats"]]
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
    out_md = args.out_md or (REPO / "wiki" / f"hiv_pi_v0.1_validation_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
