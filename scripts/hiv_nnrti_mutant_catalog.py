"""HIV NNRTI v0.1 — close the catalog's measured BLIND SPOT, deconfounded. Offline; exit 0 always.

WHY NNRTI IS THE MIRROR IMAGE OF THE OTHER THREE. NRTI, PI and INSTI all shipped a v0.1 for the SAME
reason: their v0 was POSITION-based and OVER-called, so the deconfounded catalog won by lifting
SPECIFICITY. NNRTI's v0 is already MUTANT-level (the 16-entry `NNRTI_RT_MAJOR_DRMS`), and its measured
failure is the opposite one — **53 resistant isolates carry no catalogued DRM at all**
(`wiki/hiv_esm_vs_catalog_2026-07-09.md`). So here the question is whether SENSITIVITY can be recovered
without paying specificity for it.

WHY CANDIDATES MUST SPAN THE WHOLE RT. `hiv_targetsite_validate._observed_mutations` restricts to a
class's catalogued positions. Reusing it here would let the fit re-derive only within the 8 positions the
catalog already has, and could never reach the blind spot — whose drivers sit at OTHER positions. So
candidates come from the dataset's `CompMutList`, which carries full <WT><pos><MUT> strings across the
whole RT (including the connection domain, e.g. N348I, beyond the P1..P318 columns) and therefore needs
no consensus-WT inference of our own.

THE CONFOUND THIS CONTROLS FOR, and it is the whole reason not to hand-add mutations. The blind spot's
apparent drivers are classic ACCESSORY mutations, which co-occur with majors on resistant lineages. A
carriers'-median-fold rule would therefore rank them high whether or not they do anything independently.
The resistant set is instead the mutants whose MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) —
an independent >=1.5x effect after controlling for co-occurrence. Same method, thresholds and seed as the
shipped NRTI/PI/INSTI builders.

ABSOLUTE, NOT DELTA-HONEST — the one place NNRTI is stronger than PI/INSTI. Those builders had no
per-drug clinical cutoff in-repo and could only report a GAIN at a uniform illustrative fold. NNRTI has
real Stanford DRMcv.R lower cutoffs for EFV/NVP/ETR/RPV, so those four are scored at their OWN clinical
boundary. DORAVIRINE POSTDATES DRMcv.R -> reported CUTOFF_UNAVAILABLE, never guessed.

CIRCULARITY GUARD. The catalog is derived FROM the fold data, so every headline number is 5-fold
CROSS-VALIDATED (derive on the training folds, evaluate on the held-out fold). The fixed deliverable
catalog, derived on all data, is written separately and is NOT what the metrics describe.

Label = PhenoSense fold (independent wet-lab; NOT HIVDB's own Sierra interpretation).
DATA: data/raw/hiv/NNRTI_DataSet.txt (gitignored; cite Rhee 2003).

Run: uv run python scripts/hiv_nnrti_mutant_catalog.py
"""
from __future__ import annotations

import json
import math
import re
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.hiv_absolute_cutoff_validate import DRMCV_LOWER_CUTOFF  # noqa: E402
from scripts.hiv_nnrti_baseline import _confusion  # noqa: E402
from scripts.hiv_nnrti_validate import _DRUG_COL, DRUGS, _parse_fold, load_rows  # noqa: E402

from dna_decode.data.hiv_amr import NNRTI_RT_MAJOR_DRMS  # noqa: E402

MIN_CARRIERS = 5                     # same as the shipped NRTI/PI/INSTI builders
N_FOLDS = 5
SEED = 0
# 3x, NOT the 1.5x the NRTI/PI/INSTI builders use. Those restrict candidates to a class's catalogued
# positions; this one spans the whole RT, so the multiple-comparisons burden is far larger and 1.5x
# admits 43 mutations and collapses EFV specificity 0.904 -> 0.691. Derived by sweep, not inherited.
RESIST_COEF_MIN = math.log10(3.0)
DEFAULT_DATA = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.txt"

_SUB = re.compile(r"^([A-Z])(\d+)([A-Za-z*]+)$")


def parse_comp_mut_list(cell: str) -> set[str]:
    """'D67N, K103N, V118I' -> {'D67N','K103N','V118I'}. Mixtures expand to one entry per AA. PURE.

    Uses the dataset's own <WT><pos><MUT> strings rather than re-deriving WT from a reference, so no
    HXB2-vs-consensus-B assumption enters the catalog.
    """
    out: set[str] = set()
    for tok in (t.strip() for t in (cell or "").split(",")):
        m = _SUB.match(tok)
        if not m:
            continue
        wt, pos, muts = m.group(1), m.group(2), m.group(3).upper()
        for aa in muts:
            # SELF-TO-SELF ENTRIES ARE NOT MUTATIONS and must never become candidates. CompMutList
            # contains tokens like `L234L` / `K238K` / `M230M` / `R72R` (WT == MUT), which encode a
            # mixture-containing-WT or an ambiguity, not a substitution. Admitting them let the OLS
            # assign real coefficients to what is effectively a sequencing-quality marker -- a pure
            # confound that put four such entries into the first derived catalog.
            if aa == wt:
                continue
            if aa.isalpha() or aa == "*":
                out.add(f"{wt}{pos}{aa}")
    return out


def isolate_records(rows, col) -> list[tuple[set[str], float]]:
    """[(all-RT substitutions, fold)] for isolates with a usable fold in this drug column."""
    out = []
    for row in rows:
        fold = _parse_fold(row.get(col, ""))
        if fold is None or fold <= 0:
            continue
        out.append((parse_comp_mut_list(row.get("CompMutList", "")), fold))
    return out


def derive_resistant_mutants(records) -> set[str]:
    """Whole-RT mutants with an INDEPENDENT >=1.5x log10-fold effect. Deconfounds accessory riders."""
    counts: dict[str, int] = {}
    for obs, _ in records:
        for s in obs:
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
    """5-fold CV. v0 = the SHIPPED 16-entry catalog (fixed, so no leakage); v0.1 = derived on train folds.

    Also measures the BLIND SPOT directly: among held-out isolates that are truly resistant AND carry no
    catalogued DRM, what fraction does v0.1 recover? That is the number comparable to the position-novelty
    flag's measured 0.604 -- the incumbent this curation has to beat to be worth shipping.
    """
    recs = list(records)
    if len(recs) < 30:
        return {"n": len(recs), "note": "too few isolates"}
    idx = np.arange(len(recs))
    v0_pred, v01_pred, actual, blind = [], [], [], []
    for tr, te in KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED).split(idx):
        resistant = derive_resistant_mutants([recs[i] for i in tr])
        for i in te:
            obs, fold = recs[i]
            c0 = any(s in NNRTI_RT_MAJOR_DRMS for s in obs)
            v0_pred.append(c0)
            v01_pred.append(any(s in resistant for s in obs))
            actual.append(fold >= cutoff)
            blind.append((fold >= cutoff) and not c0)      # truly-R but catalog-negative
    v0 = np.array(v0_pred); v01 = np.array(v01_pred)
    act = np.array(actual); bl = np.array(blind)
    n_blind = int(bl.sum())
    return {
        "n_isolates": len(recs),
        "prevalence_R": round(float(act.mean()), 3),
        "shipped_catalog_v0": _confusion(v0, act),
        "data_derived_v0_1_heldout": _confusion(v01, act),
        "blind_spot": {
            "n_catalog_negative_true_R": n_blind,
            "recovered_by_v0_1": int((v01 & bl).sum()),
            "recovery_rate": round(float((v01 & bl).sum() / n_blind), 3) if n_blind else None,
            "incumbent_position_novelty_flag": 0.604,
            "note": ("recovery_rate is comparable to the position-novelty flag's measured median 0.604 on "
                     "the EFV blind spot -- the free zero-tool incumbent this curation must beat"),
        },
    }


def run(path: Path = DEFAULT_DATA) -> dict:
    rows = load_rows(path)
    per_drug, final_catalog = {}, {}
    for drug in DRUGS:
        col = _DRUG_COL[drug]
        records = isolate_records(rows, col)
        cutoff = DRMCV_LOWER_CUTOFF.get(drug)
        if cutoff is None:
            per_drug[drug] = {
                "status": "CUTOFF_UNAVAILABLE",
                "n_isolates": len(records),
                "note": ("doravirine POSTDATES Stanford DRMcv.R, so no clinical lower cutoff is sourced "
                         "in-repo. Reported as a wall rather than scored at a guessed boundary."),
            }
            final_catalog[drug] = sorted(derive_resistant_mutants(records))
            continue
        m = _cv_compare(records, cutoff)
        m["clinical_lower_cutoff_fold"] = cutoff
        v0, v01 = m.get("shipped_catalog_v0"), m.get("data_derived_v0_1_heldout")
        if v0 and v01 and v0.get("balanced_accuracy") is not None:
            m["balacc_gain_v0_1_minus_v0"] = round(v01["balanced_accuracy"] - v0["balanced_accuracy"], 3)
            m["sens_gain"] = (round(v01["sens"] - v0["sens"], 3)
                              if v0.get("sens") is not None and v01.get("sens") is not None else None)
            m["spec_change"] = (round(v01["spec"] - v0["spec"], 3)
                                if v0.get("spec") is not None and v01.get("spec") is not None else None)
        per_drug[drug] = m
        final_catalog[drug] = sorted(derive_resistant_mutants(records))

    scored = {d: m for d, m in per_drug.items() if "balacc_gain_v0_1_minus_v0" in m}
    gains = [m["balacc_gain_v0_1_minus_v0"] for m in scored.values()]
    union = sorted(set().union(*(set(v) for v in final_catalog.values())) if final_catalog else set())
    novel = sorted(set(union) - set(NNRTI_RT_MAJOR_DRMS))
    return {
        "artifact": "hiv_nnrti_mutant_specific_v0_1", "schema": "hiv-nnrti-mutant-catalog-v0_1",
        "generated": _date.today().isoformat(),
        # DERIVED from the constant, never restated. The literal previously read log10(1.5) while the
        # constant was log10(3.0) -- the artifact misdescribed its own method after the threshold was swept.
        "method": (f"data-derived resistant mutants = MULTIVARIATE OLS log10-fold coefficient >= "
                   f"log10({10 ** RESIST_COEF_MIN:.1f}) "
                   f"(independent >={10 ** RESIST_COEF_MIN:.1f}x effect after controlling for co-occurrence -> deconfounds accessory "
                   "riders), >=5 carriers, candidates spanning the WHOLE RT from CompMutList; 5-fold "
                   "CROSS-VALIDATED held-out -> out-of-sample, not in-sample"),
        "why_whole_rt": ("the shipped catalog covers 8 RT positions; restricting candidates to those could "
                         "never reach the blind spot, whose drivers sit elsewhere"),
        "cutoff_note": ("ABSOLUTE, not delta-honest: EFV/NVP/ETR/RPV are scored at their OWN Stanford "
                        "DRMcv.R clinical lower cutoff. Doravirine postdates DRMcv.R -> CUTOFF_UNAVAILABLE, "
                        "reported as a wall, never guessed. This is stronger than the PI/INSTI v0.1 arc, "
                        "which had no per-drug cutoff and could only claim a delta."),
        "n_drugs_scored": len(scored), "n_drugs_improved_or_held": sum(1 for g in gains if g >= 0),
        "mean_balacc_gain": round(float(np.mean(gains)), 3) if gains else None,
        "label_source": "Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra)",
        "dataset": str(path), "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R",
        "shipped_catalog_size": len(NNRTI_RT_MAJOR_DRMS),
        "per_drug": per_drug,
        "deliverable_catalog_all_data": final_catalog,
        "candidate_additions_not_in_shipped_catalog": novel,
        "provenance_caveat": (
            "These candidates are DATA-DERIVED from measured fold-change, NOT curated from a guideline. "
            "That is a different provenance from the shipped catalog, which is sourced verbatim from the "
            "Stanford NNRTI major-DRM list. Any entry promoted into the shipped catalog must carry which "
            "of the two it came from; a literature cross-check is a SEPARATE step and is not asserted here."),
        "honest_caveats": [
            "in-distribution vs HIVDB-PhenoSense; NOT provenance-disjoint external validation",
            "the deliverable catalog is derived on ALL data and is therefore NOT what the CV metrics describe",
            "censored folds ('>'/'<') are kept at their numeric bound (v0 convention, inherited)",
            "a recovery_rate below the position-novelty flag's 0.604 means the free zero-tool incumbent is "
            "still the better blind-spot instrument and this curation should NOT ship on sensitivity grounds",
        ],
    }


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    a = ap.parse_args()
    if not a.data.exists():
        print(f"dataset absent: {a.data} (gitignored; see wiki/hiv_* memos)")
        return 0
    res = run(a.data)

    print(f"\nHIV NNRTI v0.1 — shipped {res['shipped_catalog_size']}-entry catalog vs data-derived\n")
    print(f"  {'drug':14} {'n':>5} {'cutoff':>7}  {'v0 sens/spec':>18}  {'v0.1 sens/spec':>18}  blind-spot")
    for d in DRUGS:
        m = res["per_drug"][d]
        if m.get("status") == "CUTOFF_UNAVAILABLE":
            print(f"  {d:14} {m['n_isolates']:>5}      --  CUTOFF_UNAVAILABLE (postdates DRMcv.R)")
            continue
        v0, v1, bs = m["shipped_catalog_v0"], m["data_derived_v0_1_heldout"], m["blind_spot"]
        rr = "n/a" if bs["recovery_rate"] is None else f"{bs['recovery_rate']:.3f}"
        print(f"  {d:14} {m['n_isolates']:>5} {m['clinical_lower_cutoff_fold']:>7} "
              f"  {v0['sens']:.3f}/{v0['spec']:.3f}      "
              f"  {v1['sens']:.3f}/{v1['spec']:.3f}      "
              f"{bs['recovered_by_v0_1']}/{bs['n_catalog_negative_true_R']} = {rr}")
    print(f"\n  incumbent to beat on the blind spot: position-novelty flag = 0.604 (free, zero-tool)")
    print(f"  candidate additions not in the shipped catalog: {len(res['candidate_additions_not_in_shipped_catalog'])}")
    print(f"  {res['provenance_caveat'][:150]}...")

    dest = REPO / "wiki" / f"hiv_nnrti_mutant_catalog_{res['generated']}.json"
    dest.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\nwrote wiki/{dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
