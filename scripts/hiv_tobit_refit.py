"""Does censoring-aware fitting change which mutants the v0.1 catalogs select?

THE QUESTION. `wiki/hiv_fold_censoring_audit_2026-09-22.md` measured that the Stanford PhenoSense fold
tables are right-censored at 100 (nothing above it anywhere; 3TC 45% and NVP 33% of observations at the
cap). The NRTI / PI / INSTI v0.1 catalogs were each selected by thresholding a multivariate-OLS
log10-fold coefficient at >= log10(1.5) -- fit on that censored response. The audit deliberately did NOT
claim any published figure is wrong; it named this as the decidable follow-on. This is that follow-on.

WHAT IT ASKS, AND THE NARROWER THING IT DOES NOT. Only whether the SELECTION DECISION moves: does any
mutant cross RESIST_COEF_MIN when the same design is refit with a Tobit (censored-normal) likelihood?
It does NOT re-run the 5-fold CV, so it produces NO new balanced-accuracy gain figures. Re-deriving
those is a separate job, gated on a MATERIAL verdict.

THE NEGATIVE CONTROL IS SUPPLIED BY THE DATA AND IS LOAD-BEARING. Tobit REDUCES to OLS as the censored
fraction goes to zero -- there is nothing to correct. Five drugs here have <2% censoring (tenofovir and
bictegravir have literally none), so if their coefficients move, the TOBIT IMPLEMENTATION is wrong, not
the finding. That check runs first and can veto the whole run (`IMPLEMENTATION_SUSPECT`).

LIKE-FOR-LIKE BY CONSTRUCTION. The candidate set, design matrix, MIN_CARRIERS and RESIST_COEF_MIN are
IMPORTED from the three shipped catalog builders rather than restated, so OLS and Tobit are compared on
ONE design. A restated constant drifts from the thing it mirrors and then the comparison is not the one
the catalogs were built on.

Bar: wiki/hiv_tobit_refit_acceptance_bar.json (frozen AFTER the mechanism analysis, before any fit).
Offline, CPU, seconds. NO GPU is needed or used -- this is a censored-regression MLE on ~2,000 rows.
Frozen AMR surface untouched; hiv_amr.py unmodified; no deployed catalog is changed by any verdict.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from sklearn.linear_model import LinearRegression

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.hiv_nrti_mutant_catalog import (  # noqa: E402
    MIN_CARRIERS, RESIST_COEF_MIN, _isolate_records as _nrti_records,
)
from scripts.hiv_nrti_validate import DEFAULT_DATA as NRTI_DATA, NRTI_DRUGS, _NRTI_COL  # noqa: E402
from scripts.hiv_nnrti_validate import load_rows  # noqa: E402
from scripts.hiv_pi_mutant_catalog import _isolate_records as _ts_records  # noqa: E402
from scripts.hiv_targetsite_validate import _CLASS_SPECS  # noqa: E402

CENSOR_LOG10 = math.log10(100.0)      # the measured ceiling, on the scale the catalogs are fit on
CONTROL_MAX_CENSORED = 0.02           # a drug below this is a Tobit-reduces-to-OLS control
CONTROL_MAX_DRIFT = 0.01              # ... and its coefficients must not move by more than this
HIGH_CENSORING = 0.10                 # a drug at or above this is where an effect is possible

MATERIAL, NEGLIGIBLE, SUSPECT, INDETERMINATE = (
    "MATERIAL", "NEGLIGIBLE", "IMPLEMENTATION_SUSPECT", "INDETERMINATE")


def build_design(records) -> tuple[np.ndarray, np.ndarray, list[str]] | tuple[None, None, None]:
    """The EXACT design the shipped catalog builders fit, rebuilt once so both models share it.

    Mirrors `derive_resistant_mutants`: candidates are mutants with >= MIN_CARRIERS carriers, columns
    are binary presence, the response is log10(fold). Kept here (rather than imported) only because the
    shipped function returns the thresholded SET and never exposes the matrix -- the constants it uses
    ARE imported, so the two cannot drift on anything that decides membership.
    """
    counts: dict[str, int] = {}
    for obs, _ in records:
        for s in obs:
            counts[s] = counts.get(s, 0) + 1
    cands = sorted(m for m, n in counts.items() if n >= MIN_CARRIERS)
    if not cands or len(records) < 30:
        return None, None, None
    fidx = {m: j for j, m in enumerate(cands)}
    X = np.zeros((len(records), len(cands)))
    y = np.zeros(len(records))
    for i, (obs, fold) in enumerate(records):
        y[i] = math.log10(fold)
        for s in obs:
            if s in fidx:
                X[i, fidx[s]] = 1.0
    return X, y, cands


def fit_tobit(X: np.ndarray, y: np.ndarray, censor: float = CENSOR_LOG10) -> dict:
    """Right-censored normal regression (Tobit type I) by maximum likelihood.

    An observation at or above `censor` is known only to be >= censor, so it contributes the SURVIVAL
    probability rather than a density. OLS instead treats it as exactly `censor`, which understates
    every such response and therefore biases downward the coefficients of the predictors that drive the
    response up -- the attenuation this run exists to measure.

    Optimised in (beta, log sigma) so sigma stays positive without a constraint, started from the OLS
    solution, using `logpdf`/`logsf` rather than logs of densities so the censored tail does not
    underflow to -inf and silently flatten the objective.
    """
    Xc = np.hstack([np.ones((X.shape[0], 1)), X])          # explicit intercept
    cens = y >= censor - 1e-9
    ols = LinearRegression(fit_intercept=True).fit(X, y)
    resid = y - ols.predict(X)
    start = np.concatenate([[ols.intercept_], ols.coef_,
                            [math.log(max(resid.std(ddof=1), 1e-3))]])

    def nll(theta):
        beta, log_s = theta[:-1], theta[-1]
        s = math.exp(min(log_s, 20.0))
        mu = Xc @ beta
        z = (y - mu) / s
        ll = np.where(cens,
                      norm.logsf((censor - mu) / s),
                      norm.logpdf(z) - log_s)
        if not np.all(np.isfinite(ll)):
            return 1e12
        return -float(ll.sum())

    res = minimize(nll, start, method="L-BFGS-B",
                   options={"maxiter": 20000, "maxfun": 40000, "ftol": 1e-12})
    beta = res.x[:-1]
    return {
        "converged": bool(res.success),
        "message": str(res.message),
        "intercept": float(beta[0]),
        "coef": beta[1:],
        "sigma": float(math.exp(res.x[-1])),
        "n_censored": int(cens.sum()),
        "censored_fraction": float(cens.mean()),
    }


def compare_drug(drug: str, cls: str, records) -> dict:
    X, y, cands = build_design(records)
    if X is None:
        return {"drug": drug, "class": cls, "status": "too_few_records"}
    ols_coef = LinearRegression().fit(X, y).coef_
    tob = fit_tobit(X, y)
    t_coef = tob["coef"]

    ols_set = {cands[j] for j in range(len(cands)) if ols_coef[j] >= RESIST_COEF_MIN}
    tob_set = {cands[j] for j in range(len(cands)) if t_coef[j] >= RESIST_COEF_MIN}
    delta = t_coef - ols_coef

    # KNIFE-EDGE MARGIN. A mutant whose OLS coefficient already sat a hair under the threshold flips on
    # a rounding-scale nudge -- that is a coin toss, not evidence that censoring reached the decision.
    # Without this a change of 0.0002 on a coefficient of 0.1759 reads identically to a real one.
    flips = []
    for j, m in enumerate(cands):
        if (m in tob_set) != (m in ols_set):
            flips.append({
                "mutant": m,
                "direction": "gained" if m in tob_set else "lost",
                "ols_coef": round(float(ols_coef[j]), 5),
                "tobit_coef": round(float(t_coef[j]), 5),
                "ols_margin_to_threshold": round(float(abs(ols_coef[j] - RESIST_COEF_MIN)), 5),
                "knife_edge": bool(abs(ols_coef[j] - RESIST_COEF_MIN) < 0.01),
            })
    return {
        "drug": drug, "class": cls, "status": "ok",
        "n_records": int(len(y)),
        "n_candidates": len(cands),
        "censored_fraction": round(tob["censored_fraction"], 4),
        "n_censored": tob["n_censored"],
        "tobit_converged": tob["converged"],
        "tobit_message": tob["message"],
        "mean_coef_delta": round(float(delta.mean()), 5),
        "max_abs_coef_delta": round(float(np.abs(delta).max()), 5),
        "catalog_ols_n": len(ols_set),
        "catalog_tobit_n": len(tob_set),
        "gained_under_tobit": sorted(tob_set - ols_set),
        "lost_under_tobit": sorted(ols_set - tob_set),
        "membership_changed": bool(ols_set != tob_set),
        "flips": flips,
        "n_flips_knife_edge": sum(1 for f in flips if f["knife_edge"]),
        "n_flips_decisive": sum(1 for f in flips if not f["knife_edge"]),
    }


def verdict(per_drug: list[dict]) -> dict:
    """The FROZEN rule, applied mechanically. The control is evaluated FIRST and can veto."""
    ok = [d for d in per_drug if d.get("status") == "ok"]
    if not ok:
        return {"verdict": INDETERMINATE, "reason": "no drug produced a fit"}

    controls = [d for d in ok if d["censored_fraction"] < CONTROL_MAX_CENSORED]
    failed = [d["drug"] for d in controls if d["max_abs_coef_delta"] >= CONTROL_MAX_DRIFT]
    if failed:
        return {"verdict": SUSPECT,
                "reason": (f"Tobit must reduce to OLS at ~zero censoring, but {failed} moved by "
                           f">= {CONTROL_MAX_DRIFT}. The implementation is suspect; interpret nothing."),
                "control_drugs": [d["drug"] for d in controls], "control_failures": failed}

    high = [d for d in ok if d["censored_fraction"] >= HIGH_CENSORING]
    unconverged = [d["drug"] for d in high if not d["tobit_converged"]]
    if unconverged:
        return {"verdict": INDETERMINATE,
                "reason": f"Tobit did not converge on high-censoring drug(s) {unconverged}",
                "control_drugs": [d["drug"] for d in controls]}

    changed = [d["drug"] for d in ok if d["membership_changed"]]
    return {
        "verdict": MATERIAL if changed else NEGLIGIBLE,
        "reason": (f"catalog membership changes on {changed}" if changed else
                   "the negative control holds and NO drug's catalog membership changes at "
                   "RESIST_COEF_MIN -- the censoring is real but did not reach the selection decision"),
        "control_drugs": [d["drug"] for d in controls],
        "control_max_drift": round(max((d["max_abs_coef_delta"] for d in controls), default=0.0), 6),
        "drugs_with_membership_change": changed,
        "high_censoring_drugs": [d["drug"] for d in high],
        "directional_check": {
            d["drug"]: {"censored_fraction": d["censored_fraction"],
                        "mean_coef_delta": d["mean_coef_delta"],
                        "gained": len(d["gained_under_tobit"]), "lost": len(d["lost_under_tobit"])}
            for d in high},
    }


def collect(data_dir: Path) -> list[dict]:
    out = []
    rows_n = load_rows(data_dir / NRTI_DATA.name)
    for d in NRTI_DRUGS:
        out.append(compare_drug(d, "NRTI", _nrti_records(rows_n, _NRTI_COL[d])))
    for cls_key, fname in (("PI", "PI_DataSet.txt"), ("INSTI", "INI_DataSet.txt")):
        cls, _, drug_cols = _CLASS_SPECS[cls_key]
        rows = load_rows(data_dir / fname)
        for d, col in drug_cols.items():
            out.append(compare_drug(d, cls_key, _ts_records(rows, col, cls)))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=REPO / "data" / "raw" / "hiv")
    ap.add_argument("--out", type=Path,
                    default=REPO / "wiki" / f"hiv_tobit_refit_{_date.today()}.json")
    a = ap.parse_args(argv)

    if not (a.data_dir / "PI_DataSet.txt").exists():
        print("HIV datasets not present (they are gitignored); nothing to refit")
        return 2

    per_drug = collect(a.data_dir)
    v = verdict(per_drug)
    doc = {
        "schema": "hiv-tobit-refit-v1",
        "date": str(_date.today()),
        "bar": "wiki/hiv_tobit_refit_acceptance_bar.json",
        "question": ("does refitting the v0.1 catalog design with a censoring-aware (Tobit) likelihood "
                     "change WHICH mutants clear RESIST_COEF_MIN?"),
        "censoring_point_log10": CENSOR_LOG10,
        "selection_threshold": RESIST_COEF_MIN,
        "constants_imported_not_restated": {
            "MIN_CARRIERS": MIN_CARRIERS, "RESIST_COEF_MIN": RESIST_COEF_MIN,
            "source": "scripts/hiv_nrti_mutant_catalog.py",
        },
        "results": per_drug,
        **v,
        "verdict_is_mechanical": "verdict() applied to the frozen rule; not authored",
        "honest_limits": [
            "Does NOT re-run the 5-fold CV and therefore produces NO new balanced-accuracy gain "
            "figures. It asks only whether the SELECTION DECISION moves.",
            "NOT a decoder change under any verdict: no v0.1 catalog is deployed -- hiv_amr.py routes "
            "NRTI/PI/INSTI through POSITION-based classes.",
            "The censoring point is INFERRED from a pile-up with an empty tail above it, not read off "
            "assay documentation.",
            "Tobit assumes normal homoscedastic latent errors on the log10 scale -- the same "
            "assumption OLS already makes here, not extra evidence for it.",
            "No GPU is used or needed; this is a censored-regression MLE on ~2,000 rows.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"{'class':<6}{'drug':<15}{'%cens':>7}{'meanD':>9}{'maxD':>9}{'OLS':>5}{'Tob':>5}  change")
    for d in per_drug:
        if d.get("status") != "ok":
            print(f"{d['class']:<6}{d['drug']:<15}  {d['status']}")
            continue
        ch = ("+" + ",".join(d["gained_under_tobit"]) if d["gained_under_tobit"] else "") + \
             (" -" + ",".join(d["lost_under_tobit"]) if d["lost_under_tobit"] else "")
        print(f"{d['class']:<6}{d['drug']:<15}{d['censored_fraction']*100:>6.1f}%"
              f"{d['mean_coef_delta']:>9.4f}{d['max_abs_coef_delta']:>9.4f}"
              f"{d['catalog_ols_n']:>5}{d['catalog_tobit_n']:>5}  {ch or '-'}")
    print(f"\ncontrol drugs (<2% censored): {v.get('control_drugs')}")
    print(f"control max drift: {v.get('control_max_drift')} (must be < {CONTROL_MAX_DRIFT})")
    print(f"\n-> {v['verdict']}: {v['reason']}\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
