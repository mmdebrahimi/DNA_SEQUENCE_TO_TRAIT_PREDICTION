"""HIV quantitative calibration — Family B of the genome world-model plan (2026-07-11).

QUESTION: the deployed decoder emits binary R/S; Family A showed the additive model predicts log10
fold-change well (Spearman 0.75-0.91). Can we turn that point estimate into a CALIBRATED quantitative
output — a fold-change PREDICTION INTERVAL with honest, verified coverage? I.e. does a 90% interval
actually contain the true fold 90% of the time on held-out isolates?

METHOD (split-conformal — coverage-valid by construction, distribution-free):
  * Reuse Family A's harness (`hiv_epistasis`): additive binary mutation-presence features + a nested-CV
    ElasticNet whose OUT-OF-FOLD (OOF) predictions are honest held-out point estimates of log10 fold.
  * Conformal residuals r = |y - oof| are held-out (each oof is from a model that didn't see that isolate).
  * SPLIT-CONFORMAL: shuffle the residuals, use HALF as the calibration set -> the interval half-width is the
    finite-sample conformal quantile q = the ceil((m+1)(1-alpha))/m empirical quantile of the calib |r|; test
    empirical coverage on the OTHER half (never used to set q). Repeat over REPEATS shuffles, average.
  * A calibrated interval = oof +/- q. Report held-out coverage at target 0.90 and 0.80, and the interval
    HALF-WIDTH in FOLD units (10^q) = the interpretable "fold x/÷ W" payload.

PRE-REGISTERED BAR: CALIBRATED_INTERVALS iff |coverage_90 - 0.90| <= COVER_TOL (=0.05) on >= PASS_FRACTION
(=0.5) of the powered drug-cells (n>=N_MIN, both-fold spread). Else MISCALIBRATED (over/under-covering).

Honest scope: split-conformal gives MARGINAL coverage (not conditional on the mutation profile); censored
folds ('>'/'<') kept at the bound bias the tail slightly (disclosed). This is Family A's per-cell R² turned
into a coverage-valid interval — the natural depth follow-on. NOT a sequence embedding. Frozen surface READ-only.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import hiv_epistasis as he  # noqa: E402  (reuse the Family A harness)
from scripts.hiv_nnrti_validate import _parse_fold  # noqa: E402

COVER_TOL = 0.05
PASS_FRACTION = 0.5
REPEATS = 20
TARGETS = (0.90, 0.80)
SEED = 0


def _conformal_q(calib_abs_res: np.ndarray, alpha: float) -> float:
    """Finite-sample split-conformal quantile: ceil((m+1)(1-alpha))/m empirical quantile of |residuals|."""
    m = len(calib_abs_res)
    if m == 0:
        return float("nan")
    k = math.ceil((m + 1) * (1 - alpha))
    if k > m:                              # not enough calib points for this coverage -> widest
        return float(np.max(calib_abs_res))
    return float(np.sort(calib_abs_res)[k - 1])


def calibrate_drug(y: np.ndarray, oof: np.ndarray, repeats=REPEATS, seed=SEED):
    """Split-conformal held-out coverage at each TARGET; averaged over `repeats` shuffles."""
    res = np.abs(y - oof)
    n = len(res)
    rng = np.random.default_rng(seed)
    out = {}
    for target in TARGETS:
        alpha = 1 - target
        covs, qs = [], []
        for _ in range(repeats):
            idx = rng.permutation(n)
            half = n // 2
            calib, test = idx[:half], idx[half:]
            q = _conformal_q(res[calib], alpha)
            covs.append(float(np.mean(res[test] <= q)))
            qs.append(q)
        cov = float(np.mean(covs)); q = float(np.mean(qs))
        out[f"cover_{int(target*100)}"] = round(cov, 4)
        out[f"halfwidth_log10_{int(target*100)}"] = round(q, 4)
        out[f"fold_factor_{int(target*100)}"] = round(10 ** q, 2)   # "true fold within x/÷ this factor"
    out["calibrated_90"] = bool(abs(out["cover_90"] - 0.90) <= COVER_TOL)
    return out


def run_all(data_dir: Path, classes=None, seed=SEED):
    classes = classes or list(he._CLASSES)
    per_class = {}
    for cname in classes:
        p = data_dir / he._CLASSES[cname]
        if not p.exists():
            per_class[cname] = {"class": cname, "note": f"absent {p}"}
            continue
        rows, drug_cols, pos_cols = he.load_class(p)
        X_all, feat_names, _ = he.build_presence(rows, pos_cols)
        per_drug = {}
        for dc in drug_cols:
            keep = [i for i, r in enumerate(rows)
                    if (_parse_fold(r.get(dc, "")) or 0) > 0]
            if len(keep) < he.N_MIN:
                per_drug[dc] = {"n": len(keep), "powered": False, "note": "too few"}
                continue
            idx = np.array(keep)
            y = np.array([math.log10(_parse_fold(rows[i][dc])) for i in keep], float)
            if float(np.var(y)) < 1e-6:
                per_drug[dc] = {"n": len(keep), "powered": False, "note": "degenerate"}
                continue
            oof = he.nested_oof(X_all[idx], y, seed)
            cal = calibrate_drug(y, oof, seed=seed)
            per_drug[dc] = {"n": len(keep), "powered": True,
                            "r2_oof": round(he._r2(y, oof), 4), **cal}
        per_class[cname] = {"class": cname, "per_drug": per_drug}

    powered = [(c["class"], dc, m) for c in per_class.values()
               for dc, m in c.get("per_drug", {}).items() if m.get("powered")]
    n_pow = len(powered)
    n_cal = sum(1 for _, _, m in powered if m["calibrated_90"])
    frac = (n_cal / n_pow) if n_pow else 0.0
    verdict = ("CALIBRATED_INTERVALS" if (n_pow and frac >= PASS_FRACTION)
               else ("MISCALIBRATED" if n_pow else "NO_POWERED_CELLS"))
    return {
        "artifact": "hiv_quantitative_calibration",
        "schema": "hiv-quantitative-calibration-v1",
        "question": "Does a split-conformal prediction interval on HIV fold-change achieve its nominal "
                    "held-out coverage (is the quantitative decoder honestly calibrated)?",
        "method": "split-conformal on the additive nested-CV OOF residuals (Family A harness); finite-sample "
                  "quantile; coverage tested on a held-out residual half; averaged over %d shuffles" % REPEATS,
        "label_source": "Stanford HIVDB PhenoSense fold-change (free, independent wet-lab)",
        "prereg": {"COVER_TOL": COVER_TOL, "PASS_FRACTION": PASS_FRACTION, "REPEATS": REPEATS,
                   "targets": list(TARGETS), "N_MIN": he.N_MIN, "seed": seed},
        "verdict": verdict, "n_powered": n_pow, "n_calibrated_90": n_cal, "fraction_calibrated": round(frac, 3),
        "honest_caveats": [
            "Split-conformal gives MARGINAL coverage (not conditional on the mutation profile).",
            "Censored folds ('>'/'<') kept at the numeric bound bias the tail slightly.",
            "Reuses Family A's additive OOF (no interaction terms — A showed they don't help rank; the point "
            "estimate is the additive model's).",
        ],
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; conformal per Vovk/Lei split-conformal",
        "per_class": per_class,
    }


def render_md(res, generated):
    L = [f"# HIV quantitative calibration — are the fold-change prediction intervals honest? ({generated})", "",
         f"**Verdict: {res['verdict']}** — {res['n_calibrated_90']}/{res['n_powered']} powered drug-cells have "
         f"a 90% interval within {res['prereg']['COVER_TOL']} of nominal coverage (fraction "
         f"{res['fraction_calibrated']}; PASS bar {res['prereg']['PASS_FRACTION']}).", "",
         f"{res['question']} Label = {res['label_source']}. {res['method']}.", "",
         "`cover_90` = held-out fraction of isolates whose true fold falls in the 90% interval (target 0.90). "
         "`fold_factor_90` = the interval half-width in FOLD units (true fold within x/÷ this factor).", "",
         "| class:drug | n | R2(oof) | **cover_90** | fold± (90%) | cover_80 | calibrated |",
         "|---|---|---|---|---|---|---|"]
    for c in res["per_class"].values():
        for dc, m in c.get("per_drug", {}).items():
            if not m.get("powered"):
                L.append(f"| {c['class']}:{dc} | {m.get('n')} | — | {m.get('note','')} | — | — | — |")
                continue
            L.append(f"| {c['class']}:{dc} | {m['n']} | {m['r2_oof']} | **{m['cover_90']}** | "
                     f"x/÷{m['fold_factor_90']} | {m['cover_80']} | {'YES' if m['calibrated_90'] else 'no'} |")
    L += ["", "## Honest caveats"] + [f"- {x}" for x in res["honest_caveats"]]
    L += ["", f"Citation: {res['citation']}."]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=REPO / "data" / "raw" / "hiv")
    ap.add_argument("--classes", default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if not any((a.data_dir / f).exists() for f in he._CLASSES.values()):
        print(f"ERROR: no HIV DataSets under {a.data_dir} (gitignored)", file=sys.stderr)
        return 2
    classes = [c.strip().upper() for c in a.classes.split(",")] if a.classes else None
    today = _date.today().isoformat()
    res = run_all(a.data_dir, classes=classes)
    out = a.out or (REPO / "wiki" / f"hiv_quantitative_calibration_{today}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
