"""F3 Step 4/5: calibrate the phylogeny de-confounder against a DOUBLE CONTROL.

The deliverable of F3 is an INSTRUMENT, and it does not ship as calibrated until it is shown
to respond differentially to real structure. This harness is what decides that, and every
verdict below was named BEFORE the run.

CONTROL PAIR, SELECTED BY CONSTRUCTION
--------------------------------------
  negative control = the panel whose ancestry is randomized BY CONSTRUCTION (a bi-parental
                     segregant cross). A de-confounder must NOT destroy signal here: there
                     is no structure to remove.
  positive control = the panel drawn from a NATURAL population. A de-confounder must respond
                     MORE here than on the negative control.

Selection is from each panel's `design` field in `f3_substrate_census.PANELS`, declared before
any number was read. It is NOT by pc1_variance_fraction, which the census MEASURED to be unfit
for the job: pc1 conflates linkage disequilibrium with population structure, and a bi-parental
cross maximises LD by construction, so the plan's original "highest pc1" rule ranked the
NEGATIVE control first (excess 7.91 vs the natural population's 2.00). See the census artifact's
`pc1_selection_verdict`.

WHY THE CRITERION IS DIFFERENTIAL AND NOT "THE POSITIVE CONTROL FALLS TO ITS OWN NULL"
--------------------------------------------------------------------------------------
A natural panel may carry structure AND real functional signal, and a CORRECT de-confounder
removes the former while leaving the latter. An absolute "falls to its null" bar would
therefore condemn a working instrument. What a structure-de-confounder must do is respond MORE
where structure exists than where it is randomized by construction -- so the criterion is the
CONTRAST between the two controls, which is what makes them a pair rather than two runs.

This is instrument CALIBRATION against known answers, which is a legitimate positive+negative
control pair. It is NOT the same act as citing an experiment whose outcome was foregone: both
adverse verdicts below are genuinely reachable.

Run: uv run python scripts/f3_deconfound_calibration.py
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

N_SPLITS = 5
N_PERM = 100
MIN_TRAITS_FOR_SIGN_TEST = 5          # asserted, not derived
MAX_TRAITS = 12                       # bound the run; traits are taken in file order


# ---------------------------------------------------------------- frozen verdict machinery

VERDICTS = (
    "blocked:no_intact_panel",
    "blocked:no_negative_control",
    "blocked:no_positive_control",
    "blocked:positive_control_not_fold_compatible",
    "INDETERMINATE_UNDERPOWERED",
    "INVALID_OVER_AGGRESSIVE",
    "INVALID_INSENSITIVE",
    "VALID_DIFFERENTIAL_RESPONSE",
)

VERDICT_MEANING = {
    "blocked:no_intact_panel": "the census found no intact panel; nothing can be calibrated",
    "blocked:no_negative_control": "no bi-parental cross is cached, so over-aggression is untestable",
    "blocked:no_positive_control": "no natural-population panel is cached, so sensitivity is untestable",
    "blocked:positive_control_not_fold_compatible":
        "the natural-population panel has no grouping that is BOTH fold-compatible "
        "(largest cluster <= 1/k) and not mostly singletons, so grouped CV cannot be run on "
        "it at all. The instrument is built and NOT calibrated. Unblock = a natural-population "
        "panel whose relatedness is not concentrated in one clump (a USER decision: it means "
        "acquiring a substrate, e.g. the 1001-Genomes genotype matrix).",
    "INDETERMINATE_UNDERPOWERED": "a control has too few traits for a sign test",
    "INVALID_OVER_AGGRESSIVE":
        "the negative control STOPS beating its own grouped null: the de-confounder destroys "
        "signal where there is no structure to remove",
    "INVALID_INSENSITIVE":
        "the positive control's drop does NOT exceed the negative control's: the de-confounder "
        "does not respond differentially to measured structure. NOTE: this may be the GROUPING's "
        "blindness rather than the instrument's insensitivity -- IBS-greedy-rep is one grouping "
        "of many, and this verdict does not distinguish them.",
    "VALID_DIFFERENTIAL_RESPONSE":
        "the negative control keeps beating its own grouped null AND the positive control's "
        "drop exceeds it",
}


def compute_verdict(*, n_intact: int, has_negative: bool, has_positive: bool,
                    positive_fold_compatible: bool,
                    neg: dict | None, pos: dict | None) -> str:
    """Pure, ORDERED, first match wins. Frozen before the run."""
    if n_intact == 0:
        return "blocked:no_intact_panel"
    if not has_negative:
        return "blocked:no_negative_control"
    if not has_positive:
        return "blocked:no_positive_control"
    if not positive_fold_compatible:
        return "blocked:positive_control_not_fold_compatible"
    if (neg is None or pos is None
            or neg["n_traits"] < MIN_TRAITS_FOR_SIGN_TEST
            or pos["n_traits"] < MIN_TRAITS_FOR_SIGN_TEST):
        return "INDETERMINATE_UNDERPOWERED"
    if not neg["grouped_beats_own_null"]:
        return "INVALID_OVER_AGGRESSIVE"
    if not (pos["median_drop"] > neg["median_drop"]):
        return "INVALID_INSENSITIVE"
    return "VALID_DIFFERENTIAL_RESPONSE"


def sign_test(drops: list[float]) -> dict:
    """Two-sided exact sign test on per-trait drops. PAIRED per trait -- never a difference of
    aggregate medians, which is a different and weaker quantity."""
    from math import comb

    pos = sum(1 for d in drops if d > 0)
    neg = sum(1 for d in drops if d < 0)
    n = pos + neg
    if n == 0:
        return {"n_nonzero": 0, "n_positive": 0, "p": None}
    k = min(pos, neg)
    p = min(1.0, 2.0 * sum(comb(n, i) for i in range(k + 1)) / (2 ** n))
    return {"n_nonzero": n, "n_positive": pos, "p": p}


# ---------------------------------------------------------------- the arms

def _select_controls(panels) -> tuple[object | None, object | None]:
    neg = next((s for s in panels if s.design == "constructed_biparental"), None)
    pos = next((s for s in panels if s.design == "natural_population"), None)
    return neg, pos


def run_arm(spec, census_row: dict, *, max_traits: int = MAX_TRAITS,
            n_splits: int = N_SPLITS, n_perm: int = N_PERM, seed: int = 0) -> dict:
    """Random-CV vs grouped-CV on IDENTICAL rows, per trait. Each run carries its own
    permutation null under its OWN fold structure (the `folds` seam does that by construction).
    """
    from dna_decode.eval.genomic_prediction import cv_ridge_gp
    from dna_decode.eval.panel_structure import panel_groups

    mod = importlib.import_module(spec.arm_module)
    if spec.kind == "npz":
        d = np.load(spec.npz_geno, allow_pickle=True)
        X = np.asarray(d["X"], dtype=float)
        ids = [str(s) for s in d["samples"]]
        pheno = json.loads(Path(spec.npz_pheno).read_text(encoding="utf-8"))
        traits = {k: np.asarray(v, dtype=float) for k, v in list(pheno.items())[:max_traits]} \
            if isinstance(pheno, dict) else {}
    else:
        from scripts.f3_substrate_census import _unpack_loaded
        X, ids = _unpack_loaded(mod.load_genotype(str(mod._GENO), 1))
        ph = mod.load_phenotype(str(mod._PHENO))
        traits = _phenotype_traits(ph, len(ids), max_traits)

    pg = panel_groups(X, ids, max_largest=1.0 / n_splits)
    rows = []
    for name, y in traits.items():
        if y.shape[0] != X.shape[0] or np.sum(np.isfinite(y)) < n_splits * 4:
            continue
        try:
            a = cv_ridge_gp(X, y, trait=name, n_splits=n_splits, seed=seed, n_perm=n_perm)
            b = cv_ridge_gp(X, y, trait=name, n_splits=n_splits, seed=seed, n_perm=n_perm,
                            groups=pg.groups)
        except Exception as e:                        # noqa: BLE001
            rows.append({"trait": name, "error": f"{type(e).__name__}: {e}"})
            continue
        rows.append({
            "trait": name, "n": a.n,
            "r_random": round(a.predictive_r, 4), "r_grouped": round(b.predictive_r, 4),
            "drop": round(a.predictive_r - b.predictive_r, 4),
            "random_beats_own_null": bool(a.beats_null),
            "grouped_beats_own_null": bool(b.beats_null),
            "null_p95_random": round(a.null_r_p95, 4),
            "null_p95_grouped": round(b.null_r_p95, 4),
        })
    scored = [r for r in rows if "drop" in r]
    drops = [r["drop"] for r in scored]
    return {
        "panel": spec.name, "design": spec.design,
        "grouping": pg.as_dict(),
        "n_traits": len(scored),
        "per_trait": rows,
        "median_drop": (float(np.median(drops)) if drops else None),
        "n_dropped": sum(1 for d in drops if d > 0),
        "sign_test": sign_test(drops),
        # the negative control's absolute condition: grouping must not destroy the signal
        "grouped_beats_own_null": (
            sum(1 for r in scored if r["grouped_beats_own_null"]) > len(scored) / 2
            if scored else False),
    }


def _phenotype_traits(ph, n_rows: int, max_traits: int) -> dict:
    """Accept the shapes the arms' `load_phenotype` actually return."""
    if isinstance(ph, dict):
        out = {}
        for k, v in ph.items():
            a = np.asarray(v, dtype=float)
            if a.ndim == 1 and a.shape[0] == n_rows:
                out[str(k)] = a
            if len(out) >= max_traits:
                break
        return out
    if isinstance(ph, tuple):
        names, M = None, None
        for part in ph:
            a = np.asarray(part, dtype=object)
            if a.ndim == 2:
                M = np.asarray(part, dtype=float)
            elif a.ndim == 1 and M is None:
                names = [str(x) for x in part]
        if M is None:
            return {}
        if M.shape[0] != n_rows and M.shape[1] == n_rows:
            M = M.T
        cols = range(min(max_traits, M.shape[1]))
        if names and len(names) == M.shape[1]:
            return {names[j]: M[:, j] for j in cols}
        return {f"trait{j}": M[:, j] for j in cols}
    a = np.asarray(ph, dtype=float)
    if a.ndim == 2:
        if a.shape[0] != n_rows and a.shape[1] == n_rows:
            a = a.T
        return {f"trait{j}": a[:, j] for j in range(min(max_traits, a.shape[1]))}
    return {}


def main(argv=None) -> int:
    from scripts.f3_substrate_census import PANELS

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--census", type=Path,
                    default=REPO / "wiki" / f"f3_substrate_census_{date.today().isoformat()}.json")
    ap.add_argument("--out", type=Path,
                    default=REPO / "wiki" / f"f3_deconfound_calibration_{date.today().isoformat()}.json")
    ap.add_argument("--max-traits", type=int, default=MAX_TRAITS)
    a = ap.parse_args(argv)

    if not a.census.exists():
        hits = sorted((REPO / "wiki").glob("f3_substrate_census_*.json"))
        if not hits:
            print("no census artifact; run scripts/f3_substrate_census.py first")
            return 2
        a.census = hits[-1]
    cen = json.loads(a.census.read_text(encoding="utf-8"))
    by_name = {r["panel"]: r for r in cen["panels"]}

    neg_spec, pos_spec = _select_controls(PANELS)
    neg_row = by_name.get(neg_spec.name) if neg_spec else None
    pos_row = by_name.get(pos_spec.name) if pos_spec else None

    def _ok(row):
        return bool(row) and row.get("integrity") == "intact"

    def _foldable(row):
        return bool((row or {}).get("fold_compatible_k5", {}).get("compatible"))

    has_neg, has_pos = _ok(neg_row), _ok(pos_row)
    pos_foldable = _foldable(pos_row)

    neg_res = pos_res = None
    # Only run the arms when nothing upstream already blocks: a blocked verdict must emit NO
    # metrics, so computing them would be wasted and tempting to quote.
    pre = compute_verdict(n_intact=cen["n_intact"], has_negative=has_neg, has_positive=has_pos,
                          positive_fold_compatible=pos_foldable, neg=None, pos=None)
    if not pre.startswith("blocked:"):
        neg_res = run_arm(neg_spec, neg_row, max_traits=a.max_traits)
        pos_res = run_arm(pos_spec, pos_row, max_traits=a.max_traits)

    verdict = compute_verdict(n_intact=cen["n_intact"], has_negative=has_neg,
                              has_positive=has_pos, positive_fold_compatible=pos_foldable,
                              neg=neg_res, pos=pos_res)

    out = {
        "schema": "f3-deconfound-calibration-v1",
        "generated": date.today().isoformat(),
        "census": a.census.name,
        "verdict": verdict,
        "verdict_meaning": VERDICT_MEANING[verdict],
        "all_possible_verdicts": list(VERDICTS),
        "control_selection": {
            "rule": "BY CONSTRUCTION from each panel's `design` field, declared before any "
                    "number was read; NOT by pc1 (see the census pc1_selection_verdict)",
            "negative_control": (neg_spec.name if neg_spec else None),
            "negative_basis": (neg_spec.design_basis if neg_spec else None),
            "positive_control": (pos_spec.name if pos_spec else None),
            "positive_basis": (pos_spec.design_basis if pos_spec else None),
            "positive_fold_compatible_k5": pos_foldable,
            "positive_fold_incompatibility_reason":
                (pos_row or {}).get("fold_compatible_k5", {}).get("reason"),
        },
        "criterion": "DIFFERENTIAL: negative control keeps beating its own grouped null AND "
                     "the positive control's paired per-trait drop exceeds the negative "
                     "control's. Asserts no magnitude.",
        "honest_limits": [
            "The two controls have DIFFERENT traits, so the between-control comparison is "
            "two-sample, not paired; the within-control sign tests are paired.",
            "The positive control's structure is MEASURED while the negative control's "
            "randomization is BY CONSTRUCTION -- a deliberate asymmetry in evidence strength.",
            "IBS-greedy-representative is ONE grouping. An INVALID_INSENSITIVE verdict could "
            "be the grouping's blindness rather than the instrument's insensitivity.",
            f"n_perm={N_PERM} bounds the resolution of every beats_null decision.",
            "MIN_TRAITS_FOR_SIGN_TEST and MAX_TRAITS are ASSERTED, not derived.",
        ],
    }
    if not verdict.startswith("blocked:"):
        out["negative_control_result"] = neg_res
        out["positive_control_result"] = pos_res
    else:
        out["metrics_withheld"] = (
            "This run is blocked, so NO calibration metrics are emitted. A blocked verdict is "
            "a pre-registered outcome of this plan, not a partial result to be read as one."
        )

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"F3 calibration -> {a.out}")
    print(f"  VERDICT: {verdict}")
    print(f"  {VERDICT_MEANING[verdict]}")
    print(f"  negative control = {out['control_selection']['negative_control']}")
    print(f"  positive control = {out['control_selection']['positive_control']} "
          f"(fold-compatible k=5: {pos_foldable})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
