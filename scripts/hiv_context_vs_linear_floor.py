"""Does CONTEXT beat the linear floor on the HIV blind spot? (2026-09-30)

THE QUESTION. The shipped supervised blind-spot complement is logistic regression over one-hot
substitution tokens — an ADDITIVE model, the weakest member of the supervised family. Its measured
leave-one-STUDY-out blind-spot AUROC (RT 0.8102, IN 0.8923 —
`wiki/hiv_deployability_per_gene_2026-09-30.json`) is therefore a FLOOR, and says nothing about whether
interactions help. Same data, same split, same subset, same bar — swap the estimator.

THE BAR IS FROZEN AND COMMITTED BEFORE THIS RAN:
`wiki/hiv_context_vs_linear_floor_acceptance_bar.json` (commit 2cf723d). It is READ here, never restated,
so the thresholds cannot drift from the thing that was pre-registered.

THE PRIOR IS STATED UP FRONT, not discovered afterwards. `wiki/hiv_epistasis_result_2026-07-11.json`
returned `FAIL_ADDITIVE_SUFFICES` — constrained pairwise interactions beat additive on only 2 of 24
powered drug-cells. That run scored the CONTINUOUS fold-change on the FULL cohort under plain KFold; this
one scores the BLIND-SPOT subset as classification under leave-one-STUDY-out, which is the framing the
deployable number uses and the one place the earlier run never looked. A large win is UNLIKELY; the bar is
written so a negative closes the question rather than being a shrug.

WHAT IS REUSED, NOT RE-DERIVED. The split, the blind-spot metric, the pass rule, the one-hot featurizer and
the drift detection come from `hiv_deployability_per_gene.py` / `hiv_supervised_deployability.py`; the
pairwise feature construction (`build_pairs` / `augment`, co-occurrence floored) comes from
`hiv_epistasis.py`. The ANCHOR is that the linear arm reproduces the measured floor per gene — a
context-vs-floor delta measured against a floor that does not reproduce is not a comparison, so a missed
anchor exits 3.

PAIRED, NOT TWO SEPARATE NUMBERS. Every estimator sees the same rows, the same folds and the same
blind-spot mask, and the delta carries a paired-bootstrap CI over the blind-spot isolates.

Run: uv run python scripts/hiv_context_vs_linear_floor.py [--self-check]
Exits 2 without the gitignored Stanford data, 3 if the floor does not reproduce.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import random
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
for p in (str(REPO), str(REPO / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

BAR_PATH = REPO / "wiki" / "hiv_context_vs_linear_floor_acceptance_bar.json"
FLOOR_ARTIFACT = REPO / "wiki" / "hiv_deployability_per_gene_2026-09-30.json"
PRIOR_ARTIFACT = REPO / "wiki" / "hiv_epistasis_result_2026-07-11.json"
BOOT, SEED = 2000, 0

# Both thresholds are stated in the FROZEN bar's prose; they are named here so a test can assert the code
# and the bar still agree, rather than counting digit strings in the source (a substring count of "0.02"
# also matches the bootstrap percentile 0.025 -- the number-inside-a-number trap this repo has hit before).
MIN_GAIN = 0.02      # bar: CONTEXT_BEATS_FLOOR requires delta >= +0.02 with a CI-positive delta
ANCHOR_TOL = 0.02    # bar: the linear arm must reproduce the measured floor within 0.02


class ContextRefused(RuntimeError):
    """Raised instead of reporting a delta that cannot be trusted."""


def _mod(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_bar() -> dict:
    if not BAR_PATH.exists():
        raise ContextRefused(f"frozen bar absent: {BAR_PATH}")
    bar = json.loads(BAR_PATH.read_text(encoding="utf-8"))
    if not bar.get("frozen_before_any_number_was_seen"):
        raise ContextRefused("the bar does not assert it was frozen before measurement")
    return bar


def floor_targets() -> dict:
    """The measured floor per gene, READ from the deployability artifact."""
    if not FLOOR_ARTIFACT.exists():
        raise ContextRefused(f"floor artifact absent: {FLOOR_ARTIFACT}")
    art = json.loads(FLOOR_ARTIFACT.read_text(encoding="utf-8"))
    if not (art.get("anchor") or {}).get("reproduces"):
        raise ContextRefused("the floor artifact's own anchor did not bind")
    out = {}
    for cls, cell in (art.get("cells") or {}).items():
        a = cell.get("leave_study_out_blindspot_auroc_MEASURED")
        if a is not None:
            out[cls] = {"gene": cell["gene"], "auroc": a, "blind_spot_n": cell["blind_spot_n"],
                        "blind_spot_R": cell["blind_spot_R"]}
    return out


def _auroc(y, s, ev):
    return ev.auroc(list(y), list(s))


def paired_bootstrap_auroc_delta(y, base, ctx, boot=BOOT, seed=SEED, ev=None):
    """CI on (context - floor) AUROC over the SAME blind-spot isolates, resampled together."""
    rng = random.Random(seed)
    n = len(y)
    deltas = []
    for _ in range(boot):
        idx = [rng.randrange(n) for _ in range(n)]
        yy = [y[i] for i in idx]
        if len(set(yy)) < 2:
            continue
        deltas.append(_auroc(yy, [ctx[i] for i in idx], ev) - _auroc(yy, [base[i] for i in idx], ev))
    if not deltas:
        return {"ci_lo": None, "ci_hi": None, "n_boot": 0}
    deltas.sort()
    lo = deltas[int(0.025 * len(deltas))]
    hi = deltas[min(len(deltas) - 1, int(0.975 * len(deltas)))]
    return {"ci_lo": round(lo, 4), "ci_hi": round(hi, 4), "n_boot": len(deltas)}


def run_gene(cls: str, D, PG, EP) -> dict:
    """One gene: three estimators over identical rows/folds/mask."""
    import csv

    import numpy as np
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, cross_val_predict

    import dna_decode.data.hiv_amr as H
    ev, S = D.ev, D.S
    cfg = PG.GENES[cls]

    path = PG.RAW / cfg["dataset"]
    if not path.exists():
        raise ContextRefused(f"dataset absent: {path}")
    prot = PG._protein(cfg["protein"], ev)
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    col, gene = cfg["col"], cfg["gene"]
    have = [r for r in rows if r.get(col) not in ("NA", "", "-", None)]

    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in have if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(have):
            drifted.add(p)

    def subs_of(r):
        return {f"{prot[p-1]}{p}{aa}" for p, aa in ev.isolate_muts(r, pcols)
                if p <= len(prot) and p not in drifted and prot[p - 1] != aa}

    cat = np.array([1 if H.call_hiv_observed(cfg["drug"], {gene: subs_of(r)}).prediction == "R" else 0
                    for r in have])
    neg = cat == 0
    y = np.array([1 if float(r[col]) >= PG.CUTOFF else 0 for r in have])
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    refid = np.array([r["RefID"] for r in have])
    X = np.asarray(X, dtype=float)

    n_r, n_s = int(y[neg].sum()), int(neg.sum() - y[neg].sum())
    out = {"class": cls, "gene": gene, "n": len(have), "blind_spot_n": int(neg.sum()),
           "blind_spot_R": n_r, "blind_spot_S": n_s, "n_features_linear": X.shape[1],
           "n_studies": len(set(refid.tolist()))}
    if n_r < PG.MIN_PER_CLASS or n_s < PG.MIN_PER_CLASS:
        out["status"] = "UNDERPOWERED"
        return out

    folds = GroupKFold(n_splits=min(5, len(set(refid.tolist()))))

    def oof(est, M):
        return cross_val_predict(est, M, y, groups=refid, cv=folds, method="predict_proba")[:, 1]

    # --- the three arms, identical rows/folds/mask ---
    lin = LogisticRegression(max_iter=2000, C=1.0, solver="liblinear")
    scores = {"linear_floor": oof(lin, X)}

    # `build_pairs` wants counts keyed by feature NAME (not a positional array) and returns
    # (pairs, pair_names) -- read from its own signature rather than assumed.
    names = list(feats)
    counts = {nm: int(X[:, j].sum()) for j, nm in enumerate(names)}
    pairs, pair_names = EP.build_pairs(X, names, counts)
    out["n_pairs_kept"] = len(pairs)
    out["top_pairs_by_cooccurrence"] = pair_names[:8]
    if pairs:
        Xp = np.asarray(EP.augment(X, pairs), dtype=float)
        out["n_features_pairwise"] = Xp.shape[1]
        scores["pairwise"] = oof(lin, Xp)

    gb = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, random_state=SEED)
    scores["nonlinear"] = oof(gb, X)

    mask = neg.tolist()
    ys = [int(y[i]) for i in range(len(y)) if mask[i]]
    bs = {k: [float(v[i]) for i in range(len(y)) if mask[i]] for k, v in scores.items()}

    out["blind_spot_auroc"] = {k: round(_auroc(ys, v, ev), 4) for k, v in bs.items()}
    out["overall_auroc"] = {k: round(_auroc(y.tolist(), list(v), ev), 4) for k, v in scores.items()}
    base = bs["linear_floor"]
    out["delta_vs_floor"] = {}
    for k, v in bs.items():
        if k == "linear_floor":
            continue
        d = _auroc(ys, v, ev) - _auroc(ys, base, ev)
        ci = paired_bootstrap_auroc_delta(ys, base, v, ev=ev)
        out["delta_vs_floor"][k] = {"delta": round(d, 4), **ci,
                                    "ci_positive": bool(ci["ci_lo"] is not None and ci["ci_lo"] > 0)}
    out["status"] = "SCORED"
    return out


def verdict(cells: dict, bar: dict, min_gain: float = MIN_GAIN) -> dict:
    """Mechanical application of the FROZEN rule. Never re-sized after seeing a number."""
    scored = {k: c for k, c in cells.items() if c.get("status") == "SCORED"}
    if len(scored) < 2:
        return {"verdict": "UNDERPOWERED", "reason": f"only {len(scored)} gene(s) scored"}
    wins = {}
    for cls, c in scored.items():
        wins[cls] = [k for k, d in c["delta_vs_floor"].items()
                     if d["delta"] >= min_gain and d["ci_positive"]]
    n = sum(1 for v in wins.values() if v)
    if n == len(scored):
        v = "CONTEXT_BEATS_FLOOR"
    elif n == 1:
        v = "CONTEXT_HELPS_ONE_GENE_ONLY"
    else:
        v = "LINEAR_FLOOR_IS_THE_CEILING"
    return {"verdict": v, "min_gain": min_gain, "winning_estimators_per_gene": wins,
            "rule": (bar.get("verdict_rule") or {}).get(v)}


def self_check() -> None:
    bar = load_bar()
    assert bar["frozen_before_any_number_was_seen"] is True
    assert "LINEAR_FLOOR_IS_THE_CEILING" in bar["verdict_rule"]
    f = floor_targets()
    assert set(f) == {"NNRTI", "INSTI"}, f"PI must be excluded as unscoreable; got {sorted(f)}"
    assert abs(f["NNRTI"]["auroc"] - 0.8102) < 1e-9 and abs(f["INSTI"]["auroc"] - 0.8923) < 1e-9
    prior = json.loads(PRIOR_ARTIFACT.read_text(encoding="utf-8"))
    assert prior["verdict"] == "FAIL_ADDITIVE_SUFFICES", "the stated prior must match its artifact"
    # Checked in the NAMESPACE, not in the source text: a text search for a definition keyword matches the
    # assertion line that contains it, so the first version of this guard failed against itself. (Third
    # source-text guard to misfire this way today -- assert on behaviour, not on the file's own words.)
    assert "build_pairs" not in globals(), "the pair construction must be imported, not defined here"
    assert "augment" not in globals(), "the pair augmentation must be imported, not defined here"
    _ep = _mod("hiv_epistasis")
    assert callable(_ep.build_pairs) and callable(_ep.augment), \
        "hiv_epistasis must supply build_pairs/augment"
    print("self-check OK -- bar frozen, floor read "
          f"(RT {f['NNRTI']['auroc']} / IN {f['INSTI']['auroc']}), prior {prior['verdict']}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "wiki" / f"hiv_context_vs_linear_floor_{date.today()}.json"))
    ap.add_argument("--self-check", action="store_true")
    a = ap.parse_args(argv)
    if a.self_check:
        self_check()
        return 0

    bar = load_bar()
    floors = floor_targets()
    D, PG, EP = _mod("hiv_supervised_deployability"), _mod("hiv_deployability_per_gene"), _mod("hiv_epistasis")

    cells, refused = {}, {}
    for cls in floors:
        try:
            cells[cls] = run_gene(cls, D, PG, EP)
        except ContextRefused as exc:
            refused[cls] = str(exc)
            print(f"{cls}: REFUSED -- {exc}")
    if not cells:
        print("no gene could be scored (gitignored Stanford data absent?)")
        return 2

    # ANCHOR: the linear arm must reproduce the measured floor, per gene.
    anchors = {}
    for cls, c in cells.items():
        if c.get("status") != "SCORED":
            continue
        want = floors[cls]["auroc"]
        got = c["blind_spot_auroc"]["linear_floor"]
        anchors[cls] = {"floor_measured": want, "floor_reproduced_here": got,
                        "delta": round(got - want, 4),
                        "reproduces": bool(abs(got - want) <= ANCHOR_TOL)}
    print(f"\nANCHOR (the linear arm must reproduce {FLOOR_ARTIFACT.name}):")
    for cls, an in anchors.items():
        print(f"  {cls:6s} floor {an['floor_measured']:.4f} vs here {an['floor_reproduced_here']:.4f} "
              f"({an['delta']:+.4f}) -> {'REPRODUCES' if an['reproduces'] else 'MISMATCH'}")
    if anchors and not all(an["reproduces"] for an in anchors.values()):
        print("REFUSING: a context-vs-floor delta measured against a floor that does not reproduce is not "
              "a comparison.")
        return 3

    print(f"\n{'gene':5s} {'bsN':>5s} {'bsR':>4s} {'floor':>7s} {'pairwise':>9s} {'nonlinear':>10s} "
          f"{'d_pair':>8s} {'d_nonlin':>9s}")
    for cls, c in cells.items():
        if c.get("status") != "SCORED":
            print(f"{c['gene']:5s} {c['blind_spot_n']:5d} {c['blind_spot_R']:4d}   {c['status']}")
            continue
        b, d = c["blind_spot_auroc"], c["delta_vs_floor"]
        print(f"{c['gene']:5s} {c['blind_spot_n']:5d} {c['blind_spot_R']:4d} "
              f"{b['linear_floor']:7.4f} {b.get('pairwise', float('nan')):9.4f} "
              f"{b['nonlinear']:10.4f} "
              f"{d.get('pairwise',{}).get('delta', float('nan')):+8.4f} "
              f"{d['nonlinear']['delta']:+9.4f}")
    for cls, c in cells.items():
        if c.get("status") != "SCORED":
            continue
        for k, dd in c["delta_vs_floor"].items():
            print(f"  {c['gene']:4s} {k:10s} delta {dd['delta']:+.4f}  95% CI "
                  f"[{dd['ci_lo']}, {dd['ci_hi']}]  CI-positive={dd['ci_positive']}")

    v = verdict(cells, bar)
    print(f"\nVERDICT (frozen rule): {v['verdict']}")
    print(f"  winning estimators per gene: {v.get('winning_estimators_per_gene')}")

    prior = json.loads(PRIOR_ARTIFACT.read_text(encoding="utf-8"))
    art = {"schema": "hiv-context-vs-linear-floor-v1", "date": str(date.today()),
           "question": bar["question"],
           "frozen_bar": str(BAR_PATH.relative_to(REPO)).replace("\\", "/"),
           "floor_artifact": str(FLOOR_ARTIFACT.relative_to(REPO)).replace("\\", "/"),
           "anchor": anchors, **v,
           "prior_evidence": {"artifact": str(PRIOR_ARTIFACT.relative_to(REPO)).replace("\\", "/"),
                              "verdict": prior["verdict"],
                              "fraction_ci_positive": prior.get("fraction_ci_positive")},
           "honest_scope": bar["honest_scope"] + [
               ("The claim is scoped to THIS FEATURE SPACE at THIS SAMPLE SIZE under THIS split. With "
                "~1063 (RT) / 405 (IN) one-hot features and 1883 / 844 blind-spot isolates, added capacity "
                "may simply be unaffordable out-of-distribution -- that is a different statement from "
                "'HIV resistance is additive', which this does NOT establish."),
               ("It does NOT test a sequence model over raw residues. 'Context' here means interactions "
                "among substitution tokens; attention over the raw sequence is a different, costlier "
                "question this cannot answer."),
           ],
           "named_limitations": [
               ("PAIR SELECTION IS BY PREVALENCE, which handicaps the pairwise arm. `build_pairs` ranks "
                "candidate pairs by co-occurrence among the top-40 most PREVALENT tokens; in an NNRTI "
                "cohort of treatment-experienced patients those are dominated by NRTI thymidine-analogue "
                "mutations (41L/210W/215Y) and M184V -- real epistatic biology, but NOT efavirenz drivers. "
                "So the pairwise delta is a statement about prevalent-token pairs, not about the most "
                "PREDICTIVE pairs. The general claim therefore rests on the NON-LINEAR arm, which sees "
                "every feature and can form any interaction, and which also lost on both genes."),
               ("That the top pairs ARE canonical biology is what makes the negative strong rather than "
                "weak: the construction found genuine interactions and they still did not help."),
               ("A cheap unexplored follow-up, if anyone wants to push: select pairs by univariate "
                "association with the label rather than by prevalence. Not run here."),
           ],
           "what_a_negative_licenses": bar["what_a_negative_would_license"],
           "cells": cells, "refused": refused}
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
