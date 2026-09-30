"""GLM headroom: does a SUPERVISED model over genotype tokens beat a ZERO-SHOT protein LM on the
curated catalog's own blind spot? (2026-09-30)

WHY THIS EXISTS. Two results in this repo look contradictory and the difference decides whether a
genotype language model has anything to learn:

  * 2026-07-09 CLOSED NEGATIVE -- "a learned variant scorer does NOT fill the curated AMR catalog's blind
    spot". ESM2-650M, ZERO-SHOT, on the catalog-negative subset: AUROC 0.4485 (n=1111, 53R) against a
    PRE-REGISTERED bar of `>=0.65 AND > mutation-burden AND null < 0.55`. FAILED.
  * 2026-09-30 -- the SHIPPED supervised complement flags that same blind spot at 15-16x enrichment on
    two independent drug classes.

Same blind spot. Same label. The methods differ in exactly one way: **zero-shot likelihood vs supervised
fitting over position-resolved substitution tokens.** That is the GLM's central question, so this script
puts both on the IDENTICAL isolate set and judges them by the SAME pre-registered bar.

WHAT IT DOES NOT DO. It does not test a transformer, and it does not test whether attention/context helps
-- it tests whether the SUPERVISED DIRECTION has headroom the zero-shot direction does not. A linear model
over one-hot substitution tokens is the WEAKEST possible member of that family, so a pass is a floor on
the family, not a ceiling (and a fail would be strong evidence against it).

THE BAR IS NOT MINE AND IS NOT MOVED. `PREREGISTERED` below is copied verbatim from the committed
`wiki/hiv_esm_vs_catalog_2026-07-09.json` `pass_bar` field -- the bar the zero-shot arm was held to and
failed. Judging a new method by a fresh bar is how a negative gets laundered into a positive.

HONEST DESIGN NOTES, each load-bearing:
  * `.Full` is PRIMARY because it alone carries `RefID`/`PtID`, so the supervised model can be scored
    out-of-fold by LEAVE-ONE-STUDY-OUT (GroupKFold by RefID) -- predicting studies it never trained on.
    Plain k-fold would leak: the metadata shows 1.35x patient duplication.
  * The non-Full dataset is the PIPELINE ANCHOR: ESM2 there must reproduce the committed 0.4485. If it
    does not, this script's subset construction disagrees with the 2026-07-09 run and NOTHING else it
    prints can be trusted -- so that check gates the whole comparison.
  * ESM2 is zero-shot by construction (nothing is fit), so its AUROC needs no held-out split. The
    supervised side gets the strictest split available. That asymmetry FAVOURS ESM2, deliberately.
  * Both arms see ONLY the isolate's own substitutions. Neither sees the catalog.

Run: uv run python scripts/glm_alphabet_headroom.py [--self-check]
Needs the gitignored Stanford data + the cached masked-marginals (`data/processed/hiv_rt_esm650m_masked_
marginals.json`, ~30 min CPU to regenerate); exits 2 when absent rather than reporting an empty result.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / "scripts"))

RAW = REPO / "data" / "raw" / "hiv"
LOGP_CACHE = REPO / "data" / "processed" / "hiv_rt_esm650m_masked_marginals.json"
CUTOFF = 3.0
DRUG = "EFV"
SEED = 0

# VERBATIM from wiki/hiv_esm_vs_catalog_2026-07-09.json -- the bar the ZERO-SHOT arm failed.
PREREGISTERED = {
    "source": "wiki/hiv_esm_vs_catalog_2026-07-09.json :: pass_bar",
    # VERBATIM as committed, kept exact so the citation is traceable (a test asserts byte-equality with
    # that artifact). My own test caught a PARAPHRASE here -- the thresholds were identical but 'esm' had
    # become 'score', and a paraphrased bar cannot be verified against its source.
    "rule_verbatim": "esm>=0.65 AND esm>burden AND null<0.55",
    # The ONLY generalisation: the subject is the method under test, not ESM specifically. Every
    # threshold is unchanged. Stated explicitly rather than folded silently into the quote.
    "rule_as_applied": "score>=0.65 AND score>burden AND null<0.55",
    "generalisation": "subject widened from 'esm' to any scoring method; NO threshold altered",
    "min_auroc": 0.65,
    "must_beat_burden": True,
    "max_null_auroc": 0.55,
    "zero_shot_result_being_compared": {"esm_mean_damage_auroc": 0.4485, "n": 1111, "n_resistant": 53,
                                        "passed": False},
}
ANCHOR_TOL = 0.02          # the non-Full ESM AUROC must land within this of the committed 0.4485


class HeadroomRefused(RuntimeError):
    """The comparison cannot be trusted -- raised instead of returning a number."""


def verdict(score: float, burden: float, null: float) -> str:
    """The pre-registered bar, applied mechanically. No clause is re-weighted after seeing a result."""
    if score < PREREGISTERED["min_auroc"]:
        return "FAIL_BELOW_MIN_AUROC"
    if PREREGISTERED["must_beat_burden"] and not score > burden:
        return "FAIL_NOT_BETTER_THAN_MUTATION_BURDEN"
    if null >= PREREGISTERED["max_null_auroc"]:
        return "FAIL_NULL_TOO_HIGH"
    return "PASS"


def _mods():
    spec = importlib.util.spec_from_file_location(
        "hiv_supervised_vs_catalog", REPO / "scripts" / "hiv_supervised_vs_catalog.py")
    S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
    return S, S.ev


def _load_logp() -> dict[int, dict[str, float]]:
    """Cached ESM2 masked marginals. JSON round-trips position keys as STRINGS -- coercing them is the
    documented trap: an uncoerced table silently scores nothing and prints a clean-looking zero."""
    if not LOGP_CACHE.exists():
        raise HeadroomRefused(f"masked-marginal cache absent: {LOGP_CACHE}")
    raw = json.loads(LOGP_CACHE.read_text(encoding="utf-8"))
    out = {int(k): v for k, v in raw.items()}
    if not out:
        raise HeadroomRefused("masked-marginal cache is empty")
    return out


def build(dataset: str):
    """Catalog-negative subset + its genotype-token matrix, labels, groups, ESM + burden scores."""
    import numpy as np

    import dna_decode.data.hiv_amr as H
    S, ev = _mods()
    path = RAW / ("NNRTI_DataSet" + (".Full.txt" if dataset == "full" else ".txt"))
    if not path.exists():
        raise HeadroomRefused(f"dataset absent: {path}")

    prot = ev.rt_protein()
    for pos, wt in H._RT_WT.items():                      # integrity gate, same as the 07-09 run
        if prot[pos - 1] != wt:
            raise HeadroomRefused(f"reference/catalog WT mismatch at {pos}")

    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = {int(c[1:]) for c in pcols
               if int(c[1:]) <= len(prot)
               and sum(1 for r in rows if (r[c] or "").strip() == prot[int(c[1:]) - 1]) > 0.01 * len(rows)}
    if len(drifted) >= 20:
        raise HeadroomRefused(f"drift detector over-fired ({len(drifted)} positions)")

    have = [r for r in rows if r.get(DRUG) not in ("NA", "", "-", None)]
    majors = H.NNRTI_RT_MAJOR_DRMS

    def catalog_positive(r):
        return any(f"{H._RT_WT.get(p, '?')}{p}{aa}" in majors for p, aa in ev.isolate_muts(r, pcols))

    sub = [r for r in have if not catalog_positive(r)]
    if not sub:
        raise HeadroomRefused("catalog-negative subset is empty -- subset construction is broken")
    y = np.array([1 if float(r[DRUG]) >= CUTOFF else 0 for r in sub])
    if y.sum() == 0 or y.sum() == len(y):
        raise HeadroomRefused("the blind spot has only one class -- no AUROC is defined")

    X, feats = S.build_onehot(sub, pcols, drifted, prot)
    logp = _load_logp()

    def esm(r):
        d = [logp[p][prot[p - 1]] - logp[p][aa] for p, aa in ev.isolate_muts(r, pcols)
             if p in logp and p <= len(prot) and prot[p - 1] != aa and aa in logp[p]
             and prot[p - 1] in logp[p]]
        return float(np.mean(d)) if d else 0.0

    esm_s = np.array([esm(r) for r in sub])
    burden = np.array([float(len([1 for p, aa in ev.isolate_muts(r, pcols)
                                  if p <= len(prot) and prot[p - 1] != aa])) for r in sub])
    groups = {}
    for key in ("RefID", "PtID"):
        if key in rows[0]:
            groups[key] = np.array([r[key] for r in sub])
    return dict(sub=sub, X=X, n_features=len(feats),
                feature_names=[f"{p}{aa}" for (p, aa) in feats],
                y=y, esm=esm_s, burden=burden, groups=groups,
                n=len(sub), dataset=dataset, n_drifted=len(drifted))


def _auroc(y, s) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))


def _oof(X, y, groups=None):
    """Out-of-fold predictions. GROUPED when a group key is available -- plain k-fold leaks (1.35x
    patient duplication), so an ungrouped number would be optimistic and is only used when nothing else
    exists (and is labelled `ungrouped` in the output)."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
    clf = LogisticRegression(max_iter=3000, C=1.0, solver="liblinear")
    if groups is not None:
        n = min(5, len(set(groups)))
        if n < 2:
            raise HeadroomRefused("fewer than 2 groups -- grouped CV is undefined")
        return cross_val_predict(clf, X, y, groups=groups, cv=GroupKFold(n_splits=n),
                                 method="predict_proba")[:, 1]
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    return cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]


def run(dataset: str) -> dict:
    import numpy as np
    d = build(dataset)
    y, X = d["y"], d["X"]
    rng = np.random.default_rng(SEED)

    out = {"dataset": dataset, "n_blind_spot": d["n"], "n_resistant": int(y.sum()),
           "n_susceptible": int(len(y) - y.sum()), "n_features": d["n_features"],
           "base_rate": round(float(y.mean()), 4), "cutoff_fold": CUTOFF,
           "zero_shot_esm2_650m_auroc": round(_auroc(y, d["esm"]), 4),
           "mutation_burden_auroc": round(_auroc(y, d["burden"]), 4)}

    # label-shuffled null on the SUPERVISED pipeline (the arm under test), not on a bare score
    nulls = []
    for _ in range(5):
        yp = rng.permutation(y)
        if 0 < yp.sum() < len(yp):
            nulls.append(_auroc(yp, _oof(X, yp, d["groups"].get("RefID"))))
    out["shuffled_null_auroc"] = round(float(np.mean(nulls)), 4) if nulls else None

    sup = {}
    for label, gkey in (("leave_one_study_out", "RefID"), ("patient_grouped", "PtID")):
        if gkey in d["groups"]:
            sup[label] = round(_auroc(y, _oof(X, y, d["groups"][gkey])), 4)
    if not sup:
        sup["ungrouped_stratified_5fold_OPTIMISTIC"] = round(_auroc(y, _oof(X, y)), 4)
    out["supervised_genotype_tokens_auroc"] = sup

    strictest = sup.get("leave_one_study_out") or sup.get("patient_grouped") or next(iter(sup.values()))
    out["strictest_supervised_auroc"] = strictest
    out["strictest_split"] = ("leave_one_study_out" if "leave_one_study_out" in sup else
                              ("patient_grouped" if "patient_grouped" in sup else
                               "ungrouped_stratified_5fold_OPTIMISTIC"))
    out["verdict_supervised"] = verdict(strictest, out["mutation_burden_auroc"],
                                        out["shuffled_null_auroc"] or 1.0)
    out["verdict_zero_shot"] = verdict(out["zero_shot_esm2_650m_auroc"], out["mutation_burden_auroc"],
                                       out["shuffled_null_auroc"] or 1.0)
    out["delta_supervised_minus_zero_shot"] = round(strictest - out["zero_shot_esm2_650m_auroc"], 4)

    # WHAT is it using? The subset is catalog-NEGATIVE by construction, so no catalogued major DRM is
    # present in ANY isolate -- the model therefore cannot be rediscovering the catalog, and its top
    # features name the non-catalogued substitutions carrying the blind spot's signal.
    from sklearn.linear_model import LogisticRegression
    full = LogisticRegression(max_iter=3000, C=1.0, solver="liblinear").fit(X, y)
    order = np.argsort(full.coef_[0])[::-1][:12]
    out["top_resistance_features_full_fit"] = [
        {"token": d["feature_names"][i], "coef": round(float(full.coef_[0][i]), 3),
         "n_carriers": int(X[:, i].sum())} for i in order]
    return out


def self_check() -> None:
    assert verdict(0.81, 0.45, 0.50) == "PASS"
    assert verdict(0.60, 0.45, 0.50) == "FAIL_BELOW_MIN_AUROC"
    assert verdict(0.70, 0.75, 0.50) == "FAIL_NOT_BETTER_THAN_MUTATION_BURDEN"
    assert verdict(0.70, 0.45, 0.60) == "FAIL_NULL_TOO_HIGH"
    # the committed zero-shot result must FAIL the bar it was held to -- if it passes, the bar drifted
    z = PREREGISTERED["zero_shot_result_being_compared"]["esm_mean_damage_auroc"]
    assert verdict(z, 0.4514, 0.4974) != "PASS", "the 2026-07-09 bar no longer fails its own result"
    print("self-check OK")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "wiki" / f"glm_alphabet_headroom_{date.today()}.json"))
    ap.add_argument("--self-check", action="store_true")
    a = ap.parse_args(argv)
    if a.self_check:
        self_check()
        return 0

    cells, refused = [], []
    for ds in ("validation", "full"):
        try:
            cells.append(run(ds))
        except HeadroomRefused as exc:
            refused.append({"dataset": ds, "refused": str(exc)})
            print(f"REFUSED {ds}: {exc}")
    if not cells:
        print("no trustworthy cell (gitignored data or ESM cache absent)")
        return 2

    # PIPELINE ANCHOR: on the non-Full set our ESM must reproduce the committed 0.4485.
    anchor = next((c for c in cells if c["dataset"] == "validation"), None)
    anchor_note = "not evaluated (validation dataset refused)"
    if anchor:
        got = anchor["zero_shot_esm2_650m_auroc"]
        exp = PREREGISTERED["zero_shot_result_being_compared"]["esm_mean_damage_auroc"]
        ok = abs(got - exp) <= ANCHOR_TOL
        anchor_note = (f"ESM2 on the non-Full catalog-negative subset = {got} vs committed {exp} "
                       f"({'REPRODUCES' if ok else 'DISAGREES'}, tol {ANCHOR_TOL})")
        print(f"\nPIPELINE ANCHOR: {anchor_note}")
        if not ok:
            print("  REFUSING to publish a comparison whose subset disagrees with the committed run.")
            return 3

    print()
    for c in cells:
        print(f"[{c['dataset']}] blind spot n={c['n_blind_spot']}  R={c['n_resistant']}  "
              f"feats={c['n_features']}")
        print(f"    zero-shot ESM2-650M       {c['zero_shot_esm2_650m_auroc']:.4f}  "
              f"-> {c['verdict_zero_shot']}")
        for k, v in c["supervised_genotype_tokens_auroc"].items():
            print(f"    supervised ({k}) {v:.4f}")
        print(f"    mutation burden           {c['mutation_burden_auroc']:.4f}")
        print(f"    shuffled null             {c['shuffled_null_auroc']}")
        print(f"    STRICTEST supervised      {c['strictest_supervised_auroc']:.4f} "
              f"({c['strictest_split']}) -> {c['verdict_supervised']}   "
              f"delta vs zero-shot {c['delta_supervised_minus_zero_shot']:+.4f}\n")

    art = {"schema": "glm-alphabet-headroom-v1", "date": str(date.today()),
           "question": ("does a SUPERVISED model over position-resolved genotype tokens clear the SAME "
                        "pre-registered bar that ZERO-SHOT ESM2-650M failed, on the SAME catalog blind "
                        "spot?"),
           "preregistered_bar": PREREGISTERED, "pipeline_anchor": anchor_note,
           "honest_scope": ("A linear model over one-hot substitution tokens is the WEAKEST member of the "
                            "supervised family, so a pass is a FLOOR on that family, not a ceiling, and "
                            "says nothing about whether context/attention helps. In-distribution to the "
                            "Stanford knowledge base. One drug (EFV), one gene (RT), one label cutoff."),
           "cells": cells, "refused": refused}
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
