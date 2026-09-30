"""GLM headroom across ALL RT drugs — generalizing the EFV result off one drug (2026-09-30).

THE GAP THIS CLOSES. `scripts/glm_alphabet_headroom.py` showed that a SUPERVISED model over
position-resolved genotype tokens clears the pre-registered bar that ZERO-SHOT ESM2 failed, on the
catalog's own blind spot — but on ONE drug (EFV) and ONE gene (RT). "One drug, one gene" was the largest
honest limit on that result.

WHY THIS IS NEARLY FREE. NNRTI and NRTI drugs SHARE the RT protein, so the cached masked-marginals
(`data/processed/hiv_rt_esm650m_masked_marginals.json`) already cover every drug here, and there is a
COMMITTED zero-shot all-RT baseline to compare against: `wiki/hiv_esm_vs_catalog_allrt_2026-07-09.json`
(6 scoreable RT drugs, median subset AUROC 0.4558, **0 of 6 genuine passes**).

REUSE, NOT RE-DERIVATION — this is what lets the anchors bind. The catalog-negative subset is built with
the SAME `DRIFT` / `CUTOFF` / `subs_of` recipe and the SAME deployed comparator
(`hiv_amr.call_hiv_observed`) as the committed run, imported from those modules rather than re-written.
The committed run's own `UNDERPOWERED` guard (<10 per class in the subset) is reused verbatim, which is
why 6 of 11 RT drugs are scoreable and the other 5 are reported as underpowered rather than scored.

PER-DRUG PIPELINE ANCHOR. For each drug we recompute the zero-shot subset AUROC and require it to match
the committed value. A drug whose anchor misses is reported `ANCHOR_MISMATCH` and its supervised number is
WITHHELD — a head-to-head on a differently-built subset is not a head-to-head. Six anchors, not one.

TWO DATASETS, because grouping needs metadata the anchored files do not carry:
  * non-Full (`*_DataSet.txt`) is what the committed baseline used -> the ANCHOR, and an ungrouped
    supervised number that is honestly labelled OPTIMISTIC.
  * `.Full` carries `RefID`/`PtID` -> LEAVE-ONE-STUDY-OUT, the honest supervised headline.

THE BAR IS IMPORTED, NOT RESTATED — `PREREGISTERED` / `verdict` come from `glm_alphabet_headroom.py`, so
the bar has ONE definition across both scripts and cannot drift between them.

Run: uv run python scripts/glm_headroom_allrt.py [--self-check]
Exits 2 when the gitignored Stanford data or the ESM cache is absent.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import statistics as st
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / "scripts"))

COMMITTED = REPO / "wiki" / "hiv_esm_vs_catalog_allrt_2026-07-09.json"
ANCHOR_TOL = 0.02
MIN_PER_CLASS = 10          # the committed run's own UNDERPOWERED guard, reused verbatim
SEED = 0

# NOMINAL-vs-GENUINE, the trap this repo already documented. The committed zero-shot run reports
# `n_drugs_passing_bar_nominal: 1` beside `n_drugs_passing_bar_genuine: 0` -- the nominal pass is DOR,
# which clears >=0.65 on a 37-isolate subset with mutation-burden at 0.783 and a near-chance catalog.
# A raw pass-count therefore launders a small-subset artefact into a result, so the SAME split is applied
# to the supervised arm and the excluded drugs are NAMED with their n rather than quietly dropped.
SMALL_SUBSET_N = 100


def summarize(scored: list[dict]) -> dict:
    """Nominal AND genuine pass counts. Never report one without the other."""
    big = [c for c in scored if c["subset_n"] >= SMALL_SUBSET_N]
    small = [c for c in scored if c["subset_n"] < SMALL_SUBSET_N]
    return {
        "n_scored": len(scored), "n_genuine": len(big),
        "median_zero_shot": st.median([c["zero_shot_esm2_auroc"] for c in scored]),
        "median_supervised": st.median([c["supervised_auroc"] for c in scored]),
        "supervised_pass_nominal": sum(1 for c in scored if c["verdict_supervised"] == "PASS"),
        "zero_shot_pass_nominal": sum(1 for c in scored if c["verdict_zero_shot"] == "PASS"),
        "supervised_pass_genuine": sum(1 for c in big if c["verdict_supervised"] == "PASS"),
        "zero_shot_pass_genuine": sum(1 for c in big if c["verdict_zero_shot"] == "PASS"),
        "small_subsets": ", ".join(f"{c['drug']}(n={c['subset_n']})" for c in small),
        "median_supervised_genuine": (st.median([c["supervised_auroc"] for c in big]) if big else None),
        "median_zero_shot_genuine": (st.median([c["zero_shot_esm2_auroc"] for c in big]) if big else None),
    }


def _bar():
    spec = importlib.util.spec_from_file_location(
        "glm_alphabet_headroom", REPO / "scripts" / "glm_alphabet_headroom.py")
    G = importlib.util.module_from_spec(spec); spec.loader.exec_module(G)
    return G


class AllRtRefused(RuntimeError):
    """Cannot be trusted -- raised instead of returning a number."""


def _auroc(y, s) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))


def _oof(X, y, groups=None):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
    clf = LogisticRegression(max_iter=3000, C=1.0, solver="liblinear")
    if groups is not None and len(set(groups)) >= 2:
        return cross_val_predict(clf, X, y, groups=groups,
                                 cv=GroupKFold(n_splits=min(5, len(set(groups)))),
                                 method="predict_proba")[:, 1]
    return cross_val_predict(clf, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                             method="predict_proba")[:, 1]


def score_drug(col: str, drug: str, cls: str, path: Path, committed: dict | None) -> dict:
    """One drug. Anchors zero-shot against the committed run, then scores the supervised arm."""
    import numpy as np

    import dna_decode.data.hiv_amr as H
    from scripts.hiv_esm_vs_catalog import LOGP_CACHE, isolate_muts, rt_protein
    from scripts.hiv_esm_vs_catalog_allrt import CUTOFF, DRIFT
    spec = importlib.util.spec_from_file_location(
        "hiv_supervised_vs_catalog", REPO / "scripts" / "hiv_supervised_vs_catalog.py")
    S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)

    prot = rt_protein()
    cache = Path(LOGP_CACHE)
    if not cache.is_absolute():
        cache = REPO / cache
    if not cache.exists():
        raise AllRtRefused(f"masked-marginal cache absent: {cache}")
    logp = {int(k): v for k, v in json.loads(cache.read_text(encoding="utf-8")).items()}
    if not path.exists():
        raise AllRtRefused(f"dataset absent: {path}")

    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    if col not in rows[0]:
        raise AllRtRefused(f"fold column {col} absent in {path.name}")

    def subs_of(r):                                   # the committed recipe, verbatim
        return {f"{prot[p-1]}{p}{aa}" for p, aa in isolate_muts(r, pcols)
                if p <= len(prot) and p not in DRIFT and prot[p - 1] != aa}

    def esm_score(r):
        d = [logp[p][prot[p - 1]] - logp[p][aa] for p, aa in isolate_muts(r, pcols)
             if p in logp and p <= len(prot) and prot[p - 1] != aa and aa in logp[p]
             and prot[p - 1] in logp[p]]
        return float(np.mean(d)) if d else 0.0

    have = [r for r in rows if r[col] not in ("NA", "", "-", None)]
    cat = [H.call_hiv_observed(drug, {"RT": subs_of(r)}).prediction == "R" for r in have]
    sub = [r for r, c in zip(have, cat) if not c]
    ys = np.array([1 if float(r[col]) >= CUTOFF else 0 for r in sub])

    out = {"drug": col, "drug_full": drug, "cls": cls, "dataset": path.name,
           "n_with_fold": len(have), "subset_n": len(sub), "subset_r": int(ys.sum())}
    if ys.sum() < MIN_PER_CLASS or (len(ys) - ys.sum()) < MIN_PER_CLASS:
        out["status"] = "UNDERPOWERED"      # the committed run's own guard, so the drug sets agree
        return out

    out["zero_shot_esm2_auroc"] = round(_auroc(ys, np.array([esm_score(r) for r in sub])), 4)
    out["burden_auroc"] = round(_auroc(ys, np.array(
        [float(len(isolate_muts(r, pcols))) for r in sub])), 4)

    # --- per-drug anchor against the committed zero-shot run ---
    if committed is not None:
        dn, dr = committed["subset_n"], committed["subset_r"]
        de = committed["esm_subset_auroc"]
        ok = (dn == len(sub) and dr == int(ys.sum())
              and abs(out["zero_shot_esm2_auroc"] - de) <= ANCHOR_TOL)
        out["anchor"] = {"committed_subset_n": dn, "committed_subset_r": dr,
                         "committed_esm_subset_auroc": de, "reproduces": bool(ok)}
        if not ok:
            out["status"] = "ANCHOR_MISMATCH"
            return out                      # supervised number WITHHELD on a mismatched subset

    X, feats = S.build_onehot(sub, pcols, set(DRIFT), prot)
    out["n_features"] = len(feats)
    groups = np.array([r["RefID"] for r in sub]) if "RefID" in rows[0] else None
    out["split"] = "leave_one_study_out" if groups is not None else \
        "ungrouped_stratified_5fold_OPTIMISTIC"
    out["supervised_auroc"] = round(_auroc(ys, _oof(X, ys, groups)), 4)

    rng = np.random.default_rng(SEED)
    nulls = []
    for _ in range(3):
        yp = rng.permutation(ys)
        if 0 < yp.sum() < len(yp):
            nulls.append(_auroc(yp, _oof(X, yp, groups)))
    out["shuffled_null_auroc"] = round(float(st.mean(nulls)), 4) if nulls else None

    G = _bar()
    out["verdict_supervised"] = G.verdict(out["supervised_auroc"], out["burden_auroc"],
                                          out["shuffled_null_auroc"] or 1.0)
    out["verdict_zero_shot"] = G.verdict(out["zero_shot_esm2_auroc"], out["burden_auroc"],
                                         out["shuffled_null_auroc"] or 1.0)
    out["delta"] = round(out["supervised_auroc"] - out["zero_shot_esm2_auroc"], 4)
    out["status"] = "SCORED"
    return out


def self_check() -> None:
    G = _bar()
    assert G.PREREGISTERED["min_auroc"] == 0.65, "the bar must be the imported one"
    assert G.verdict(0.81, 0.45, 0.50) == "PASS"
    z = G.PREREGISTERED["zero_shot_result_being_compared"]["esm_mean_damage_auroc"]
    assert G.verdict(z, 0.4514, 0.4974) != "PASS", "the compared zero-shot number must still fail"
    assert MIN_PER_CLASS == 10, "the underpowered guard must match the committed run's"
    print("self-check OK")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "wiki" / f"glm_headroom_allrt_{date.today()}.json"))
    ap.add_argument("--self-check", action="store_true")
    a = ap.parse_args(argv)
    if a.self_check:
        self_check()
        return 0

    from scripts.hiv_esm_vs_catalog_allrt import ARMS
    if not COMMITTED.exists():
        print(f"committed zero-shot baseline absent: {COMMITTED}")
        return 2
    base = json.loads(COMMITTED.read_text(encoding="utf-8"))
    by_drug = {r["drug"]: r for r in base["per_drug"]}

    cells, refused = [], []
    for cls, rel, drugs in ARMS:
        for variant in ("anchor", "full"):
            p = REPO / (rel if variant == "anchor" else rel.replace(".txt", ".Full.txt"))
            for col, drug in drugs:
                try:
                    c = score_drug(col, drug, cls, p,
                                   by_drug.get(col) if variant == "anchor" else None)
                    c["variant"] = variant
                    cells.append(c)
                except AllRtRefused as exc:
                    refused.append({"drug": col, "variant": variant, "refused": str(exc)})
    if not cells:
        print("no trustworthy cell (data or ESM cache absent)")
        return 2

    anchors = [c for c in cells if c["variant"] == "anchor" and "anchor" in c]
    nbind = sum(1 for c in anchors if c["anchor"]["reproduces"])
    print(f"\nPER-DRUG ANCHORS: {nbind}/{len(anchors)} reproduce the committed zero-shot run "
          f"(tol {ANCHOR_TOL})")
    if anchors and nbind == 0:
        print("  REFUSING: no anchor binds -- subset construction disagrees with the committed run.")
        return 3

    for variant, title in (("anchor", "NON-FULL (anchored to the committed run; ungrouped = OPTIMISTIC)"),
                           ("full", ".FULL (leave-one-STUDY-out -- the honest supervised headline)")):
        rows = [c for c in cells if c["variant"] == variant]
        print(f"\n=== {title} ===")
        print(f"{'drug':5s} {'cls':6s} {'subN':>5s} {'subR':>5s} {'zero':>7s} {'super':>7s} "
              f"{'burden':>7s} {'null':>6s} {'delta':>7s}  status")
        for c in sorted(rows, key=lambda r: (r["cls"], r["drug"])):
            if c.get("status") != "SCORED":
                print(f"{c['drug']:5s} {c['cls']:6s} {c['subset_n']:5d} {c['subset_r']:5d} "
                      f"{'':>7s} {'':>7s} {'':>7s} {'':>6s} {'':>7s}  {c.get('status')}")
                continue
            print(f"{c['drug']:5s} {c['cls']:6s} {c['subset_n']:5d} {c['subset_r']:5d} "
                  f"{c['zero_shot_esm2_auroc']:7.4f} {c['supervised_auroc']:7.4f} "
                  f"{c['burden_auroc']:7.4f} {c['shuffled_null_auroc'] or 0:6.3f} "
                  f"{c['delta']:+7.4f}  {c['verdict_supervised']}")
        sc = [c for c in rows if c.get("status") == "SCORED"]
        if sc:
            summ = summarize(sc)
            print(f"  median zero-shot {summ['median_zero_shot']:.4f}  "
                  f"median supervised {summ['median_supervised']:.4f}")
            print(f"  NOMINAL  supervised {summ['supervised_pass_nominal']}/{summ['n_scored']}"
                  f"   zero-shot {summ['zero_shot_pass_nominal']}/{summ['n_scored']}")
            print(f"  GENUINE (subset_n >= {SMALL_SUBSET_N})  supervised "
                  f"{summ['supervised_pass_genuine']}/{summ['n_genuine']}"
                  f"   zero-shot {summ['zero_shot_pass_genuine']}/{summ['n_genuine']}"
                  f"   [small, excluded: {summ['small_subsets'] or 'none'}]")

    G = _bar()
    art = {"schema": "glm-headroom-allrt-v1", "date": str(date.today()),
           "question": ("does the supervised genotype-token result generalize off ONE drug, judged by the "
                        "SAME pre-registered bar, against the COMMITTED zero-shot all-RT baseline?"),
           "preregistered_bar": G.PREREGISTERED,
           "committed_zero_shot_baseline": {
               "artifact": str(COMMITTED.relative_to(REPO)).replace("\\", "/"),
               "median_esm_subset_auroc": base.get("median_esm_subset_auroc"),
               "n_drugs_passing_bar_genuine": base.get("n_drugs_passing_bar_genuine")},
           "anchors_reproducing": f"{nbind}/{len(anchors)}",
           "summary_by_variant": {v: summarize([c for c in cells
                                                if c["variant"] == v and c.get("status") == "SCORED"])
                                  for v in ("anchor", "full")
                                  if any(c["variant"] == v and c.get("status") == "SCORED"
                                         for c in cells)},
           "nominal_vs_genuine": (
               "A raw pass-count is a TRAP this repo already documented: the committed zero-shot run "
               f"reports nominal 1 / genuine 0, the nominal pass being DOR on a {37}-isolate subset with "
               "burden 0.783. The same nominal/genuine split is applied to the supervised arm and every "
               f"drug with subset_n < {SMALL_SUBSET_N} is NAMED, not silently dropped."),
           "honest_scope": ("Linear model over one-hot substitution tokens = the WEAKEST supervised "
                            "family member, so a pass is a FLOOR on that family and says nothing about "
                            "whether attention/context helps. IN-DISTRIBUTION to Stanford. All drugs "
                            "share ONE gene (RT), so this generalizes across DRUGS, NOT across genes. "
                            "The 5 underpowered drugs are reported, never scored."),
           "cells": cells, "refused": refused}
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
