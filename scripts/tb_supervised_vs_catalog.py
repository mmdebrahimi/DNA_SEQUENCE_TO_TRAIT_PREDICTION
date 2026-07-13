"""Does the supervised blind-spot rescue TRANSFER to bacteria? TB rifampicin (2026-07-12).

Recommendation 3, with the user's caveat: HIV and bacterial resistance behave differently for sequence
models, so the HIV win doesn't automatically transfer. This runs the SAME experiment on TB RIF (CRyPTIC):
a supervised model on the FULL set of rpoB variants (not just the WHO catalogue's grade-1/2 determinants)
vs the WHO catalogue, on the free measured CRyPTIC binary phenotype, with the bacterial de-confound =
LEAVE-ONE-LINEAGE-OUT (the documented bacterial trap is learning lineage, not mechanism).

Per split (plain 5-fold + leave-one-lineage-out): overall AUROC + BLIND-SPOT AUROC (isolates the WHO
catalogue calls susceptible = no grade-1/2 determinant). PASS bar (same as HIV): blind-spot AUROC >= 0.65
AND > mutation-burden AND null < 0.55, under LEAVE-LINEAGE-OUT.

VERDICT:
  TRANSFERS_TO_BACTERIA         — leave-lineage-out blind-spot PASSES (rescue is real beyond lineage).
  IN_DISTRIBUTION_ONLY          — plain 5-fold passes but leave-lineage-out collapses (learned lineage).
  DOES_NOT_TRANSFER             — no blind-spot signal even in-distribution (the user's caveat holds).

Streams VARIANTS.parquet once for rpoB + determinant + barcode positions (cached). Frozen decoder surface
untouched; reuses the TB scoring infra (tb_amr / tb_who_catalogue / tb_lineage) READ-only.
"""
from __future__ import annotations

import importlib.util
import json
import random
import statistics as st
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))

import score_tb_cryptic_parquet as P  # noqa: E402
from dna_decode.data import tb_who_catalogue, tb_lineage_barcode  # noqa: E402
from dna_decode.organism_rules import tb_amr, tb_lineage  # noqa: E402

DRUG, CODE = "rifampicin", "RIF"
RPOB = range(759807, 763326)          # H37Rv NC_000962.3 rpoB gene span
CACHE = REPO / "data" / "processed" / "tb_rif_rpob_cache.json"
MIN_CARRIERS, SEED = 5, 0
DUMP = P.DEFAULT_DUMP
REUSE = P.DEFAULT_REUSE


def build_cache(force=False):
    if CACHE.exists() and not force:
        return json.loads(CACHE.read_text())
    tb_who_catalogue.verify_pins()
    dets = tb_who_catalogue.load_determinants(DRUG)
    barcode = tb_lineage_barcode.load_barcode()
    det_pos = {d.pos for d in dets}
    bc_pos = P.barcode_positions(barcode)
    wanted = det_pos | bc_pos | set(RPOB)
    labels = P.load_labels(REUSE, CODE)
    print(f"[tb-sup] {len(labels)} HIGH {CODE} labels; streaming for {len(wanted)} positions "
          f"({len(det_pos)} det + {len(bc_pos)} barcode + rpoB)...", flush=True)
    calls, _ = P.load_calls_by_strain(DUMP / "VARIANTS.parquet", wanted)
    rows = {}
    for uid, lab in labels.items():
        c = calls.get(uid, {})
        rpob_vars = sorted(f"{pos}{vc.alt}" for pos, vc in c.items() if pos in RPOB)
        cat = tb_amr.score_drug(DRUG, c, dets).prediction
        lin = tb_lineage.assign_lineage(c, barcode)
        rows[uid] = {"label": lab, "rpob": rpob_vars, "catalog": cat, "lineage": lin, "n_var": len(c)}
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(rows))
    print(f"[tb-sup] cached {len(rows)} isolates -> {CACHE}", flush=True)
    return rows


def _auroc(y, s):
    y = np.asarray(y); s = np.asarray(s, float)
    pos, neg = y == 1, y == 0
    if pos.sum() == 0 or neg.sum() == 0:
        return None
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    return float((ranks[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum()))


def _blind(y, sc, neg, burden, seed=SEED):
    ys = [y[i] for i in range(len(y)) if neg[i]]
    ss = [sc[i] for i in range(len(y)) if neg[i]]
    bs = [burden[i] for i in range(len(y)) if neg[i]]
    if len(set(ys)) < 2:
        return {"n": len(ys), "R": int(sum(ys)), "auroc": None, "pass": False, "note": "single-class blind spot"}
    a, ab = _auroc(ys, ss), _auroc(ys, bs)
    rng = random.Random(seed); nulls = []
    for _ in range(200):
        yy = ys[:]; rng.shuffle(yy); nulls.append(_auroc(yy, ss))
    an = st.median([x for x in nulls if x is not None])
    return {"n": len(ys), "R": int(sum(ys)), "auroc": round(a, 4), "burden": round(ab, 4),
            "null": round(an, 4), "pass": bool(a >= 0.65 and a > ab and an < 0.55)}


def main(force=False):
    rows = build_cache(force)
    uids = [u for u, r in rows.items() if r["label"] in ("R", "S")]
    y = np.array([1 if rows[u]["label"] == "R" else 0 for u in uids])
    cat = np.array([1 if rows[u]["catalog"] == "R" else 0 for u in uids])
    neg = (cat == 0)
    lineage = np.array([rows[u]["lineage"] or "UNASSIGNED" for u in uids])
    burden = [rows[u]["n_var"] for u in uids]
    # feature vocab: rpoB variants in >= MIN_CARRIERS isolates
    from collections import Counter
    cnt = Counter(v for u in uids for v in rows[u]["rpob"])
    feats = sorted(v for v, n in cnt.items() if n >= MIN_CARRIERS)
    fidx = {v: j for j, v in enumerate(feats)}
    X = np.zeros((len(uids), len(feats)), np.float32)
    for i, u in enumerate(uids):
        for v in rows[u]["rpob"]:
            if v in fidx:
                X[i, fidx[v]] = 1.0
    clf = LogisticRegression(max_iter=2000, C=1.0, solver="liblinear")
    out = {}
    # plain 5-fold
    oof5 = cross_val_predict(clf, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                             method="predict_proba")[:, 1]
    out["plain_5fold"] = {"overall_auroc": round(_auroc(y.tolist(), oof5.tolist()), 4),
                          "blind_spot": _blind(y.tolist(), oof5.tolist(), neg.tolist(), burden)}
    # leave-one-lineage-out
    nsp = min(5, len(set(lineage.tolist())))
    oofl = cross_val_predict(clf, X, y, groups=lineage, cv=GroupKFold(nsp), method="predict_proba")[:, 1]
    out["leave_lineage_out"] = {"overall_auroc": round(_auroc(y.tolist(), oofl.tolist()), 4),
                                "blind_spot": _blind(y.tolist(), oofl.tolist(), neg.tolist(), burden)}
    bl = out["leave_lineage_out"]["blind_spot"]
    p5 = out["plain_5fold"]["blind_spot"]
    verdict = ("TRANSFERS_TO_BACTERIA" if bl.get("pass")
               else ("IN_DISTRIBUTION_ONLY" if p5.get("pass") else "DOES_NOT_TRANSFER"))
    res = {"artifact": "tb_supervised_vs_catalog", "schema": "tb-supervised-vs-catalog-v1",
           "date": str(_date.today()), "drug": DRUG, "gene": "rpoB",
           "cohort_n": len(uids), "R": int(y.sum()), "n_features": len(feats),
           "n_lineages": len(set(lineage.tolist())),
           "catalog_negative_n": int(neg.sum()), "catalog_negative_R": int(y[neg].sum()),
           "results": out, "verdict": verdict,
           "honest_note": ("Bacterial analog of the HIV test: supervised on ALL rpoB variants (not just WHO "
                           "grade-1/2) vs the WHO catalogue, free CRyPTIC measured phenotype, de-confounded by "
                           "LEAVE-ONE-LINEAGE-OUT (barcode lineage). A leave-lineage-out blind-spot PASS = the "
                           "technique transfers to bacteria; a collapse from plain-5fold to leave-lineage-out = "
                           "it learned lineage not mechanism (the user's caveat). rpoB-only target-gene scope; "
                           "efflux/promoter/copy-number resistance is invisible to point-variant features.")}
    out_path = REPO / "wiki" / f"tb_supervised_vs_catalog_{_date.today()}.json"
    out_path.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\ncohort n={len(uids)} R={int(y.sum())} | catalog-negative n={int(neg.sum())} "
          f"R={int(y[neg].sum())} | features={len(feats)} | lineages={len(set(lineage.tolist()))}")
    for k, v in out.items():
        b = v["blind_spot"]
        print(f"{k}: overall {v['overall_auroc']} | blind-spot n={b['n']} R={b['R']} AUROC={b.get('auroc')} "
              f"burden={b.get('burden')} null={b.get('null')} PASS={b['pass']}"
              + (f" [{b['note']}]" if b.get('note') else ""))
    print(f"VERDICT: {verdict}\n[wrote {out_path}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(force="--force" in sys.argv))
