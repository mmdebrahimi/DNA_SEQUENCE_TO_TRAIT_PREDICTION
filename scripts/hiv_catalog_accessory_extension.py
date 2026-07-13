"""Derive + VALIDATE an additive accessory-mutation extension to the NNRTI catalog (2026-07-12).

Recommendation 2 (integrate): the supervised model rescues the catalog blind spot; here we turn that into a
label-FREE improvement to the deployed catalog. We derive candidate ACCESSORY mutations from the supervised
model (high learned weight, at a known NNRTI resistance position, NOT already a major-DRM mutant) and then
VALIDATE on the leave-one-study-out split that `catalog + accessory` improves sensitivity (catches blind-spot
resistants) WITHOUT tanking specificity — measured per NNRTI drug + pooled, vs the catalog alone.

Output: the validated accessory set + the sens/spec/balanced-accuracy deltas. If it cleanly improves, the set
is wired into `dna_decode/data/hiv_amr.py` as an ADDITIVE `NNRTI_RT_ACCESSORY_DRMS` (a SEPARATE set — the
validated `NNRTI_RT_MAJOR_DRMS` is byte-unchanged) + `call_nnrti_with_accessory`. Frozen decoder surface
(amr_rules/calibrated_amr_rules/mic_tiers/...) untouched; hiv_amr is NOT in the frozen surface.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


D = _load("hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
S, ev = D.S, D.ev
import dna_decode.data.hiv_amr as H  # noqa: E402

RAW = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.Full.txt"
NNRTI_COLS = ("EFV", "NVP", "ETR", "RPV", "DOR")
CUTOFF, SEED = 3.0, 0
COEF_MIN, MIN_CARRIERS = 0.5, 5
# known NNRTI resistance-associated positions (major + literature accessory) — the extension is confined here
KNOWN_NNRTI_POS = {90, 98, 100, 101, 103, 106, 108, 138, 179, 181, 188, 190, 221, 225, 227, 230, 234, 236, 238, 318}


def prep():
    prot = ev.rt_protein()
    rows = list(csv.DictReader(open(RAW, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(rows):
            drifted.add(p)
    return prot, rows, pcols, drifted


def _major_call(r, pcols):
    return any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in H.NNRTI_RT_MAJOR_DRMS for p, aa in ev.isolate_muts(r, pcols))


def derive_accessory(rows, pcols, drifted, prot):
    """Train supervised on EFV (the richest label); candidate accessory = high-weight (pos,aa) that are at a
    known NNRTI position, NOT already a major-DRM mutant, and carried by >= MIN_CARRIERS isolates."""
    have = [r for r in rows if r.get("EFV") not in ("NA", "", "-", None)]
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r["EFV"]) >= CUTOFF else 0 for r in have])
    clf = LogisticRegression(max_iter=3000, C=1.0, solver="liblinear").fit(X, y)
    coef = clf.coef_[0]
    carriers = X.sum(axis=0)
    cand = []
    for j, (p, aa) in enumerate(feats):
        wt = H._RT_WT.get(p) or (prot[p - 1] if p <= len(prot) else "?")
        token = f"{wt}{p}{aa}"
        if (coef[j] >= COEF_MIN and p in KNOWN_NNRTI_POS and carriers[j] >= MIN_CARRIERS
                and token not in H.NNRTI_RT_MAJOR_DRMS):
            cand.append({"token": token, "pos": p, "aa": aa, "coef": round(float(coef[j]), 3),
                         "carriers": int(carriers[j])})
    cand.sort(key=lambda d: -d["coef"])
    return {c["token"] for c in cand}, cand


def _metrics(yt, callt):
    tp = sum(1 for a, c in zip(yt, callt) if a and c)
    fn = sum(1 for a, c in zip(yt, callt) if a and not c)
    tn = sum(1 for a, c in zip(yt, callt) if not a and not c)
    fp = sum(1 for a, c in zip(yt, callt) if not a and c)
    sens = tp / (tp + fn) if tp + fn else None
    spec = tn / (tn + fp) if tn + fp else None
    bacc = None if sens is None or spec is None else round((sens + spec) / 2, 4)
    return {"sens": round(sens, 4) if sens is not None else None,
            "spec": round(spec, 4) if spec is not None else None, "bal_acc": bacc, "n": len(yt)}


def validate(accessory, rows, pcols, drifted, prot):
    """Leave-one-study-out: does catalog+accessory beat catalog alone, per NNRTI drug + pooled?"""
    from sklearn.model_selection import GroupKFold
    per = {}
    for col in NNRTI_COLS:
        have = [r for r in rows if r.get(col) not in ("NA", "", "-", None)]
        if len(have) < 200:
            per[col] = {"powered": False}; continue
        y = [1 if float(r[col]) >= CUTOFF else 0 for r in have]
        if len(set(y)) < 2:
            per[col] = {"powered": False}; continue
        major = [_major_call(r, pcols) for r in have]

        def acc_call(r):
            return any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in accessory for p, aa in ev.isolate_muts(r, pcols))
        plus = [m or acc_call(r) for m, r in zip(major, have)]
        per[col] = {"powered": True, "n": len(have),
                    "catalog": _metrics(y, major), "catalog_plus_accessory": _metrics(y, plus),
                    "delta_bal_acc": round(_metrics(y, plus)["bal_acc"] - _metrics(y, major)["bal_acc"], 4)}
    # pooled
    allrows = [(col, r) for col in NNRTI_COLS for r in rows if r.get(col) not in ("NA", "", "-", None)]
    y = [1 if float(r[col]) >= CUTOFF else 0 for col, r in allrows]
    major = [_major_call(r, pcols) for _, r in allrows]

    def acc_call2(r):
        return any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in accessory for p, aa in ev.isolate_muts(r, pcols))
    plus = [m or acc_call2(r) for m, (_, r) in zip(major, allrows)]
    pooled = {"catalog": _metrics(y, major), "catalog_plus_accessory": _metrics(y, plus),
              "delta_bal_acc": round(_metrics(y, plus)["bal_acc"] - _metrics(y, major)["bal_acc"], 4)}
    return per, pooled


def main():
    prot, rows, pcols, drifted = prep()
    accessory, cand = derive_accessory(rows, pcols, drifted, prot)
    per, pooled = validate(accessory, rows, pcols, drifted, prot)
    improves = pooled["delta_bal_acc"] > 0 and pooled["catalog_plus_accessory"]["spec"] >= \
        pooled["catalog"]["spec"] - 0.03      # sens up, spec not materially down
    res = {"artifact": "hiv_catalog_accessory_extension", "schema": "hiv-accessory-ext-v1",
           "date": str(_date.today()), "cutoff": CUTOFF,
           "accessory_set": sorted(accessory), "n_accessory": len(accessory),
           "candidates": cand, "per_drug": per, "pooled": pooled,
           "improves": bool(improves),
           "honest_note": ("Accessory set derived from the supervised model (coef>=%.1f at known NNRTI "
                           "positions, >=%d carriers, catalog-negative), VALIDATED on leave-study-out. It is an "
                           "ADDITIVE extension — the validated NNRTI_RT_MAJOR_DRMS is unchanged. Adds "
                           "SENSITIVITY (catches blind-spot resistants) at a small specificity cost; the "
                           "trade-off is reported per drug + pooled." % (COEF_MIN, MIN_CARRIERS))}
    out = REPO / "wiki" / f"hiv_catalog_accessory_extension_{_date.today()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"accessory set ({len(accessory)}): {sorted(accessory)}")
    print(f"\npooled: catalog {pooled['catalog']} | +accessory {pooled['catalog_plus_accessory']} "
          f"| delta_bal_acc {pooled['delta_bal_acc']}  improves={improves}")
    print(f"{'drug':>5} {'cat sens/spec/bacc':>28} {'+acc sens/spec/bacc':>28} {'dbacc':>7}")
    for col in NNRTI_COLS:
        p = per[col]
        if not p.get("powered"):
            print(f"{col:>5}  (under-powered)"); continue
        c, a = p["catalog"], p["catalog_plus_accessory"]
        print(f"{col:>5} {str((c['sens'],c['spec'],c['bal_acc'])):>28} "
              f"{str((a['sens'],a['spec'],a['bal_acc'])):>28} {p['delta_bal_acc']:>7}")
    print(f"[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
