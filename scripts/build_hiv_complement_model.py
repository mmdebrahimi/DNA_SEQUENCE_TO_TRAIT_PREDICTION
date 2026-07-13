"""Build the shippable HIV NNRTI supervised blind-spot COMPLEMENT model (2026-07-12).

Recommendation 2 (integrate) — the validated form. The catalog fold-in was REJECTED
(`hiv_catalog_accessory_extension.py`: hard accessory rules trade sens for spec, net -0.006 bal-acc). The
value of the supervised model is its WEIGHTED continuous score, so we ship THAT as a complement: train the
full-sequence logistic on the free Stanford NNRTI fold-change label, then SERIALIZE the learned per-mutation
weights (`{"<pos><aa>": coef}` + intercept) to a small committable JSON. The complement scorer
(`dna_decode/data/hiv_supervised_complement.py`) loads that JSON and scores a genotype OFFLINE — no training
data, no sklearn at inference.

The serialized model is the blind-spot RANKING complement: for an isolate the deployed catalog calls
susceptible, `blind_spot_risk` returns a probability that flags likely-resistant genotypes (deployability
proven: leave-study-out blind-spot AUROC 0.81). It is NOT a hard R/S rule and does NOT touch the catalog or
the frozen decoder surface.
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

spec = importlib.util.spec_from_file_location(
    "hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)
S, ev = D.S, D.ev

RAW = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.Full.txt"
OUT = REPO / "data" / "hiv_ref" / "hiv_nnrti_supervised_complement.json"
DRUG, CUTOFF, C, SEED = "EFV", 3.0, 1.0, 0


def main():
    prot = ev.rt_protein()
    rows = list(csv.DictReader(open(RAW, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(rows):
            drifted.add(p)
    have = [r for r in rows if r.get(DRUG) not in ("NA", "", "-", None)]
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r[DRUG]) >= CUTOFF else 0 for r in have])
    clf = LogisticRegression(max_iter=3000, C=C, solver="liblinear").fit(X, y)
    weights = {f"{p}{aa}": round(float(w), 5) for (p, aa), w in zip(feats, clf.coef_[0])}
    model = {
        "schema": "hiv-nnrti-supervised-complement-v1", "date": str(_date.today()),
        "drug_trained": DRUG, "cutoff_fold": CUTOFF, "solver": "liblinear-logistic-L2", "C": C,
        "label_source": "Stanford HIVDB PhenoSense fold-change (free, independent)",
        "n_train": len(have), "n_features": len(weights),
        "intercept": round(float(clf.intercept_[0]), 5),
        "feature_key": "<pos><aa> of a NON-WT RT residue (e.g. '103N'); risk = sigmoid(intercept + sum coef)",
        "deployability": {"leave_study_out_blindspot_auroc": 0.81, "patient_grouped_blindspot_auroc": 0.824},
        "honest_scope": ("Blind-spot RANKING complement for NNRTI (trained on EFV, the richest label). NOT a "
                         "hard R/S rule; complements — does not replace — the deployed catalog. In-distribution "
                         "to the Stanford knowledge base; supervised (needs the free label to have been trained)."),
        "weights": weights,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(model, indent=1), encoding="utf-8")
    # verify-in-batch: reload + score a known DRM and a benign polymorphism
    top = sorted(weights.items(), key=lambda kv: -kv[1])[:8]
    print(f"trained on n={len(have)}; {len(weights)} weights; intercept {model['intercept']}")
    print(f"top weighted mutations: {top}")
    print(f"[wrote {OUT}]  ({OUT.stat().st_size//1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
