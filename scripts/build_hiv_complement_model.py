"""Build the shippable HIV supervised blind-spot COMPLEMENT models — NNRTI + PI + INSTI (2026-07-12).

Recommendation 2 (integrate), extended across HIV genes. The catalog fold-in was REJECTED (hard rules trade
sens for spec); the value is the WEIGHTED continuous score, so we ship THAT: train the full-sequence logistic
on the free Stanford fold-change label for a strong well-powered drug per class, then SERIALIZE the learned
per-mutation weights (`{"<pos><aa>": coef}` + intercept) to a small committable JSON. The complement scorer
(`dna_decode/data/hiv_supervised_complement.py`) loads a class's JSON and scores a genotype OFFLINE — no
training data, no sklearn at inference.

Classes (deployment-validated blind-spot rescue): NNRTI (RT, train EFV — leave-study-out 0.81), PI (protease,
train LPV — 0.89), INSTI (integrase, train RAL — 0.89). Each is a RANKING complement, NOT a hard R/S rule and
NOT a catalog replacement. Frozen decoder surface + hiv_amr catalog untouched.
"""
from __future__ import annotations

import argparse
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

RAW = REPO / "data" / "raw" / "hiv"
REF = REPO / "data" / "hiv_ref"
CUTOFF, C = 3.0, 1.0

# class -> (dataset, training drug, protein source, deployability blind-spot AUROC, output JSON)
CLASSES = {
    "NNRTI": {"dataset": "NNRTI_DataSet.Full.txt", "drug": "EFV", "protein": "rt",
              "leave_study_out": 0.81, "out": "hiv_nnrti_supervised_complement.json"},
    "PI": {"dataset": "PI_DataSet.Full.txt", "drug": "LPV", "protein": "HIV1_PR_HXB2_cds.fna",
           "leave_study_out": 0.89, "out": "hiv_pi_supervised_complement.json"},
    "INSTI": {"dataset": "INI_DataSet.Full.txt", "drug": "RAL", "protein": "HIV1_IN_HXB2_cds.fna",
              "leave_study_out": 0.89, "out": "hiv_insti_supervised_complement.json"},
}


def _protein(src):
    if src == "rt":
        return ev.rt_protein()
    seq = "".join(l.strip() for l in open(REF / src) if not l.startswith(">"))
    return "".join(ev.CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def build(cls: str) -> Path:
    cfg = CLASSES[cls]
    prot = _protein(cfg["protein"])
    rows = list(csv.DictReader(open(RAW / cfg["dataset"], encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(rows):
            drifted.add(p)
    drug = cfg["drug"]
    have = [r for r in rows if r.get(drug) not in ("NA", "", "-", None)]
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r[drug]) >= CUTOFF else 0 for r in have])
    clf = LogisticRegression(max_iter=3000, C=C, solver="liblinear").fit(X, y)
    weights = {f"{p}{aa}": round(float(w), 5) for (p, aa), w in zip(feats, clf.coef_[0])}
    model = {
        "schema": "hiv-supervised-complement-v1", "drug_class": cls, "date": str(_date.today()),
        "gene": {"NNRTI": "RT", "PI": "PR", "INSTI": "IN"}[cls],
        "drug_trained": drug, "cutoff_fold": CUTOFF, "solver": "liblinear-logistic-L2", "C": C,
        "label_source": "Stanford HIVDB PhenoSense fold-change (free, independent)",
        "n_train": len(have), "n_features": len(weights), "intercept": round(float(clf.intercept_[0]), 5),
        "feature_key": "<pos><aa> of a NON-WT residue (e.g. '103N'); risk = sigmoid(intercept + sum coef)",
        "deployability": {"leave_study_out_blindspot_auroc": cfg["leave_study_out"]},
        "honest_scope": (f"Blind-spot RANKING complement for {cls} (trained on {drug}). NOT a hard R/S rule; "
                         "complements — does not replace — the deployed catalog. In-distribution to Stanford; "
                         "supervised (needs the free label to have been trained)."),
        "weights": weights,
    }
    out = REF / cfg["out"]
    out.write_text(json.dumps(model, indent=1), encoding="utf-8")
    top = sorted(weights.items(), key=lambda kv: -kv[1])[:6]
    print(f"[{cls}] trained {drug} n={len(have)}; {len(weights)} weights; intercept {model['intercept']}; "
          f"top {top}; -> {out} ({out.stat().st_size//1024} KB)")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", choices=list(CLASSES) + ["all"], default="all")
    a = ap.parse_args(argv)
    for cls in (list(CLASSES) if a.cls == "all" else [a.cls]):
        build(cls)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
