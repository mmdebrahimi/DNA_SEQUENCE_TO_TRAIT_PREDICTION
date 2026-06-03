"""Interpret the winning presence/absence k-mer model + probe the study confound (EP-4).

Context: full-spectrum k-mer PRESENCE/ABSENCE LOSO = 0.729 on ExPEC(Salipante) vs
EPEC(Hazen). Before building the v0 gene-targeted resolver on top of that result we
must know what the 0.729 actually IS.

STRUCTURAL CONFOUND (discovered by inspecting the substrate): in this 24-genome
subset study == class EXACTLY (all 12 ExPEC are Salipante, all 12 EPEC are Hazen).
So no feature can separate biology from study/assembler batch here — the two
partitions are identical. We CANNOT resolve the confound on this subset; we can only
characterize the SHAPE of the signal, which tells us which hopeful/worried regime
we're in:

  - SPARSE  (a few k-mers recover ~0.729; importance highly concentrated) -> the
    discrimination rests on a handful of accessory loci. Consistent with real
    virulence-gene presence (the v0-resolver-friendly, hopeful case) — though still
    study-confounded until validated cross-study.
  - DIFFUSE (need hundreds of k-mers; importance spread genome-wide) -> discrimination
    is whole-genome composition = lineage/assembler batch. Worrying: the 0.729 is
    likely batch, and a gene-targeted resolver won't inherit it.

Outputs: importance concentration, an honest within-fold top-K selection AUROC curve
(no leakage — top-K chosen inside each LOSO fold), and the top-k-mer per-class presence
table. NO external deps beyond the existing bake-off helpers + sklearn/xgboost.
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score

from scripts.pathotype_kmer_bakeoff import (
    load_strains, precompute_counts, _vocab_from_cache, _matrix_from_cache, F1, WGS_MASTER,
)
from dna_decode.models.classifiers import predict_proba, train_xgboost_classifier

K = 8
FULL_TOP_N = 4 ** K
TOPK_GRID = [10, 25, 50, 100, 250]


def study_of(ids: list[str]) -> dict[str, str]:
    """Map strain id -> originating study (Source col) from the F1 metadata."""
    rows = {r["ID"]: r["Source"].split()[0] for r in csv.DictReader(open(F1, encoding="utf-8"))}
    return {sid: rows.get(sid.split("|")[0], "?") for sid in ids}


def binary_matrix(strains, cache, vocab):
    return (_matrix_from_cache(strains, vocab, cache) > 0).astype(np.float32)


def loso_full(X, y, drug):
    n = len(y); proba = np.zeros(n)
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        clf = train_xgboost_classifier(X[tr], y[tr], drug_name=drug, calibrate=False)
        proba[i] = float(predict_proba(clf, X[i:i+1])[0])
    return float(roc_auc_score(y, proba)) if len(set(y.tolist())) == 2 else float("nan")


def loso_topk_within_fold(X, y, drug, k_top):
    """Honest top-K: select the K most-important k-mers INSIDE each fold (no leakage)."""
    n = len(y); proba = np.zeros(n)
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        full = train_xgboost_classifier(X[tr], y[tr], drug_name=drug, calibrate=False)
        imp = full.model.feature_importances_
        top = np.argsort(imp)[::-1][:k_top]
        sub = train_xgboost_classifier(X[tr][:, top], y[tr], drug_name=drug, calibrate=False)
        proba[i] = float(predict_proba(sub, X[i:i+1][:, top])[0])
    return float(roc_auc_score(y, proba)) if len(set(y.tolist())) == 2 else float("nan")


def main() -> int:
    seqs, labels, ids = load_strains()
    y = np.array([labels[s] for s in ids], dtype=int)
    study = study_of(ids)
    # --- structural confound check: study vs class ---
    by = {}
    for s in ids:
        by.setdefault((labels[s], study[s]), 0)
        by[(labels[s], study[s])] += 1
    study_eq_class = all(
        (labels[s] == 0) == (study[s].startswith("Salipante")) for s in ids
    )
    print(f"[interp] {len(ids)} genomes; class x study cells = {by}")
    print(f"[interp] study == class (structural confound): {study_eq_class}")

    print("[interp] building full-spectrum presence/absence matrix...", flush=True)
    cache = precompute_counts(seqs, k=K)
    vocab = _vocab_from_cache(ids, cache, FULL_TOP_N)
    X = binary_matrix(ids, cache, vocab)
    print(f"[interp] X = {X.shape} (binary presence/absence)")

    # --- all-data model: importance concentration ---
    full = train_xgboost_classifier(X, y, drug_name="expec_vs_epec_pa", calibrate=False)
    imp = full.model.feature_importances_
    order = np.argsort(imp)[::-1]
    nz = int((imp > 0).sum())
    csum = np.cumsum(imp[order]) / imp.sum()
    def kmers_to(frac):
        return int(np.searchsorted(csum, frac) + 1)
    conc = {f"kmers_to_{int(f*100)}pct": kmers_to(f) for f in (0.5, 0.8, 0.95)}
    print(f"[interp] non-zero-importance k-mers: {nz} / {X.shape[1]}; concentration {conc}")

    # --- honest within-fold top-K AUROC curve ---
    auroc_full = loso_full(X, y, "expec_vs_epec_pa")
    print(f"[interp] full presence/absence LOSO AUROC = {auroc_full:.4f} (should ~match 0.729)")
    curve = {}
    for k_top in TOPK_GRID:
        a = loso_topk_within_fold(X, y, "expec_vs_epec_pa", k_top)
        curve[k_top] = a
        print(f"[interp] within-fold top-{k_top:<4d} LOSO AUROC = {a:.4f}", flush=True)

    # --- top-10 k-mer per-class presence table (descriptive markers) ---
    top10 = []
    for idx in order[:10]:
        km = vocab[idx]
        present_expec = int(X[y == 0, idx].sum())
        present_epec = int(X[y == 1, idx].sum())
        top10.append({"kmer": km, "importance": float(imp[idx]),
                      "in_expec_n": present_expec, "in_epec_n": present_epec})
    print("[interp] top-10 discriminative k-mers (presence by class, /12 each):")
    for t in top10:
        print(f"    {t['kmer']}  imp={t['importance']:.3f}  ExPEC={t['in_expec_n']}/12  EPEC={t['in_epec_n']}/12")

    # --- verdict: sparse (gene-like) vs diffuse (lineage/batch) ---
    best_small = max(curve[k] for k in (10, 25, 50) if k in curve)
    if best_small >= auroc_full - 0.03:
        shape = "SPARSE: <=50 k-mers recover the full signal -> a few accessory loci drive it (gene-presence-like; v0-resolver-friendly, hopeful)"
    elif curve.get(250, 0) >= auroc_full - 0.03 and best_small < auroc_full - 0.10:
        shape = "DIFFUSE: needs hundreds of k-mers -> whole-genome composition = lineage/assembler batch-like (worrying for biology)"
    else:
        shape = "INTERMEDIATE: signal partly concentrated; inconclusive shape on N=24"

    res = {
        "contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "n": len(ids),
        "representation": "kmer_k8_presence_absence_full_spectrum",
        "study_equals_class": study_eq_class,
        "class_x_study_cells": {f"{k[0]}|{k[1]}": v for k, v in by.items()},
        "confound_note": "study == class on this subset -> biology and batch are the SAME partition; the 0.729 cannot be attributed to biology here. Need within-study OR lineage/ST-matched genomes to resolve.",
        "full_pa_loso_auroc": auroc_full,
        "nonzero_importance_kmers": nz, "importance_concentration": conc,
        "within_fold_topk_auroc": curve,
        "top10_kmers": top10,
        "signal_shape_verdict": shape,
        "next_action": "Obtain ExPEC+EPEC from a SINGLE study (or ST/lineage-matched cross-study pairs) to break study==class; only then is an AUROC interpretable as pathotype biology. In parallel, the v0 deterministic virulence-gene-cluster resolver (ledger-locked) is the interpretable path that does not depend on this confounded learned signal.",
        "ids": ids,
    }
    out = REPO / "research_outputs/pathotype_model_interpret_confound_2026-06-02.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[interp] SIGNAL SHAPE: {shape}")
    print(f"[interp] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
