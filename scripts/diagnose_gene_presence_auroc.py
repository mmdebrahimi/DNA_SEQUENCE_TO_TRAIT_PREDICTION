"""Diagnostic: gene-presence + XGBoost LOSO returns AUROC=0.000 at N=12.

Reproduces the smoke-gate gene-presence variant per-fold and prints:
- Training set R/S balance per fold
- Test strain identity + true label
- Feature matrix shape (N_train, n_features) per fold
- Train feature density (fraction non-zero)
- Test row density
- Test-vocab-overlap (fraction of test strain's genes also in train vocab)
- XGBoost predicted probability per test strain
- Symmetric two-value pattern check

Goal: distinguish between (a) anti-predictive symmetric output (similar to
the calibration bug — would suggest a deeper plumbing issue), (b) high-dim
overfitting (XGBoost memorizes training fold and predicts uniform value on
held-out), or (c) test rows mostly all-zero (features so train-vocab-biased
that test strains look like empty rows to the classifier).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from dna_decode.data.annotations import parse_gff3
from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import gff_path
from dna_decode.models.classical_baselines import build_gene_presence_matrix
from dna_decode.models.classifiers import (
    ClassifierTrainingError,
    predict_proba,
    train_xgboost_classifier,
)


def _extract_gene_ids(annotations_table) -> set[str]:
    cds = annotations_table[annotations_table["type"] == "CDS"]
    gene_ids: set[str] = set()
    for _, row in cds.iterrows():
        gid = row.get("gene_id") or row.get("locus_tag")
        if gid:
            gene_ids.add(str(gid))
    return gene_ids


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_mini_cohort.parquet"))
    parser.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    parser.add_argument("--drug", default="ciprofloxacin")
    args = parser.parse_args(argv)

    cohort = load_cohort(args.cohort)
    drug_lower = args.drug.lower()

    strain_genes: dict[str, set[str]] = {}
    labels_dict: dict[str, int] = {}
    print(f"[diag] Parsing GFF3 for {len(cohort.strains)} strains...")
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        gp = gff_path(s.assembly_accession, args.refseq_cache)
        ann = parse_gff3(gp)
        strain_genes[s.strain_id] = _extract_gene_ids(ann)
        labels_dict[s.strain_id] = s.ast_labels[drug_lower]

    strain_order = sorted(strain_genes.keys())
    y = np.array([labels_dict[sid] for sid in strain_order], dtype=int)
    n = len(strain_order)

    print(f"[diag] N strains = {n}; R/S balance = {int((y==1).sum())}R / {int((y==0).sum())}S")
    print(f"[diag] gene-set sizes per strain (sorted desc): "
          f"{sorted([len(strain_genes[sid]) for sid in strain_order], reverse=True)[:5]} ...")

    all_y_score: list[float] = []
    all_y_true: list[int] = []

    print()
    print(f"{'fold':>4} {'held_id':<35} {'lbl':>3} {'tr R/S':>7} "
          f"{'n_feat':>7} {'tr_dens':>8} {'te_dens':>8} {'te_in_vocab%':>12} {'p(R)':>7}")
    print("-" * 100)

    for i in range(n):
        held_sid = strain_order[i]
        held_y = int(y[i])
        train_idx = [j for j in range(n) if j != i]
        train_gene_sets = [strain_genes[strain_order[j]] for j in train_idx]
        train_y = y[train_idx]
        test_gene_set = strain_genes[held_sid]

        X_train, vocab = build_gene_presence_matrix(train_gene_sets, gene_vocabulary=None)
        X_test, _ = build_gene_presence_matrix([test_gene_set], gene_vocabulary=vocab)

        n_feat = len(vocab)
        vocab_set = set(vocab)
        te_in_vocab = sum(1 for g in test_gene_set if g in vocab_set)
        te_in_vocab_frac = te_in_vocab / max(1, len(test_gene_set))
        tr_density = float(X_train.mean())
        te_density = float(X_test.mean())

        try:
            clf = train_xgboost_classifier(X_train, train_y, drug_name=args.drug, calibrate=False)
            p_R = float(predict_proba(clf, X_test)[0])
        except ClassifierTrainingError as e:
            p_R = float(train_y.mean())
            print(f"  fold {i}: training error: {e}; fell back to mean = {p_R:.3f}")

        all_y_score.append(p_R)
        all_y_true.append(held_y)

        tr_R = int((train_y == 1).sum())
        tr_S = int((train_y == 0).sum())

        print(f"{i:>4} {held_sid:<35} {held_y:>3} "
              f"{tr_R:>2}R/{tr_S:>2}S {n_feat:>7} {tr_density:>8.4f} {te_density:>8.4f} "
              f"{te_in_vocab_frac:>11.2%} {p_R:>7.3f}")

    print("-" * 100)
    y_true_arr = np.array(all_y_true)
    y_score_arr = np.array(all_y_score)
    print(f"\n[diag] Unique probability values returned: {sorted(set(round(s, 6) for s in all_y_score))}")
    print(f"[diag] Distribution of p(R) by true label:")
    for lbl_val, lbl_name in [(0, "S=0"), (1, "R=1")]:
        scores = y_score_arr[y_true_arr == lbl_val]
        if len(scores):
            print(f"   {lbl_name}: n={len(scores)}, mean={scores.mean():.4f}, scores={list(np.round(scores, 4))}")

    # Compute AUROC manually
    from sklearn.metrics import roc_auc_score
    try:
        auroc = roc_auc_score(y_true_arr, y_score_arr)
        print(f"\n[diag] AUROC = {auroc:.4f}")
    except Exception as e:
        print(f"\n[diag] AUROC failed: {e}")

    # Check anti-predictive pattern: if we INVERT predictions, what's AUROC?
    try:
        auroc_inv = roc_auc_score(y_true_arr, 1.0 - y_score_arr)
        print(f"[diag] AUROC with INVERTED scores = {auroc_inv:.4f}")
    except Exception as e:
        print(f"[diag] AUROC (inverted) failed: {e}")

    # Predicted-class-only AUROC (treat p>=0.5 as class 1)
    pred_class = (y_score_arr >= 0.5).astype(int)
    print(f"[diag] Class predictions at threshold 0.5: {list(pred_class)}")
    print(f"[diag] True labels:                       {list(y_true_arr.astype(int))}")
    print(f"[diag] Agreement: {(pred_class == y_true_arr).sum()}/{len(y_true_arr)}")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
