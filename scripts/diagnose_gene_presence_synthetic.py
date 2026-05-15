"""Synthetic-data check of the LOSO base-rate-inversion hypothesis.

Builds an (N=12, ~5000) sparse binary feature matrix where features carry NO
class signal (random Bernoulli, identical distribution per class). Runs the
same LOSO loop the smoke gate uses. If the hypothesis is correct:

  - Each fold's prediction lands ~train_y.mean() (≈0.455 or 0.545)
  - The two possible predictions are SYMMETRIC around 0.5
  - LOSO rank inverts → AUROC=0.000 or close to it

If we see AUROC≈0.5 instead, the hypothesis is wrong and the real-data bug
needs deeper investigation against the actual GFF3 features.
"""
from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score

from dna_decode.models.classifiers import predict_proba, train_xgboost_classifier


def main() -> None:
    rng = np.random.default_rng(42)
    n_strains = 12
    n_features = 5000
    # Balanced 6R/6S
    y = np.array([1] * 6 + [0] * 6, dtype=int)
    rng.shuffle(y)
    # No-signal features: each gene present in any strain with probability 0.6
    # (mimicking E. coli core-genome high gene overlap)
    X = (rng.random((n_strains, n_features)) < 0.6).astype(np.float32)

    print(f"Synthetic feature matrix: {X.shape}, density={float(X.mean()):.3f}")
    print(f"Labels: {list(y)}")
    print()

    all_y_true: list[int] = []
    all_y_score: list[float] = []
    print(f"{'fold':>4} {'true':>4} {'tr_R/S':>7} {'tr_mean':>8} {'p(R)':>7}")
    for i in range(n_strains):
        train_idx = [j for j in range(n_strains) if j != i]
        X_tr = X[train_idx]
        y_tr = y[train_idx]
        clf = train_xgboost_classifier(X_tr, y_tr, drug_name="synth", calibrate=False)
        p_R = float(predict_proba(clf, X[[i]])[0])
        all_y_true.append(int(y[i]))
        all_y_score.append(p_R)
        print(f"{i:>4} {int(y[i]):>4} "
              f"{int((y_tr==1).sum()):>2}R/{int((y_tr==0).sum()):>2}S "
              f"{float(y_tr.mean()):>8.4f} {p_R:>7.4f}")

    y_t = np.array(all_y_true)
    y_s = np.array(all_y_score)
    print()
    print(f"Unique p(R) values: {sorted(set(round(s, 4) for s in all_y_score))}")
    print(f"AUROC = {roc_auc_score(y_t, y_s):.4f}")
    print(f"AUROC (inverted) = {roc_auc_score(y_t, 1.0 - y_s):.4f}")


if __name__ == "__main__":
    main()
