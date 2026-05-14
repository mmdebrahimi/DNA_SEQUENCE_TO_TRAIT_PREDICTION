"""Shared LOSO runners for k-mer + NT/k-mer fusion classifiers.

Order-explicit API per `plans/Stage1_Refactor_And_Test_Hardening_Plan.md` Step 2.1.
Caller supplies `strain_ids: list[str]` to lock fold order; the runner respects it
verbatim with no internal sorting. This prevents the alignment-drift risk that
prompted the /brainstorm Round 1 critique (Stage 1 vs smoke gate had different
subsets/orders for the same logical comparator).

Used by both `scripts/smoke_gate_12strain_cipro.py` and
`scripts/stage1_n40_cipro.py`; returns `CVResult` so the strain_ids alignment
property is available downstream.

Re-raises `ClassifierTrainingError` rather than mean-fallback — silent fallback
masks the exact failures these runners are supposed to surface.
"""
from __future__ import annotations

import numpy as np

from dna_decode.eval.cv import CVResult, FoldResult
from dna_decode.models.classical_baselines import (
    CONTIG_SEPARATOR,
    build_kmer_vocabulary,
    kmers_to_feature_matrix,
)
from dna_decode.models.classifiers import (
    ClassifierTrainingError,
    predict_proba,
    train_xgboost_classifier,
)


def run_kmer_xgboost_loso(
    seqs_by_strain: dict[str, str],
    labels_by_strain: dict[str, int],
    strain_ids: list[str],
    drug: str,
    k: int = 8,
    top_n: int = 10_000,
) -> CVResult:
    """k-mer + XGBoost LOSO with WITHIN-FOLD vocabulary rebuild.

    Args:
        seqs_by_strain: strain_id -> already-N-concatenated genome string. Caller
            is responsible for concatenation using CONTIG_SEPARATOR.
        labels_by_strain: strain_id -> binary R/S (1/0) label.
        strain_ids: ordered list of strain IDs to use as LOSO folds. The runner
            does NOT sort or filter this list — fold order = input order.
        drug: drug name for `train_xgboost_classifier` reporting + CVResult.drug.
        k: k-mer length.
        top_n: vocabulary cap (top-N most-frequent k-mers from the training set).

    Returns:
        CVResult with `strategy="loso"`, one FoldResult per strain. `held_out_id`
        is the strain_id; `y_true` / `y_score` are length-1 arrays.

    Raises:
        ClassifierTrainingError: if any per-fold XGBoost training fails. NOT
            silently caught — Stage 1 needs to surface real bugs rather than
            mask them with a mean-fallback (per Step 2.4).
    """
    result = CVResult(strategy="loso", drug=drug)
    n = len(strain_ids)
    for i, held in enumerate(strain_ids):
        train_idx = [j for j in range(n) if j != i]
        train_strains = [strain_ids[j] for j in train_idx]
        train_seqs = [seqs_by_strain[s] for s in train_strains]
        train_y = np.array([labels_by_strain[s] for s in train_strains], dtype=int)
        test_seq = seqs_by_strain[held]
        test_y = int(labels_by_strain[held])

        vocab = build_kmer_vocabulary(train_seqs, k=k, top_n=top_n)
        X_train = kmers_to_feature_matrix(train_seqs, vocab, k=k)
        X_test = kmers_to_feature_matrix([test_seq], vocab, k=k)

        clf = train_xgboost_classifier(X_train, train_y, drug_name=drug, calibrate=False)
        score = float(predict_proba(clf, X_test)[0])

        result.folds.append(
            FoldResult(
                held_out_id=held,
                held_out_indices=[i],
                train_indices=train_idx,
                y_true=np.array([test_y], dtype=int),
                y_score=np.array([score], dtype=np.float32),
                n_train=len(train_idx),
                n_test=1,
            )
        )
    return result


def run_fusion_loso(
    X_nt: np.ndarray,
    seqs_by_strain: dict[str, str],
    labels_by_strain: dict[str, int],
    strain_ids: list[str],
    drug: str,
    k: int = 8,
    top_n: int = 10_000,
) -> CVResult:
    """NT + k-mer concat features → logistic regression LOSO.

    Diagnostic only at Stage 1; not gate-bearing. Caller MUST pass an `X_nt`
    matrix whose row i corresponds to strain_ids[i] (alignment contract).

    Args:
        X_nt: (len(strain_ids), nt_dim) NT mean-pooled embeddings, row-aligned
            with `strain_ids`.
        seqs_by_strain: strain_id -> N-concatenated genome string.
        labels_by_strain: strain_id -> binary R/S (1/0) label.
        strain_ids: ordered list of strain IDs (LOSO fold order; matches X_nt rows).
        drug: drug name (reporting).
        k: k-mer length.
        top_n: vocabulary cap.

    Returns:
        CVResult with `strategy="loso"`, one FoldResult per strain.

    Raises:
        ClassifierTrainingError: passthrough from `train_xgboost_classifier`-style
            logreg failures. The fusion path uses raw LogisticRegression so a
            small-N single-class fold would surface here.
    """
    from sklearn.linear_model import LogisticRegression

    if X_nt.shape[0] != len(strain_ids):
        raise ValueError(
            f"X_nt.shape[0]={X_nt.shape[0]} != len(strain_ids)={len(strain_ids)}"
        )

    result = CVResult(strategy="loso", drug=drug)
    n = len(strain_ids)
    for i, held in enumerate(strain_ids):
        train_idx = [j for j in range(n) if j != i]
        train_strains = [strain_ids[j] for j in train_idx]
        train_seqs = [seqs_by_strain[s] for s in train_strains]
        train_y = np.array([labels_by_strain[s] for s in train_strains], dtype=int)
        test_seq = seqs_by_strain[held]
        test_y = int(labels_by_strain[held])

        vocab = build_kmer_vocabulary(train_seqs, k=k, top_n=top_n)
        Xk_tr = kmers_to_feature_matrix(train_seqs, vocab, k=k)
        Xk_te = kmers_to_feature_matrix([test_seq], vocab, k=k)
        X_tr = np.hstack([X_nt[train_idx], Xk_tr]).astype(np.float32)
        X_te = np.hstack([X_nt[[i]], Xk_te]).astype(np.float32)

        clf = LogisticRegression(max_iter=2000, solver="liblinear", C=1.0, random_state=42)
        clf.fit(X_tr, train_y)
        score = float(clf.predict_proba(X_te)[0, 1])

        result.folds.append(
            FoldResult(
                held_out_id=held,
                held_out_indices=[i],
                train_indices=train_idx,
                y_true=np.array([test_y], dtype=int),
                y_score=np.array([score], dtype=np.float32),
                n_train=len(train_idx),
                n_test=1,
            )
        )
    return result
