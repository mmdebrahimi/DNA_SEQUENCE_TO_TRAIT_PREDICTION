"""Step 18 — Classical baseline benchmarks.

Foundation-model embeddings must beat classical baselines to justify the
foundation-model premise. Without this control, every Phase 1 result is
vulnerable to "would a simple k-mer baseline perform equally well?"

Three baselines:
  1. AMRFinder gene-call presence/absence + XGBoost
  2. k-mer (k=8) count vectors + logistic regression OR XGBoost
  3. Gene-family presence/absence (from Bakta annotations) + XGBoost

Each baseline reuses Step 10's CV harness for fair comparison against
foundation-model classifiers (Step 9).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any

import numpy as np

from dna_decode.models.classifiers import (
    ClassifierTrainingError,
    MIN_TRAINING_SAMPLES,
    TrainedClassifier,
    XGBParams,
)


DEFAULT_KMER_K = 8
DEFAULT_KMER_TOP_N = 10_000  # top-N most frequent k-mers

# Canonical N-run separator for contig concatenation. 100 N's is longer than any
# realistic k-mer (k<=32), so extract_kmer_counts skipping N-containing windows
# correctly prevents cross-contig k-mers. Promoted to module constant 2026-05-14
# (Stage1_Refactor_And_Test_Hardening_Plan Step 2.9) so loso_kmer + Stage 1 runner
# share one source of truth instead of three magic "N" * 100 strings.
CONTIG_SEPARATOR: str = "N" * 100


# ---- k-mer feature extraction ----


def extract_kmer_counts(
    sequence: str, k: int = DEFAULT_KMER_K, vocabulary: list[str] | None = None
) -> dict[str, int]:
    """Count k-mers in a sequence.

    If `vocabulary` is supplied, returns counts ONLY for those k-mers (zero
    for absent). Otherwise returns all observed k-mers.
    """
    counts: dict[str, int] = {}
    seq_upper = sequence.upper()
    for i in range(len(seq_upper) - k + 1):
        kmer = seq_upper[i : i + k]
        if "N" in kmer:
            continue  # skip ambiguous bases
        counts[kmer] = counts.get(kmer, 0) + 1
    if vocabulary is not None:
        return {kmer: counts.get(kmer, 0) for kmer in vocabulary}
    return counts


def build_kmer_vocabulary(
    sequences: list[str], k: int = DEFAULT_KMER_K, top_n: int = DEFAULT_KMER_TOP_N
) -> list[str]:
    """Select the top-N most frequent k-mers across a corpus.

    Reduces feature dimensionality from theoretical 4^k = 65,536 (k=8) to
    a manageable top-N for downstream classifiers.
    """
    aggregate: Counter[str] = Counter()
    for seq in sequences:
        seq_counts = extract_kmer_counts(seq, k=k)
        aggregate.update(seq_counts)
    return [kmer for kmer, _ in aggregate.most_common(top_n)]


def kmers_to_feature_matrix(
    sequences: list[str], vocabulary: list[str], k: int = DEFAULT_KMER_K
) -> np.ndarray:
    """Build (n_sequences, len(vocabulary)) count matrix."""
    n = len(sequences)
    m = len(vocabulary)
    out = np.zeros((n, m), dtype=np.float32)
    vocab_index = {kmer: i for i, kmer in enumerate(vocabulary)}
    for row_i, seq in enumerate(sequences):
        counts = extract_kmer_counts(seq, k=k, vocabulary=vocabulary)
        for kmer, count in counts.items():
            if kmer in vocab_index:
                out[row_i, vocab_index[kmer]] = count
    return out


# ---- Gene presence/absence feature extraction ----


def build_gene_presence_matrix(
    strain_gene_sets: list[set[str]],
    gene_vocabulary: list[str] | None = None,
) -> tuple[np.ndarray, list[str]]:
    """Convert per-strain gene-id sets into a (n_strains, n_genes) binary matrix.

    If `gene_vocabulary` not supplied, uses the union of all genes across
    strains.
    """
    if gene_vocabulary is None:
        all_genes = sorted(set().union(*strain_gene_sets) if strain_gene_sets else set())
    else:
        all_genes = list(gene_vocabulary)
    gene_index = {g: i for i, g in enumerate(all_genes)}

    n = len(strain_gene_sets)
    m = len(all_genes)
    out = np.zeros((n, m), dtype=np.float32)
    for row_i, gene_set in enumerate(strain_gene_sets):
        for g in gene_set:
            if g in gene_index:
                out[row_i, gene_index[g]] = 1.0
    return out, all_genes


# ---- Classifier training (shared logic) ----


def _train_baseline_xgboost(
    X: np.ndarray, y: np.ndarray, drug_name: str, params: XGBParams | None = None
) -> TrainedClassifier:
    """Train an XGBoost classifier on the given feature matrix.

    Wraps `dna_decode.models.classifiers.train_xgboost_classifier` so all
    baselines share the same calibration + error-handling discipline.
    """
    from dna_decode.models.classifiers import train_xgboost_classifier

    return train_xgboost_classifier(X, y, drug_name=drug_name, params=params, calibrate=True)


def _train_baseline_logreg(
    X: np.ndarray, y: np.ndarray, drug_name: str, *, calibrate: bool = True
) -> TrainedClassifier:
    """Train a logistic-regression classifier (for k-mer baseline).

    When `calibrate=False`, returns the raw `LogisticRegression` without
    `CalibratedClassifierCV` wrapping — used by Stage 1 to keep AUROC measurement
    a function of representation quality rather than calibration-wrapper behavior.
    Per `plans/Stage1_Refactor_And_Test_Hardening_Plan.md` Step 2.3. Calibration
    is a small-N footgun (see LESSONS_LEARNED 2026-05-14 calibration overcorrection).
    """
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.calibration import CalibratedClassifierCV
    except ImportError as e:
        raise ClassifierTrainingError(
            "scikit-learn not installed; run `uv sync` to install Phase 1 deps"
        ) from e

    if X.shape[0] < MIN_TRAINING_SAMPLES:
        raise ClassifierTrainingError(
            f"Need >= {MIN_TRAINING_SAMPLES} training samples; got {X.shape[0]}"
        )
    if len(set(y.tolist())) < 2:
        raise ClassifierTrainingError(
            f"Drug {drug_name!r}: training labels are single-class"
        )

    base = LogisticRegression(max_iter=1000, solver="liblinear", random_state=42)

    if not calibrate:
        # Raw logreg path — skip CalibratedClassifierCV entirely.
        base.fit(X, y)
        return TrainedClassifier(
            model=base,
            drug_name=drug_name,
            feature_dim=X.shape[1],
            calibrated=False,
        )

    # Wave 3.5 C7 fix: CV folds bounded by MINORITY class count
    minority_count = int(min((y == 1).sum(), (y == 0).sum()))
    if minority_count < 2:
        import warnings as _warnings
        _warnings.warn(
            f"Drug {drug_name!r}: minority class has {minority_count} sample(s); "
            f"skipping sigmoid calibration. Returning uncalibrated logistic regression.",
            RuntimeWarning,
            stacklevel=2,
        )
        base.fit(X, y)
        return TrainedClassifier(
            model=base,
            drug_name=drug_name,
            feature_dim=X.shape[1],
            calibrated=False,
        )
    cv_folds = max(2, min(3, minority_count))
    model = CalibratedClassifierCV(base, cv=cv_folds, method="sigmoid")
    model.fit(X, y)

    return TrainedClassifier(
        model=model,
        drug_name=drug_name,
        feature_dim=X.shape[1],
        calibrated=True,
    )


# ---- Three baseline entry points ----


def train_amrfinder_baseline(
    strain_amr_genes: dict[str, set[str]],
    labels: dict[str, int],
    drug_name: str,
    gene_vocabulary: list[str] | None = None,
) -> tuple[TrainedClassifier, list[str], list[str]]:
    """Train XGBoost on AMRFinder gene-call presence/absence vectors per strain.

    Args:
        strain_amr_genes: strain_id -> set of AMR gene symbols detected by
            AMRFinder on that strain (sourced from Step 4 `ResistanceCatalog`
            cross-referenced with the strain's annotations).
        labels: strain_id -> binary R/S label for the target drug.
        drug_name: for downstream attribution + reporting.
        gene_vocabulary: optional pre-computed vocabulary (e.g., for held-out
            test sets); defaults to the union across the training strains.

    Returns:
        (TrainedClassifier, gene_vocabulary, strain_id_order). The gene
        vocabulary and strain order let callers reconstruct the feature
        matrix for prediction.
    """
    strain_ids = sorted(strain_amr_genes.keys() & labels.keys())
    if not strain_ids:
        raise ClassifierTrainingError(
            f"AMRFinder baseline for {drug_name!r}: no strain overlap between "
            f"strain_amr_genes ({len(strain_amr_genes)}) and labels ({len(labels)})"
        )
    gene_sets = [strain_amr_genes[sid] for sid in strain_ids]
    X, vocab = build_gene_presence_matrix(gene_sets, gene_vocabulary)
    y = np.array([labels[sid] for sid in strain_ids], dtype=int)
    clf = _train_baseline_xgboost(X, y, drug_name=drug_name)
    return clf, vocab, strain_ids


def train_kmer_baseline(
    strain_sequences: dict[str, "str | list[str]"],
    labels: dict[str, int],
    drug_name: str,
    k: int = DEFAULT_KMER_K,
    top_n: int = DEFAULT_KMER_TOP_N,
    classifier_type: str = "logreg",
    vocabulary: list[str] | None = None,
    contig_separator: str = CONTIG_SEPARATOR,
) -> tuple[TrainedClassifier, list[str], list[str]]:
    """Train k-mer baseline on per-strain sequences.

    Wave 3.5 M5 fix: accepts EITHER a single concatenated string OR a list of
    contigs per strain (chromosome + plasmids). When a list is supplied, the
    function joins contigs with a separator of N's (default 100 N's) so k-mers
    don't span contig boundaries — k-mer extraction skips N-containing windows.

    Without this, plasmid-borne β-lactamase signal was silently lost when
    callers passed only the chromosomal sequence.

    Args:
        strain_sequences: strain_id -> str OR list[str]. Either format works.
        labels: strain_id -> binary R/S label.
        drug_name: for reporting.
        k: k-mer size (default 8).
        top_n: vocabulary cap (default 10K most-frequent k-mers).
        classifier_type: 'logreg' (default) or 'xgboost'.
        vocabulary: pre-computed k-mer vocabulary (for held-out test sets).
        contig_separator: N-run inserted between contigs in list mode; long
            enough that no real k-mer (k≤32) can span it.

    Returns:
        (TrainedClassifier, kmer_vocabulary, strain_id_order).
    """
    strain_ids = sorted(strain_sequences.keys() & labels.keys())
    if not strain_ids:
        raise ClassifierTrainingError(
            f"k-mer baseline for {drug_name!r}: no strain overlap"
        )

    def _coerce(seq_or_list: "str | list[str]") -> str:
        if isinstance(seq_or_list, str):
            return seq_or_list
        # List of contigs → N-separated concatenation (k-mer N-skip handles boundaries)
        return contig_separator.join(seq_or_list)

    sequences = [_coerce(strain_sequences[sid]) for sid in strain_ids]
    if vocabulary is None:
        vocabulary = build_kmer_vocabulary(sequences, k=k, top_n=top_n)
    X = kmers_to_feature_matrix(sequences, vocabulary, k=k)
    y = np.array([labels[sid] for sid in strain_ids], dtype=int)

    if classifier_type == "logreg":
        clf = _train_baseline_logreg(X, y, drug_name=drug_name)
    elif classifier_type == "xgboost":
        clf = _train_baseline_xgboost(X, y, drug_name=drug_name)
    else:
        raise ValueError(
            f"Unknown classifier_type {classifier_type!r}; use 'logreg' or 'xgboost'"
        )
    return clf, vocabulary, strain_ids


def train_gene_presence_baseline(
    strain_gene_calls: dict[str, set[str]],
    labels: dict[str, int],
    drug_name: str,
    gene_vocabulary: list[str] | None = None,
) -> tuple[TrainedClassifier, list[str], list[str]]:
    """Train XGBoost on per-strain gene-family presence/absence.

    Same shape as AMRFinder baseline but consumes the broader Bakta-annotation
    gene-family set rather than only resistance genes. Lets the classifier
    discover non-AMR genes that correlate with resistance (typically efflux
    pumps, porins, regulators).

    Args:
        strain_gene_calls: strain_id -> set of gene-family symbols from Bakta.
        labels: strain_id -> binary R/S.
        drug_name: reporting.
        gene_vocabulary: optional pre-computed vocabulary.

    Returns:
        (TrainedClassifier, gene_vocabulary, strain_id_order).
    """
    # Reuses the AMRFinder baseline implementation — same feature shape, just
    # a different upstream catalog.
    return train_amrfinder_baseline(
        strain_amr_genes=strain_gene_calls,
        labels=labels,
        drug_name=drug_name,
        gene_vocabulary=gene_vocabulary,
    )
