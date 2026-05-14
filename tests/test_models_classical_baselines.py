"""Tests for Step 18 — Classical baseline benchmark."""
from __future__ import annotations

import numpy as np
import pytest

xgboost = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")

from dna_decode.models.classical_baselines import (  # noqa: E402
    CONTIG_SEPARATOR,
    DEFAULT_KMER_K,
    DEFAULT_KMER_TOP_N,
    _train_baseline_logreg,
    build_gene_presence_matrix,
    build_kmer_vocabulary,
    extract_kmer_counts,
    kmers_to_feature_matrix,
    train_amrfinder_baseline,
    train_gene_presence_baseline,
    train_kmer_baseline,
)
from dna_decode.models.classifiers import ClassifierTrainingError  # noqa: E402


# ---- k-mer feature extraction ----


def test_extract_kmer_counts_basic():
    counts = extract_kmer_counts("ATGCATGC", k=3)
    assert counts["ATG"] == 2
    assert counts["TGC"] == 2
    assert counts["GCA"] == 1
    assert counts["CAT"] == 1


def test_extract_kmer_counts_skips_n_containing():
    counts = extract_kmer_counts("ATNGAT", k=3)
    # ATN, TNG, NGA all have N → excluded. GAT remains.
    assert counts == {"GAT": 1}


def test_extract_kmer_counts_with_vocabulary_zero_fills_absent():
    counts = extract_kmer_counts(
        "ATGCATGC", k=3, vocabulary=["ATG", "TGC", "MISSING"]
    )
    assert counts == {"ATG": 2, "TGC": 2, "MISSING": 0}


def test_extract_kmer_counts_uppercase_normalization():
    """Lowercase sequence is uppercased so counts match canonical form."""
    counts = extract_kmer_counts("atgcatgc", k=3)
    assert counts["ATG"] == 2


def test_extract_kmer_counts_empty_sequence():
    assert extract_kmer_counts("", k=3) == {}


# ---- k-mer vocabulary building ----


def test_build_kmer_vocabulary_top_n_caps_size():
    # 6 unique 3-mers across two sequences; cap to top-3
    vocab = build_kmer_vocabulary(["ATGCATGC", "GGGGAAAA"], k=3, top_n=3)
    assert len(vocab) == 3


def test_build_kmer_vocabulary_returns_most_frequent():
    """ATGC repeated → ATG + TGC dominate."""
    vocab = build_kmer_vocabulary(["ATGCATGCATGC"], k=3, top_n=5)
    assert "ATG" in vocab[:3]
    assert "TGC" in vocab[:3]


def test_kmers_to_feature_matrix_shape_and_counts():
    seqs = ["ATGCATGC", "ATGCATGC"]
    vocab = ["ATG", "TGC", "MISSING"]
    X = kmers_to_feature_matrix(seqs, vocab, k=3)
    assert X.shape == (2, 3)
    assert X[0, 0] == 2  # ATG appears twice in ATGCATGC
    assert X[0, 2] == 0  # MISSING absent


# ---- gene presence/absence ----


def test_build_gene_presence_matrix_union_vocab():
    sets = [{"gyrA", "parC"}, {"gyrA", "tetA"}, {"parC"}]
    X, vocab = build_gene_presence_matrix(sets)
    assert sorted(vocab) == ["gyrA", "parC", "tetA"]
    assert X.shape == (3, 3)
    # Strain 0: gyrA + parC present, tetA absent
    g_idx = vocab.index("gyrA")
    p_idx = vocab.index("parC")
    t_idx = vocab.index("tetA")
    assert X[0, g_idx] == 1
    assert X[0, p_idx] == 1
    assert X[0, t_idx] == 0


def test_build_gene_presence_matrix_with_explicit_vocab():
    sets = [{"gyrA"}, {"parC"}]
    X, vocab = build_gene_presence_matrix(sets, gene_vocabulary=["gyrA", "parC", "extra"])
    assert vocab == ["gyrA", "parC", "extra"]
    assert X[0, 0] == 1 and X[0, 1] == 0 and X[0, 2] == 0
    assert X[1, 0] == 0 and X[1, 1] == 1 and X[1, 2] == 0


def test_build_gene_presence_matrix_empty():
    X, vocab = build_gene_presence_matrix([])
    assert X.shape == (0, 0)
    assert vocab == []


# ---- AMRFinder baseline ----


def _balanced_labels(strain_ids: list[str], pattern: str = "alt") -> dict[str, int]:
    """Half resistant, half susceptible (alternating)."""
    return {sid: i % 2 for i, sid in enumerate(strain_ids)}


def test_train_amrfinder_baseline_round_trip():
    strain_ids = [f"s{i:03d}" for i in range(20)]
    strain_amr_genes = {sid: {"gyrA"} if i % 2 else {"tetA"} for i, sid in enumerate(strain_ids)}
    labels = _balanced_labels(strain_ids)

    clf, vocab, returned_ids = train_amrfinder_baseline(strain_amr_genes, labels, "cipro")
    assert sorted(vocab) == ["gyrA", "tetA"]
    assert returned_ids == strain_ids
    assert clf.drug_name == "cipro"
    assert clf.feature_dim == 2


def test_train_amrfinder_baseline_no_overlap_raises():
    with pytest.raises(ClassifierTrainingError, match="no strain overlap"):
        train_amrfinder_baseline(
            strain_amr_genes={"a": set()}, labels={"b": 1}, drug_name="cipro"
        )


# ---- k-mer baseline ----


def test_train_kmer_baseline_default_logreg():
    rng = np.random.default_rng(0)
    bases = ["A", "C", "G", "T"]
    strain_sequences = {
        f"s{i:03d}": "".join(rng.choice(bases, size=200)) for i in range(30)
    }
    strain_ids = sorted(strain_sequences.keys())
    labels = _balanced_labels(strain_ids)

    clf, vocab, _ = train_kmer_baseline(
        strain_sequences, labels, "cipro", k=4, top_n=50
    )
    assert clf.drug_name == "cipro"
    assert len(vocab) <= 50
    assert clf.feature_dim == len(vocab)


def test_train_kmer_baseline_xgboost_classifier_type():
    rng = np.random.default_rng(1)
    bases = ["A", "C", "G", "T"]
    strain_sequences = {
        f"s{i:03d}": "".join(rng.choice(bases, size=200)) for i in range(30)
    }
    strain_ids = sorted(strain_sequences.keys())
    labels = _balanced_labels(strain_ids)

    clf, _, _ = train_kmer_baseline(
        strain_sequences, labels, "cipro", k=4, top_n=20, classifier_type="xgboost"
    )
    assert clf.calibrated is True


def test_train_kmer_baseline_unknown_classifier_raises():
    with pytest.raises(ValueError, match="classifier_type"):
        train_kmer_baseline(
            {"s001": "ATGC" * 50},
            {"s001": 1},
            "cipro",
            classifier_type="random_forest",
        )


# ---- Wave 3.5 M5: train_kmer_baseline accepts list[str] per strain ----


def test_train_kmer_baseline_accepts_list_of_contigs():
    """Multi-contig strain (chromosome + plasmids) joined internally with N-separator
    so k-mers don't span contig boundaries. Verifies the API broadening."""
    rng = np.random.default_rng(2)
    bases = ["A", "C", "G", "T"]
    # Mix str + list within the SAME dict — both should work
    strain_sequences: dict[str, "str | list[str]"] = {}
    for i in range(20):
        seqs = [
            "".join(rng.choice(bases, size=100)),
            "".join(rng.choice(bases, size=80)),
        ]
        strain_sequences[f"s_list_{i:03d}"] = seqs  # list of 2 contigs
    for i in range(20):
        strain_sequences[f"s_str_{i:03d}"] = "".join(rng.choice(bases, size=200))  # plain string

    labels = {sid: i % 2 for i, sid in enumerate(sorted(strain_sequences.keys()))}
    clf, vocab, _ = train_kmer_baseline(
        strain_sequences, labels, "cipro", k=4, top_n=30
    )
    # Should produce a trained classifier without raising; both input shapes accepted
    assert clf.drug_name == "cipro"
    assert clf.feature_dim == len(vocab)


def test_train_kmer_baseline_n_separator_prevents_cross_contig_kmers():
    """N-run between contigs should prevent k-mers that span the join."""
    from dna_decode.models.classical_baselines import extract_kmer_counts

    contigs = ["ATGCATGC", "GGGGAAAA"]
    joined = ("N" * 100).join(contigs)
    counts = extract_kmer_counts(joined, k=3)
    # k=3 kmers entirely inside each contig are present
    assert "ATG" in counts
    assert "GGG" in counts
    # No k-mer should contain N
    for kmer in counts:
        assert "N" not in kmer


# ---- gene-presence baseline (delegates to amrfinder) ----


def test_train_gene_presence_baseline_round_trip():
    strain_ids = [f"s{i:03d}" for i in range(20)]
    # 'omp' gene present in resistant, 'rib' in susceptible
    strain_gene_calls = {
        sid: ({"omp", "acrA"} if i % 2 else {"rib", "lacZ"})
        for i, sid in enumerate(strain_ids)
    }
    labels = _balanced_labels(strain_ids)

    clf, vocab, _ = train_gene_presence_baseline(strain_gene_calls, labels, "cipro")
    assert clf.feature_dim == len(vocab)
    assert set(vocab) == {"omp", "acrA", "rib", "lacZ"}


# ---- Stage 1 Refactor: _train_baseline_logreg(calibrate=False) branch ----


def test_train_baseline_logreg_calibrate_false_returns_raw_logreg():
    """`_train_baseline_logreg(calibrate=False)` returns the bare `LogisticRegression`.

    Pins the Stage1 contract (Step 2.3): the NT-logreg call site disables
    CalibratedClassifierCV wrapping so AUROC measures representation quality
    rather than calibration-wrapper behavior. Regression guard against someone
    reintroducing the wrapper or flipping the default.
    """
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.linear_model import LogisticRegression

    rng = np.random.default_rng(0)
    X = rng.standard_normal((20, 6)).astype(np.float32)
    y = np.array([0] * 10 + [1] * 10)

    clf = _train_baseline_logreg(X, y, drug_name="cipro", calibrate=False)
    assert clf.calibrated is False
    assert isinstance(clf.model, LogisticRegression)
    assert not isinstance(clf.model, CalibratedClassifierCV)
    assert clf.feature_dim == 6
    # Fitted LogisticRegression exposes coef_ attribute; pin that the fit ran
    assert hasattr(clf.model, "coef_")
    assert clf.model.coef_.shape == (1, 6)


def test_train_baseline_logreg_calibrate_true_wraps_calibrated_classifier_cv():
    """Default `calibrate=True` (balanced labels, minority>=2) wraps in `CalibratedClassifierCV`.

    Sanity counterpart to the calibrate=False test — confirms the new branch
    didn't accidentally short-circuit the default path.
    """
    from sklearn.calibration import CalibratedClassifierCV

    rng = np.random.default_rng(1)
    X = rng.standard_normal((20, 6)).astype(np.float32)
    y = np.array([0] * 10 + [1] * 10)  # minority count = 10, well above the cv threshold

    clf = _train_baseline_logreg(X, y, drug_name="cipro", calibrate=True)
    assert clf.calibrated is True
    assert isinstance(clf.model, CalibratedClassifierCV)


def test_contig_separator_is_100_n_run():
    """Module constant `CONTIG_SEPARATOR` is 100 N's.

    Stage1 Step 2.9: promoted to a module-level constant so loso_kmer + Stage 1
    runner share one source of truth. Pins both the symbol and the value to
    prevent a silent shortening that could let k-mers (k<=32) span contig joins.
    """
    assert CONTIG_SEPARATOR == "N" * 100
    assert len(CONTIG_SEPARATOR) == 100
    assert set(CONTIG_SEPARATOR) == {"N"}
