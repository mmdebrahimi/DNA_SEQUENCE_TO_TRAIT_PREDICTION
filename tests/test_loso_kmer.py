"""Tests for `dna_decode/eval/loso_kmer.py` — shared LOSO runners for k-mer + fusion.

Per Stage1_Refactor_And_Test_Hardening_Plan Step 2.1: order-explicit LOSO runners
used by both the 12-strain smoke gate and the N=40 Stage 1 runner. These tests pin:

  - `CVResult.strain_ids` aligns with the caller-supplied `strain_ids` (no
    internal sort, no dedup, no reorder).
  - Within-fold vocabulary rebuild: held-out strain's k-mers MUST NOT influence
    the training-fold vocab. This is the leakage guard the factored module exists
    to enforce.
  - Error passthrough: `ClassifierTrainingError` is re-raised, not silently
    converted to a mean fallback (the per-Step-2.4 contract).
  - Fusion `X_nt` shape mismatch raises `ValueError`.
  - Fusion runner returns a `CVResult` with one fold per strain in input order.

The k-mer-XGB `calibrate=False` call-site pin lives in `tests/test_stage1_n40_cipro.py`
(`TestCalibrationDiscipline.test_kmer_xgb_call_site_passes_calibrate_false`) — not
duplicated here.
"""
from __future__ import annotations

import numpy as np
import pytest

xgboost = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")

from dna_decode.eval.cv import CVResult  # noqa: E402
from dna_decode.eval.loso_kmer import run_fusion_loso, run_kmer_xgboost_loso  # noqa: E402
from dna_decode.models.classifiers import ClassifierTrainingError  # noqa: E402


# ---- Fixtures ----


def _make_strains(n: int = 6, seed: int = 0) -> tuple[dict[str, str], dict[str, int], list[str]]:
    """n strains with deterministic random sequences + alternating R/S labels.

    Returns (seqs_by_strain, labels_by_strain, strain_ids).
    """
    rng = np.random.default_rng(seed)
    bases = ["A", "C", "G", "T"]
    strain_ids = [f"s{i:03d}" for i in range(n)]
    seqs_by_strain = {sid: "".join(rng.choice(bases, size=300)) for sid in strain_ids}
    labels_by_strain = {sid: i % 2 for i, sid in enumerate(strain_ids)}
    return seqs_by_strain, labels_by_strain, strain_ids


# ---- run_kmer_xgboost_loso: shape + alignment contract ----


class TestKmerXgboostLoSoShape:
    def test_returns_cvresult_with_one_fold_per_strain(self, monkeypatch):
        # Stub the heavy classifier path — we're testing the runner's plumbing.
        import dna_decode.eval.loso_kmer as lk_mod

        class _StubClf:
            feature_dim = 0
            calibrated = False
            model = None
            drug_name = ""

        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            stub = _StubClf()
            stub.feature_dim = X.shape[1]
            stub.drug_name = drug_name
            stub.calibrated = calibrate
            return stub

        def fake_predict(_clf, X):
            return np.full(X.shape[0], 0.5, dtype=np.float32)

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", fake_train)
        monkeypatch.setattr(lk_mod, "predict_proba", fake_predict)

        seqs, labels, ids = _make_strains(n=6)
        cv = run_kmer_xgboost_loso(seqs, labels, ids, drug="cipro_test", k=4, top_n=20)
        assert isinstance(cv, CVResult)
        assert cv.strategy == "loso"
        assert cv.drug == "cipro_test"
        assert cv.n_folds == 6
        for fold in cv.folds:
            assert fold.n_test == 1
            assert fold.n_train == 5
            assert fold.y_true.shape == (1,)
            assert fold.y_score.shape == (1,)

    def test_cvresult_strain_ids_matches_input_strain_ids_exactly(self, monkeypatch):
        """`CVResult.strain_ids` returns fold order = input `strain_ids` verbatim.

        Locks the alignment contract that the factored loso_kmer module exists to
        enforce (post-/brainstorm Round 1 alignment-drift critique).
        """
        import dna_decode.eval.loso_kmer as lk_mod

        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            class _M:
                pass
            m = _M()
            m.feature_dim = X.shape[1]
            m.drug_name = drug_name
            m.calibrated = calibrate
            return m

        def fake_predict(_clf, X):
            return np.full(X.shape[0], 0.5, dtype=np.float32)

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", fake_train)
        monkeypatch.setattr(lk_mod, "predict_proba", fake_predict)

        seqs, labels, _ = _make_strains(n=5)
        # Use a deliberately non-alphabetic order — the runner must NOT re-sort
        custom_order = ["s002", "s000", "s004", "s001", "s003"]
        cv = run_kmer_xgboost_loso(seqs, labels, custom_order, drug="d", k=4, top_n=10)
        assert cv.strain_ids == custom_order

    def test_fold_held_out_indices_match_input_position(self, monkeypatch):
        """`fold.held_out_indices` is `[i]` where i = position in input strain_ids.

        Pins the row-position contract that downstream paired comparisons assume.
        """
        import dna_decode.eval.loso_kmer as lk_mod

        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            class _M:
                pass
            m = _M()
            m.feature_dim = X.shape[1]
            m.drug_name = drug_name
            m.calibrated = calibrate
            return m

        def fake_predict(_clf, X):
            return np.full(X.shape[0], 0.5, dtype=np.float32)

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", fake_train)
        monkeypatch.setattr(lk_mod, "predict_proba", fake_predict)

        seqs, labels, ids = _make_strains(n=4)
        cv = run_kmer_xgboost_loso(seqs, labels, ids, drug="d", k=4, top_n=10)
        for i, fold in enumerate(cv.folds):
            assert fold.held_out_indices == [i]
            assert fold.held_out_id == ids[i]
            # train_indices = all positions except i
            assert sorted(fold.train_indices) == [j for j in range(4) if j != i]


# ---- run_kmer_xgboost_loso: within-fold vocab rebuild correctness ----


class TestKmerXgboostLoSoVocabLeak:
    def test_holdout_only_kmer_does_not_appear_in_train_vocab(self, monkeypatch):
        """If a k-mer appears ONLY in the held-out strain, it must not be in the
        per-fold training vocabulary.

        Constructs a synthetic case where strain `s_unique` contains a distinctive
        repeated motif `"GGGGGGGG"` that no other strain has. When LOSO holds out
        `s_unique`, the rebuilt training vocab must exclude that motif.

        Intercepts the per-fold cached vocab builder inside the runner to record
        which vocab was actually used per fold.
        """
        import dna_decode.eval.loso_kmer as lk_mod

        # All other strains: cycled ACGT, no GGGGGGGG.
        # s_unique: starts with 50 G's so "GGGGGGGG" (k=8) shows up many times.
        other_template = "ACGT" * 50  # 200 bases, no run-of-Gs
        seqs_by_strain = {
            "s_other_0": other_template,
            "s_other_1": other_template,
            "s_other_2": other_template,
            "s_other_3": other_template,
            "s_unique": ("G" * 50) + other_template,
        }
        labels_by_strain = {sid: i % 2 for i, sid in enumerate(sorted(seqs_by_strain))}
        strain_ids = ["s_other_0", "s_other_1", "s_other_2", "s_other_3", "s_unique"]

        captured_vocabs: list[tuple[str, list[str]]] = []

        # Hook the per-fold cached vocab builder (2026-06-05 perf refactor: the runner builds
        # the per-fold vocab via _vocab_from_cache(train_strains, counts_by_strain, top_n) —
        # train-only, fresh per fold — instead of build_kmer_vocabulary(train_seqs). The
        # leakage property is unchanged; we spy the new seam.
        real_vocab_from_cache = lk_mod._vocab_from_cache

        def spy_build(train_strains, counts_by_strain, top_n=10_000):
            vocab = real_vocab_from_cache(train_strains, counts_by_strain, top_n)
            captured_vocabs.append(("vocab", list(vocab)))
            return vocab

        monkeypatch.setattr(lk_mod, "_vocab_from_cache", spy_build)

        # Stub the classifier so we don't pay XGBoost cost
        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            class _M:
                pass
            m = _M()
            m.feature_dim = X.shape[1]
            m.drug_name = drug_name
            m.calibrated = calibrate
            return m

        def fake_predict(_clf, X):
            return np.full(X.shape[0], 0.5, dtype=np.float32)

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", fake_train)
        monkeypatch.setattr(lk_mod, "predict_proba", fake_predict)

        cv = run_kmer_xgboost_loso(
            seqs_by_strain, labels_by_strain, strain_ids, drug="vocab_test", k=8, top_n=50
        )
        assert cv.n_folds == 5
        assert len(captured_vocabs) == 5

        # Fold index 4 holds out s_unique. The vocab for that fold should NOT
        # contain "GGGGGGGG" — the held-out strain's signature 8-mer must be
        # excluded from the training-fold vocab.
        _, vocab_when_unique_held_out = captured_vocabs[4]
        assert "GGGGGGGG" not in vocab_when_unique_held_out, (
            "Held-out strain's signature k-mer leaked into training-fold vocab — "
            "within-fold rebuild is broken"
        )

        # Sanity counterpart: when s_unique is IN the training set (fold 0
        # holds out s_other_0), GGGGGGGG MUST appear in that fold's vocab.
        _, vocab_when_unique_in_train = captured_vocabs[0]
        assert "GGGGGGGG" in vocab_when_unique_in_train, (
            "Sanity check failed: GGGGGGGG should appear in train vocab when s_unique is in training set"
        )

    def test_vocab_rebuilt_per_fold_not_reused(self, monkeypatch):
        """The runner builds a fresh vocabulary per fold — not once outside the loop.

        Pins that the per-fold vocab builder is invoked exactly `n` times (one per
        LOSO fold), not once globally. (2026-06-05: the per-fold seam is now
        _vocab_from_cache; the per-genome COUNT cache is built once, but the
        train-only vocab is still rebuilt fresh each fold.)
        """
        import dna_decode.eval.loso_kmer as lk_mod

        real_vocab_from_cache = lk_mod._vocab_from_cache
        call_count = {"n": 0}

        def counting_build(train_strains, counts_by_strain, top_n=10_000):
            call_count["n"] += 1
            return real_vocab_from_cache(train_strains, counts_by_strain, top_n)

        monkeypatch.setattr(lk_mod, "_vocab_from_cache", counting_build)

        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            class _M:
                pass
            m = _M()
            m.feature_dim = X.shape[1]
            m.drug_name = drug_name
            m.calibrated = calibrate
            return m

        def fake_predict(_clf, X):
            return np.full(X.shape[0], 0.5, dtype=np.float32)

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", fake_train)
        monkeypatch.setattr(lk_mod, "predict_proba", fake_predict)

        seqs, labels, ids = _make_strains(n=7)
        run_kmer_xgboost_loso(seqs, labels, ids, drug="d", k=4, top_n=20)
        assert call_count["n"] == 7  # one vocab rebuild per fold


# ---- run_kmer_xgboost_loso: error passthrough ----


class TestKmerXgboostLoSoErrorPassthrough:
    def test_classifier_training_error_propagates_not_swallowed(self, monkeypatch):
        """If `train_xgboost_classifier` raises `ClassifierTrainingError`, the
        runner re-raises (NO silent mean fallback) — per Step 2.4."""
        import dna_decode.eval.loso_kmer as lk_mod

        def failing_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            raise ClassifierTrainingError("synthetic-failure-for-test")

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", failing_train)

        seqs, labels, ids = _make_strains(n=5)
        with pytest.raises(ClassifierTrainingError, match="synthetic-failure-for-test"):
            run_kmer_xgboost_loso(seqs, labels, ids, drug="d", k=4, top_n=10)


# ---- run_fusion_loso: shape + alignment ----


class TestFusionLoSoShape:
    def test_x_nt_row_count_mismatch_raises(self):
        """`X_nt.shape[0] != len(strain_ids)` -> ValueError before any training."""
        seqs, labels, ids = _make_strains(n=5)
        # X_nt has wrong number of rows (4 vs len(ids)=5)
        X_nt_bad = np.zeros((4, 16), dtype=np.float32)
        with pytest.raises(ValueError, match="X_nt.shape"):
            run_fusion_loso(X_nt_bad, seqs, labels, ids, drug="d", k=4, top_n=10)

    def test_returns_cvresult_one_fold_per_strain_in_input_order(self):
        """Fusion runs end-to-end with real LogisticRegression (cheap path).

        Pins both the per-fold shape and the strain_ids alignment with input.
        """
        seqs, labels, ids = _make_strains(n=6)
        # Need at least one R and one S per training fold; the alternating labels
        # from _make_strains guarantee this for n=6.
        X_nt = np.random.default_rng(7).standard_normal((6, 8)).astype(np.float32)
        cv = run_fusion_loso(X_nt, seqs, labels, ids, drug="fusion_test", k=4, top_n=20)
        assert cv.strategy == "loso"
        assert cv.drug == "fusion_test"
        assert cv.n_folds == 6
        assert cv.strain_ids == ids
        for i, fold in enumerate(cv.folds):
            assert fold.held_out_id == ids[i]
            assert fold.y_true.shape == (1,)
            assert fold.y_score.shape == (1,)
            # Score is a probability ∈ [0, 1]
            assert 0.0 <= float(fold.y_score[0]) <= 1.0

    def test_fusion_uses_caller_strain_order_not_sorted(self):
        """Reorder `strain_ids` and pin that `cv.strain_ids` matches input order."""
        seqs, labels, _ = _make_strains(n=6)
        custom_order = ["s003", "s001", "s005", "s000", "s004", "s002"]
        # X_nt rows aligned to custom_order (caller contract per docstring)
        rng = np.random.default_rng(11)
        X_nt = rng.standard_normal((6, 8)).astype(np.float32)
        cv = run_fusion_loso(X_nt, seqs, labels, custom_order, drug="d", k=4, top_n=20)
        assert cv.strain_ids == custom_order
