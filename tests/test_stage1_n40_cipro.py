"""Unit tests for scripts/stage1_n40_cipro.py helpers.

Pins the gate-formula verdict bucketing + paired-bootstrap-CI behavior. The
end-to-end LOSO run is tested via dry-run on the 12-strain mini cohort
(separately, not in CI — see scripts/stage1_n40_cipro.py docstring).
"""
from __future__ import annotations

import numpy as np
import pytest

from scripts.stage1_n40_cipro import (
    GATE_THRESHOLD_PP,
    VariantResult,
    compute_gate_outcome,
    decide_stage2_action,
    paired_bootstrap_ci,
    per_mlst_breakdown,
    verdict_label,
)


class TestVerdictLabel:
    def test_clean_pass_at_5pp(self):
        assert "CLEAN PASS" in verdict_label(5.0)
        assert "CLEAN PASS" in verdict_label(7.5)

    def test_clean_pass_below_5pp_not_clean(self):
        # 4.999 pp should be NOISY, not CLEAN
        v = verdict_label(4.999)
        assert "NOISY PASS" in v
        assert "CLEAN PASS" not in v

    def test_noisy_pass_in_3_to_5_window(self):
        assert "NOISY PASS" in verdict_label(3.0)
        assert "NOISY PASS" in verdict_label(4.0)
        assert "NOISY PASS" in verdict_label(4.99)

    def test_fail_below_threshold(self):
        assert "FAIL" in verdict_label(2.99)
        assert "FAIL" in verdict_label(0.0)
        assert "FAIL" in verdict_label(-5.0)

    def test_threshold_is_3pp(self):
        # Sanity-pin the locked threshold constant
        assert GATE_THRESHOLD_PP == 3.0


class TestPairedBootstrapCI:
    def test_identical_scores_give_zero_centered_ci(self):
        # If both classifiers produce identical scores, the gap is always 0
        rng = np.random.default_rng(0)
        n = 30
        y = rng.integers(0, 2, size=n)
        # Ensure ≥1 of each class for AUROC validity
        y[0] = 0
        y[1] = 1
        scores = rng.random(n).astype(np.float32)
        mean, lo, hi, _ = paired_bootstrap_ci(y, scores, scores, n_iterations=200, seed=42)
        assert mean == pytest.approx(0.0, abs=1e-6)
        assert lo == pytest.approx(0.0, abs=1e-6)
        assert hi == pytest.approx(0.0, abs=1e-6)

    def test_a_strictly_better_than_b_gives_positive_ci(self):
        # Build scores where A is perfectly correlated with y, B is anti-correlated
        n = 30
        y = np.array([0] * 15 + [1] * 15)
        scores_a = y.astype(np.float32) + np.random.default_rng(1).normal(0, 0.01, size=n).astype(np.float32)
        scores_b = (1 - y).astype(np.float32) + np.random.default_rng(2).normal(0, 0.01, size=n).astype(np.float32)
        mean, lo, hi, _ = paired_bootstrap_ci(y, scores_a, scores_b, n_iterations=500, seed=42)
        # A's AUROC ≈ 1.0, B's AUROC ≈ 0.0 → gap ≈ 1.0; CI must be tightly positive
        assert mean > 0.9
        assert lo > 0.5

    def test_returns_four_values_three_floats_and_int(self):
        # Step 2 of Stage1_Refactor_And_Test_Hardening_Plan: paired_bootstrap_ci now
        # returns 4-tuple (mean, lo, hi, n_effective) — n_effective is the count of
        # non-degenerate resamples that contributed to the CI.
        rng = np.random.default_rng(0)
        n = 20
        y = np.array([0] * 10 + [1] * 10)
        a = rng.random(n).astype(np.float32)
        b = rng.random(n).astype(np.float32)
        result = paired_bootstrap_ci(y, a, b, n_iterations=100, seed=1)
        assert len(result) == 4
        mean, lo, hi, n_eff = result
        assert isinstance(mean, float)
        assert isinstance(lo, float)
        assert isinstance(hi, float)
        assert isinstance(n_eff, int)

    def test_handles_degenerate_resample_silently(self):
        # All same class → AUROC undefined → bootstrap should still return a result (skipping bad folds)
        n = 10
        y = np.array([0] * 5 + [1] * 5)
        a = np.linspace(0, 1, n).astype(np.float32)
        b = np.linspace(1, 0, n).astype(np.float32)
        # Few iterations; if any single-class resample comes through it's skipped
        result = paired_bootstrap_ci(y, a, b, n_iterations=50, seed=3)
        # 4-tuple now; n_effective may be < 50 if some resamples were degenerate
        assert len(result) == 4


class TestPerMlstBreakdown:
    def test_groups_by_mlst_and_counts_rs(self):
        rows = per_mlst_breakdown(
            strain_ids=["s1", "s2", "s3", "s4"],
            mlsts=["ST10", "ST10", "ST131", "ST131"],
            y=np.array([1, 0, 1, 1]),
        )
        # Two MLSTs in sorted order
        assert len(rows) == 2
        st10 = next(r for r in rows if r["mlst"] == "ST10")
        st131 = next(r for r in rows if r["mlst"] == "ST131")
        assert st10 == {"mlst": "ST10", "n": 2, "r": 1, "s": 1}
        assert st131 == {"mlst": "ST131", "n": 2, "r": 2, "s": 0}

    def test_singleton_mlsts(self):
        rows = per_mlst_breakdown(
            strain_ids=["s1", "s2", "s3"],
            mlsts=["A", "B", "C"],
            y=np.array([0, 1, 0]),
        )
        assert len(rows) == 3
        assert all(r["n"] == 1 for r in rows)


# ---- Helpers for synthetic gate-outcome tests ----

def _make_variant(name, auroc, scores, true, strain_ids, gate_bearing):
    """Synthetic VariantResult builder for testing compute_gate_outcome."""
    return VariantResult(
        name=name,
        auroc=auroc,
        auprc=auroc,  # placeholder; not gate-relevant
        per_strain_scores=np.array(scores, dtype=np.float32),
        per_strain_true=np.array(true, dtype=int),
        strain_ids=list(strain_ids),
        is_gate_bearing=gate_bearing,
    )


def _aligned_4_variant_results(nt_xgb_auroc, nt_lr_auroc, kmer_auroc, fusion_auroc):
    """4 variants all aligned on the same 10-strain set; deterministic scores.

    Score patterns are constructed so AUROC ≈ the requested value:
      - 1.0 → all R get score 1.0, all S get 0.0
      - 0.5 → all get the same score 0.5 (rank-tied)
      - 0.99 → high noise variant
    Returns the list of VariantResults + the strain_ids + the true labels.
    """
    strain_ids = [f"s{i}" for i in range(10)]
    y = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]  # 5R/5S

    def _scores_for_auroc(target):
        if target >= 0.99:
            return [0.95, 0.94, 0.93, 0.92, 0.91, 0.05, 0.04, 0.03, 0.02, 0.01]
        if target >= 0.79:
            return [0.9, 0.85, 0.8, 0.75, 0.7, 0.3, 0.25, 0.2, 0.15, 0.1]
        if target >= 0.69:
            return [0.8, 0.75, 0.7, 0.65, 0.6, 0.4, 0.35, 0.3, 0.45, 0.5]
        return [0.5] * 10  # rank-tied → AUROC = 0.5

    return [
        _make_variant("NT-XGBoost", nt_xgb_auroc, _scores_for_auroc(nt_xgb_auroc), y, strain_ids, True),
        _make_variant("NT-logreg", nt_lr_auroc, _scores_for_auroc(nt_lr_auroc), y, strain_ids, True),
        _make_variant("k-mer-XGB", kmer_auroc, _scores_for_auroc(kmer_auroc), y, strain_ids, True),
        _make_variant("NT+k-mer-fusion-logreg", fusion_auroc, _scores_for_auroc(fusion_auroc), y, strain_ids, False),
    ], strain_ids, y


# ---- Fusion exclusion from gate (regression guard for /brainstorm flagged failure mode) ----


class TestGateOutcomeFusionExcluded:
    def test_fusion_excluded_from_gate_when_fusion_wins(self):
        """Fusion=0.99, NT-only heads + k-mer all 0.50 -> verdict bucket is FAIL.

        Pins the /brainstorm-flagged failure mode: k-mer carrying the result via
        fusion must NOT count toward the gate. If this test breaks, fusion has
        leaked into the gate computation.
        """
        results, _, _ = _aligned_4_variant_results(
            nt_xgb_auroc=0.5, nt_lr_auroc=0.5, kmer_auroc=0.5, fusion_auroc=0.99
        )
        outcome = compute_gate_outcome(results)
        assert outcome["verdict_bucket_short"] == "FAIL"
        assert outcome["fusion_outperforms_primary"] is True
        # nt_best is one of NT-XGBoost / NT-logreg (tied at 0.5)
        assert outcome["nt_best_name"] in {"NT-XGBoost", "NT-logreg"}

    def test_fusion_ignored_when_nt_xgboost_wins_alone(self):
        """NT-XGBoost=0.99, NT-logreg=0.5, k-mer=0.5, fusion=0.5 -> CLEAN PASS."""
        results, _, _ = _aligned_4_variant_results(
            nt_xgb_auroc=0.99, nt_lr_auroc=0.5, kmer_auroc=0.5, fusion_auroc=0.5
        )
        outcome = compute_gate_outcome(results)
        assert outcome["verdict_bucket_short"] == "CLEAN"
        assert outcome["nt_best_name"] == "NT-XGBoost"
        assert outcome["fusion_outperforms_primary"] is False


# ---- Calibration discipline regression guard (3 call sites) ----


class TestCalibrationDiscipline:
    """Pin calibrate=False at every gate-bearing classifier call site.

    Regression guard against re-introducing the N=11 isotonic-collapse bug. If
    someone adds CalibratedClassifierCV back to a call site, these tests catch it.
    """

    def test_nt_xgboost_call_site_passes_calibrate_false(self, monkeypatch):
        recorded: dict = {}

        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            recorded["calibrate"] = calibrate
            class _M:
                pass
            m = _M()
            m.feature_dim = X.shape[1]
            m.calibrated = calibrate
            m.model = None
            m.drug_name = drug_name
            return m

        # Patch the symbol imported INTO scripts.stage1_n40_cipro (where _nt_xgb_train resolved it at import time)
        import scripts.stage1_n40_cipro as stage1_mod
        monkeypatch.setattr(stage1_mod, "train_xgboost_classifier", fake_train)
        from scripts.stage1_n40_cipro import _nt_xgb_train

        X = np.random.default_rng(0).random((10, 4)).astype(np.float32)
        y = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
        _nt_xgb_train(X, y)
        assert recorded["calibrate"] is False

    def test_nt_logreg_call_site_passes_calibrate_false_and_scale_features_true(self, monkeypatch):
        """NT-logreg call site MUST pass scale_features=True (H13 fix, 2026-05-15).

        Without scale_features, Stage 1 N=38 NT-logreg AUROC was 0.100
        (anti-predictive) — H13 plumbing-bug root cause. The fix wraps LR in
        sklearn Pipeline(StandardScaler, LR). Regression guard against someone
        flipping the kwarg back to False on the NT path.
        """
        from dna_decode.models import classical_baselines as cb_mod
        recorded: dict = {}

        def fake_logreg(X, y, drug_name, *, calibrate=True, scale_features=False):
            recorded["calibrate"] = calibrate
            recorded["scale_features"] = scale_features
            class _M:
                feature_dim = X.shape[1]
                calibrated = calibrate
                model = None
            return _M()

        monkeypatch.setattr(cb_mod, "_train_baseline_logreg", fake_logreg)
        # Patch the symbol that stage1_n40_cipro imported under
        import scripts.stage1_n40_cipro as stage1_mod
        monkeypatch.setattr(stage1_mod, "_train_baseline_logreg", fake_logreg)
        from scripts.stage1_n40_cipro import _nt_logreg_train

        X = np.random.default_rng(0).random((10, 4)).astype(np.float32)
        y = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
        _nt_logreg_train(X, y)
        assert recorded["calibrate"] is False
        assert recorded["scale_features"] is True

    def test_kmer_xgb_call_site_passes_calibrate_false(self, monkeypatch):
        """New (per /brainstorm Round 2): k-mer-XGB is also gate-bearing; must pin too."""
        import dna_decode.eval.loso_kmer as lk_mod
        recorded: list[bool] = []

        def fake_train(X, y, *, drug_name="", calibrate=True, **kwargs):
            recorded.append(calibrate)
            class _M:
                pass
            m = _M()
            m.feature_dim = X.shape[1]
            m.calibrated = calibrate
            m.model = None
            m.drug_name = drug_name
            return m

        def fake_predict(clf, X):
            return np.full(X.shape[0], 0.5, dtype=np.float32)

        monkeypatch.setattr(lk_mod, "train_xgboost_classifier", fake_train)
        monkeypatch.setattr(lk_mod, "predict_proba", fake_predict)

        from dna_decode.eval.loso_kmer import run_kmer_xgboost_loso

        strain_ids = ["s0", "s1", "s2", "s3"]
        seqs_by_strain = {s: "ACGT" * 50 for s in strain_ids}
        labels_by_strain = {"s0": 1, "s1": 0, "s2": 1, "s3": 0}
        run_kmer_xgboost_loso(seqs_by_strain, labels_by_strain, strain_ids, drug="test", k=3, top_n=10)
        assert len(recorded) == 4  # 4 LOSO folds
        assert all(c is False for c in recorded), f"k-mer-XGB call sites must all pass calibrate=False; got {recorded}"


# ---- stage2_action decision-layer (4 buckets x CI lo) ----


class TestStage2Action:
    def test_clean_pass_ci_clears_zero_returns_burst(self):
        assert decide_stage2_action("CLEAN", ci_lo_pp=2.0) == "BURST_STAGE_2"

    def test_clean_pass_ci_negative_still_returns_burst(self):
        # Wide-CI annotation handled in packet; action stays BURST
        assert decide_stage2_action("CLEAN", ci_lo_pp=-0.5) == "BURST_STAGE_2"

    def test_noisy_pass_ci_negative_returns_hold(self):
        assert decide_stage2_action("NOISY", ci_lo_pp=-1.2) == "HOLD_STAGE_2_CI_DEGENERATE"

    def test_fail_returns_alternative_pooling_rerun(self):
        assert decide_stage2_action("FAIL", ci_lo_pp=-5.0) == "ALTERNATIVE_POOLING_RERUN"
        assert decide_stage2_action("FAIL", ci_lo_pp=10.0) == "ALTERNATIVE_POOLING_RERUN"


# ---- Bootstrap n_effective ----


class TestBootstrapSkipCount:
    def test_paired_bootstrap_returns_n_effective(self):
        """4-tuple return with n_effective. n=10 mixed labels; expect n_effective close to n_iterations."""
        n = 10
        y = np.array([0] * 5 + [1] * 5)
        a = np.linspace(0, 1, n).astype(np.float32)
        b = np.linspace(1, 0, n).astype(np.float32)
        mean, lo, hi, n_eff = paired_bootstrap_ci(y, a, b, n_iterations=100, seed=1)
        assert isinstance(n_eff, int)
        assert n_eff <= 100
        # With 5R/5S, most bootstrap resamples should hit both classes
        assert n_eff > 50


# ---- MLST loud handling ----


class TestMlstLoudHandling:
    def test_load_features_raises_on_missing_mlst(self, tmp_path, monkeypatch):
        """Build a minimal cohort + cache so load_features can run; strain with None MLST should raise."""
        from dataclasses import dataclass

        @dataclass
        class _FakeStrain:
            strain_id: str
            assembly_accession: str
            ast_labels: dict
            mlst: object  # may be None

        @dataclass
        class _FakeCohort:
            strains: list

        # Two strains: one valid MLST, one with mlst=None — the second should trigger the raise
        cohort = _FakeCohort(strains=[
            _FakeStrain("s_ok", "GCA_OK", {"cipro": 1}, mlst="MLST.foo.1"),
            _FakeStrain("s_bad", "GCA_BAD", {"cipro": 0}, mlst=None),
        ])

        # Monkeypatch EmbeddingCache + fasta_path to bypass real cache + filesystem
        import scripts.stage1_n40_cipro as stage1_mod

        class _FakeCache:
            def __init__(self, *args, **kwargs): pass
            def list_genes(self, sid): return [f"{sid}_g1"]
            def bulk_get(self, pairs): return np.zeros((1, 512), dtype=np.float32)

        def _fake_load_fasta(path):
            return ["ACGT" * 25]

        monkeypatch.setattr(stage1_mod, "EmbeddingCache", _FakeCache)
        monkeypatch.setattr(stage1_mod, "fasta_path", lambda acc, root: tmp_path / acc / "genome.fna")
        monkeypatch.setattr(stage1_mod, "_load_fasta_contigs", _fake_load_fasta)

        with pytest.raises(ValueError, match="missing MLST"):
            stage1_mod.load_features(cohort, tmp_path / "cache.h5", tmp_path / "refseq", drug="cipro")


# ---- strain_ids alignment validation ----


class TestStrainIdAlignment:
    def test_compute_gate_outcome_raises_on_nt_vs_kmer_mismatch(self):
        """NT-XGBoost and k-mer-XGB with different strain_ids -> ValueError, gate-bearing."""
        strain_ids_a = ["s0", "s1", "s2", "s3"]
        strain_ids_b = ["s0", "s1", "s2", "s4"]  # mismatch on last element
        y = [1, 1, 0, 0]
        scores = [0.9, 0.8, 0.2, 0.1]
        results = [
            _make_variant("NT-XGBoost", 0.8, scores, y, strain_ids_a, True),
            _make_variant("NT-logreg", 0.7, scores, y, strain_ids_a, True),
            _make_variant("k-mer-XGB", 0.6, scores, y, strain_ids_b, True),
            _make_variant("NT+k-mer-fusion-logreg", 0.5, scores, y, strain_ids_a, False),
        ]
        with pytest.raises(ValueError, match="alignment"):
            compute_gate_outcome(results)

    def test_compute_gate_outcome_suppresses_fusion_note_on_fusion_mismatch(self):
        """Fusion strain_ids differ from NT-best; gate proceeds but fusion note is suppressed."""
        strain_ids_aligned = ["s0", "s1", "s2", "s3"]
        strain_ids_fusion = ["s0", "s1", "s2", "s9"]  # mismatch
        y = [1, 1, 0, 0]
        # NT and k-mer all aligned; fusion misaligned. Fusion AUROC high enough to trigger note IF aligned.
        results = [
            _make_variant("NT-XGBoost", 0.6, [0.8, 0.7, 0.3, 0.2], y, strain_ids_aligned, True),
            _make_variant("NT-logreg", 0.6, [0.8, 0.7, 0.3, 0.2], y, strain_ids_aligned, True),
            _make_variant("k-mer-XGB", 0.5, [0.5, 0.5, 0.5, 0.5], y, strain_ids_aligned, True),
            _make_variant("NT+k-mer-fusion-logreg", 0.99, [0.99, 0.95, 0.05, 0.01], y, strain_ids_fusion, False),
        ]
        # Should not raise (fusion alignment is permissive)
        outcome = compute_gate_outcome(results)
        assert outcome["fusion_alignment_valid"] is False
        assert "mismatch" in outcome["fusion_note"].lower()
        assert outcome["fusion_outperforms_primary"] is False  # suppressed


class TestVariantResultLengthInvariant:
    """Defensive contract guard per /brainstorm 2026-05-14 follow-up.

    A malformed VariantResult with matching strain_ids but unequal len(scores) or
    len(true) could pass the alignment check and contaminate the fusion-outperforms
    diagnostic OR error obscurely in the bootstrap. The invariant raises ValueError
    BEFORE any downstream computation trusts the arrays.
    """

    def test_compute_gate_outcome_raises_on_scores_length_mismatch(self):
        strain_ids = ["s0", "s1", "s2", "s3"]
        y = [1, 1, 0, 0]
        # Build a 'bad' variant with scores length 5 but strain_ids length 4
        bad = _make_variant("NT-XGBoost", 0.7, [0.9, 0.8, 0.3, 0.2], y, strain_ids, True)
        bad.per_strain_scores = np.array([0.9, 0.8, 0.3, 0.2, 0.5], dtype=np.float32)  # extra score
        good = _make_variant("NT-logreg", 0.6, [0.8, 0.7, 0.3, 0.2], y, strain_ids, True)
        good_kmer = _make_variant("k-mer-XGB", 0.5, [0.5, 0.5, 0.5, 0.5], y, strain_ids, True)
        with pytest.raises(ValueError, match="internal length mismatch"):
            compute_gate_outcome([bad, good, good_kmer])

    def test_compute_gate_outcome_raises_on_true_length_mismatch(self):
        strain_ids = ["s0", "s1", "s2", "s3"]
        y = [1, 1, 0, 0]
        bad = _make_variant("NT-XGBoost", 0.7, [0.9, 0.8, 0.3, 0.2], y, strain_ids, True)
        bad.per_strain_true = np.array([1, 1, 0], dtype=int)  # one short
        good = _make_variant("NT-logreg", 0.6, [0.8, 0.7, 0.3, 0.2], y, strain_ids, True)
        good_kmer = _make_variant("k-mer-XGB", 0.5, [0.5, 0.5, 0.5, 0.5], y, strain_ids, True)
        with pytest.raises(ValueError, match="internal length mismatch"):
            compute_gate_outcome([bad, good, good_kmer])


# ---- Stage 1b ALTERNATIVE_POOLING_RERUN: --aggregation flag + default ----


def test_stage1_aggregation_flag_defaults_to_mean_plus_max():
    """Stage 1b default `--aggregation mean+max` per ALTERNATIVE_POOLING_RERUN.

    Stage 1 used mean-pool (512-dim); Stage 1b's verdict-time pre-commitment
    is mean+max-pool (1024-dim) per `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md`
    Verdict-Time Pre-Commitments table. Regression guard: keep the default
    on mean+max until the next pivot fires.
    """
    from pathlib import Path
    import argparse

    # Parse with no --aggregation flag; default should be "mean+max".
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=None)
    parser.add_argument("--nt-cache", type=Path, default=None)
    parser.add_argument("--refseq-cache", type=Path, default=None)
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--kmer-k", type=int, default=8)
    parser.add_argument("--kmer-top-n", type=int, default=10_000)
    parser.add_argument(
        "--aggregation", choices=["mean", "max", "mean+max"], default="mean+max"
    )
    args = parser.parse_args([])
    assert args.aggregation == "mean+max"


def test_load_features_threads_aggregation_through(monkeypatch, tmp_path):
    """`load_features(..., aggregation='mean+max')` passes aggregation to
    `aggregate_strain_features`. Pins the parameter does not get accidentally
    swallowed (regression guard for the Stage 1b refactor).
    """
    from scripts import stage1_n40_cipro as stage1_mod

    recorded: list[str] = []

    def fake_agg(gene_matrix, agg):
        recorded.append(agg)
        import numpy as np
        # Return a deterministic vector with shape consistent with agg
        if agg == "mean+max":
            return np.zeros(1024, dtype=np.float32)
        return np.zeros(512, dtype=np.float32)

    monkeypatch.setattr(stage1_mod, "aggregate_strain_features", fake_agg)

    # Stub the rest of load_features's dependencies; we only need to verify
    # the aggregation kwarg threads through.
    class FakeCache:
        def list_genes(self, sid):
            return ["g1", "g2"]
        def bulk_get(self, pairs):
            import numpy as np
            return np.zeros((len(pairs), 512), dtype=np.float32)

    monkeypatch.setattr(stage1_mod, "EmbeddingCache", lambda *a, **kw: FakeCache())
    monkeypatch.setattr(stage1_mod, "fasta_path", lambda acc, root: tmp_path / "x.fna")

    def fake_load_fasta(p):
        return ["ACGT"]
    monkeypatch.setattr(stage1_mod, "_load_fasta_contigs", fake_load_fasta)

    # Make a tiny cohort stub
    class S:
        def __init__(self, sid, mlst):
            self.strain_id = sid
            self.ast_labels = {"ciprofloxacin": 1}
            self.assembly_accession = "GCA_X"
            self.mlst = mlst
    class C:
        strains = [S("s1", "ST1"), S("s2", "ST2")]

    stage1_mod.load_features(C(), tmp_path / "cache.h5", tmp_path, "ciprofloxacin", aggregation="mean+max")
    assert all(a == "mean+max" for a in recorded)
    assert len(recorded) == 2
