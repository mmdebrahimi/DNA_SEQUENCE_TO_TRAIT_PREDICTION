"""Tests for Step 7 — Foundation model wrappers.

Tests exercise the abstract base contract via MockFoundationModel. Real
Evo/DNABERT-2/NT/GENA-LM loaders are tested at first real-data run; here we
verify the dispatch + sliding-window + batching are correct.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dna_decode.models.foundation import (
    DNABERT2Model,
    EvoModel,
    FoundationModelError,
    GenaLMModel,
    MockFoundationModel,
    ModelMetadata,
    NucleotideTransformerModel,
    model_factory,
)


# ---- MockFoundationModel: interface contract ----


def test_mock_embed_returns_correct_shape():
    m = MockFoundationModel()
    out = m.embed("ATGC")
    assert out.shape == (1, m.embedding_dim)


def test_mock_embed_is_deterministic():
    m = MockFoundationModel()
    a = m.embed("ATGC")
    b = m.embed("ATGC")
    np.testing.assert_array_equal(a, b)


def test_mock_embed_differs_by_sequence():
    m = MockFoundationModel()
    a = m.embed("ATGC")
    b = m.embed("GGGG")
    assert not np.array_equal(a, b)


def test_mock_embed_empty_raises():
    m = MockFoundationModel()
    with pytest.raises(ValueError, match="empty"):
        m.embed("")


def test_mock_embed_batch_basic():
    m = MockFoundationModel()
    out = m.embed_batch(["ATGC", "GGGG", "CCAA"])
    assert out.shape == (3, m.embedding_dim)


def test_mock_embed_batch_empty():
    m = MockFoundationModel()
    out = m.embed_batch([])
    assert out.shape == (0, m.embedding_dim)


# ---- Sliding-window aggregation ----


def test_short_sequence_is_one_window():
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    windows = m._slide_windows("A" * 50)
    assert len(windows) == 1
    assert windows[0] == "A" * 50


def test_long_sequence_slides_with_half_context_stride():
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    # 250-bp sequence, max_context=100, stride=50
    # Windows at: 0-100, 50-150, 100-200, 150-250 (last covers 150-250 fully)
    seq = "A" * 250
    windows = m._slide_windows(seq)
    assert all(len(w) <= 100 for w in windows)
    # First window starts at 0
    assert windows[0] == seq[:100]
    # Last window ends at the sequence end
    assert windows[-1].endswith(seq[-50:])


def test_long_sequence_embed_returns_multi_window():
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    out = m.embed("A" * 250)
    assert out.shape[0] >= 2
    assert out.shape[1] == 8


def test_embed_batch_mean_pools_long_sequences():
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=16, max_context=100)
    )
    out = m.embed_batch(["A" * 250, "G" * 50])
    assert out.shape == (2, 16)


# ---- model_factory ----


def test_factory_mock_no_config_needed():
    m = model_factory("mock")
    assert isinstance(m, MockFoundationModel)


def test_factory_unknown_model_raises():
    with pytest.raises(FoundationModelError, match="Unknown foundation"):
        model_factory("does_not_exist")


def test_factory_evo_from_config(project_root: Path):
    """Evo wrapper is constructed; weights are NOT loaded eagerly."""
    m = model_factory("evo", config_path=project_root / "config" / "datasources.yaml")
    assert isinstance(m, EvoModel)
    assert m.metadata.huggingface_id == "togethercomputer/evo-1-131k-base"
    assert m.metadata.embedding_dim == 4096
    assert not m._loaded  # lazy loading — no weights touched


def test_factory_dnabert2_from_config(project_root: Path):
    m = model_factory("dnabert2", config_path=project_root / "config" / "datasources.yaml")
    assert isinstance(m, DNABERT2Model)
    assert m.metadata.embedding_dim == 768


def test_factory_nucleotide_transformer_from_config(project_root: Path):
    m = model_factory(
        "nucleotide_transformer",
        config_path=project_root / "config" / "datasources.yaml",
    )
    assert isinstance(m, NucleotideTransformerModel)


def test_factory_gena_lm_from_config(project_root: Path):
    m = model_factory("gena_lm", config_path=project_root / "config" / "datasources.yaml")
    assert isinstance(m, GenaLMModel)


def test_factory_missing_config_raises(tmp_path: Path):
    with pytest.raises(FoundationModelError, match="Config not found"):
        model_factory("evo", config_path=tmp_path / "missing.yaml")


def test_factory_device_override():
    m = model_factory("mock", device="cpu")
    assert m.device == "cpu"


# ---- ModelMetadata + lazy loading guarantee ----


def test_metadata_exposes_static_fields():
    meta = ModelMetadata(name="test", huggingface_id="hub/x", embedding_dim=42, max_context=512)
    m = MockFoundationModel(meta)
    assert m.name == "test"
    assert m.embedding_dim == 42
    assert m.max_context == 512


def test_mock_does_not_eagerly_load_weights():
    m = MockFoundationModel()
    # Constructor does not flip _loaded — it's set on first embed
    assert not m._loaded
    m.embed("ATGC")
    assert m._loaded


# ---- Phase 2.5 batching fix: numerical equivalence ----


def test_mock_embed_window_batch_matches_per_sequence():
    """Default _embed_window_batch impl must match per-sequence loop exactly.

    Regression guard for Phase 2.5 cache.populate batching refactor. The
    default fallback in FoundationModel.base just stacks per-sequence calls,
    so output must be bit-exact equal.
    """
    m = MockFoundationModel()
    m._ensure_loaded()  # mock has no real weights but flips state
    seqs = ["ATCG" * 25, "GGAA" * 30, "TTTT" * 20, "ACGT" * 10]
    batched = m._embed_window_batch(seqs)
    per_seq = np.stack([m._embed_window(s) for s in seqs])
    np.testing.assert_array_equal(batched, per_seq)


def test_embed_batch_fast_path_used_for_single_window_seqs():
    """When all sequences fit in single window, embed_batch must route
    through _embed_window_batch (the fast batched path), not the
    per-sequence loop.
    """
    m = MockFoundationModel()
    seqs = ["ATCG" * 25, "GGAA" * 30]  # both < max_context=512
    # Spy on _embed_window_batch
    calls = {"batch": 0, "single": 0}
    real_batch = m._embed_window_batch
    real_single = m._embed_window

    def spy_batch(s):
        calls["batch"] += 1
        return real_batch(s)

    def spy_single(s):
        calls["single"] += 1
        return real_single(s)

    m._embed_window_batch = spy_batch
    m._embed_window = spy_single
    m.embed_batch(seqs)
    assert calls["batch"] == 1
    # Per-sequence path may still call _embed_window inside the batch fallback,
    # but the spy shows the fast path was invoked exactly once.


def test_embed_window_batch_empty_returns_zero_rows():
    """Default _embed_window_batch must short-circuit on empty input.

    Guards the `if not sequences: return np.empty((0, D))` early return —
    the np.stack fallback would raise on an empty list.
    """
    m = MockFoundationModel()
    out = m._embed_window_batch([])
    assert out.shape == (0, m.embedding_dim)
    assert out.dtype == np.float32


def test_embed_batch_slow_path_dispatches_when_any_seq_exceeds_max_context():
    """When any sequence > max_context, embed_batch must take the slow path
    (per-sequence loop with mean-pool), NOT the batched fast path.

    Guards the `all(len(s) <= max_context for s in sequences)` branch in
    embed_batch — a regression that always took the fast path would silently
    truncate long sequences via the subclass tokenizer's `truncation=True`.
    """
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    calls = {"batch": 0}
    real_batch = m._embed_window_batch

    def spy_batch(s):
        calls["batch"] += 1
        return real_batch(s)

    m._embed_window_batch = spy_batch
    # One short (fits) + one long (needs windowing) → must NOT use fast path.
    m.embed_batch(["A" * 50, "G" * 250])
    assert calls["batch"] == 0


def test_embed_batch_boundary_seq_equal_to_max_context_uses_fast_path():
    """Sequences whose length exactly equals max_context fit in one window
    and should be routed through the batched fast path (the predicate is
    `len(s) <= max_context`, inclusive).
    """
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    calls = {"batch": 0}
    real_batch = m._embed_window_batch

    def spy_batch(s):
        calls["batch"] += 1
        return real_batch(s)

    m._embed_window_batch = spy_batch
    m.embed_batch(["A" * 100, "G" * 99])  # both <= 100
    assert calls["batch"] == 1


def test_embed_batch_slow_path_matches_per_sequence_mean_pool():
    """Slow-path output must equal per-sequence embed-then-mean-pool.

    Numerical-correctness guard for the windowed path: for a long sequence,
    embed_batch's slow path computes per-window embeddings then mean-pools
    over the window axis. This must equal calling embed() once and reducing.
    """
    m = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=16, max_context=100)
    )
    sequences = ["A" * 250, "G" * 50, "C" * 320]
    batched = m.embed_batch(sequences)
    expected = np.stack([m.embed(s).mean(axis=0) for s in sequences])
    np.testing.assert_array_equal(batched, expected)
    assert batched.shape == (3, 16)


def _cuda_available() -> bool:
    try:
        import torch

        return torch.cuda.is_available()
    except ImportError:
        return False


@pytest.mark.slow
@pytest.mark.skipif(not _cuda_available(), reason="GPU + storage required")
def test_nt_embed_window_batch_matches_per_sequence():
    """NT batched forward must produce numerically-equivalent embeddings to
    per-sequence calls.

    Loads NT v2 100M once, embeds 5 fixed sequences both ways, asserts
    np.allclose(rtol=1e-4). Guards the Phase 2.5 batching refactor against
    silent numerical drift that would invalidate caches.
    """
    from dna_decode.models.foundation import model_factory

    m = model_factory("nucleotide_transformer", device="cuda")
    m._ensure_loaded()  # subclass methods don't auto-load; embed() does
    sequences = [
        "ACGT" * 50,  # 200 bp
        "GGGAAACCCTTT" * 15,  # 180 bp
        "ATATCGCGATAT" * 10,  # 120 bp
        "TTTGGGCCCAAA" * 25,  # 300 bp
        "AAACCCGGGTTTAAACCCGGGTTT" * 10,  # 240 bp
    ]
    batched = m._embed_window_batch(sequences)
    per_seq = np.stack([m._embed_window(s) for s in sequences])

    assert batched.shape == per_seq.shape, (
        f"shape mismatch: batched={batched.shape}, per_seq={per_seq.shape}"
    )
    max_diff = np.abs(batched - per_seq).max()
    assert np.allclose(batched, per_seq, rtol=1e-4, atol=1e-5), (
        f"NT batched vs per-sequence outputs diverge beyond tolerance "
        f"(max abs diff: {max_diff})"
    )
