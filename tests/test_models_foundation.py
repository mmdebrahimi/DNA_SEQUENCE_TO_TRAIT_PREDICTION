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
