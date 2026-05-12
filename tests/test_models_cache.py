"""Tests for Step 8 — HDF5 embedding cache."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# Skip the entire module if h5py is not installed (lighter test setups). The
# real test run via `uv sync && pytest` installs h5py per pyproject.toml.
h5py = pytest.importorskip("h5py")

from dna_decode.models.cache import (  # noqa: E402
    CacheMetadata,
    EmbeddingCache,
    EmbeddingCacheError,
    EmbeddingCacheVersionMismatch,
)
from dna_decode.models.foundation import MockFoundationModel, ModelMetadata  # noqa: E402


# ---- construction + metadata ----


def test_cache_creates_new_file_with_metadata(tmp_path: Path):
    cache = EmbeddingCache(
        tmp_path / "cache.h5",
        model_name="mock",
        model_version="v1",
        embedding_dim=8,
    )
    md = cache.metadata()
    assert md.model_name == "mock"
    assert md.model_version == "v1"
    assert md.embedding_dim == 8
    assert md.created_at  # non-empty ISO timestamp


def test_cache_reopens_matching_existing_file(tmp_path: Path):
    cache_path = tmp_path / "cache.h5"
    EmbeddingCache(cache_path, "mock", "v1", 8)
    # Re-open with matching params → no error
    reopened = EmbeddingCache(cache_path, "mock", "v1", 8)
    assert reopened.metadata().model_version == "v1"


def test_cache_version_mismatch_raises(tmp_path: Path):
    cache_path = tmp_path / "cache.h5"
    EmbeddingCache(cache_path, "mock", "v1", 8)
    with pytest.raises(EmbeddingCacheVersionMismatch):
        EmbeddingCache(cache_path, "mock", "v2", 8)


def test_cache_dim_mismatch_raises(tmp_path: Path):
    cache_path = tmp_path / "cache.h5"
    EmbeddingCache(cache_path, "mock", "v1", 8)
    with pytest.raises(EmbeddingCacheVersionMismatch):
        EmbeddingCache(cache_path, "mock", "v1", 16)


# ---- put / get round-trip ----


def test_put_get_round_trip(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    emb = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
    cache.put("s001", "gene_X", emb)
    out = cache.get("s001", "gene_X")
    np.testing.assert_array_almost_equal(out, emb)


def test_put_wrong_shape_raises(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 8)
    with pytest.raises(EmbeddingCacheError, match="shape"):
        cache.put("s001", "gene_X", np.zeros(4))


def test_put_overwrites_existing(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("s001", "g", np.array([1, 2, 3, 4], dtype=np.float32))
    cache.put("s001", "g", np.array([9, 9, 9, 9], dtype=np.float32))
    out = cache.get("s001", "g")
    np.testing.assert_array_almost_equal(out, [9, 9, 9, 9])


def test_get_missing_returns_none(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    assert cache.get("nonexistent", "missing") is None


def test_has_returns_bool(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("s001", "g", np.zeros(4, dtype=np.float32))
    assert cache.has("s001", "g")
    assert not cache.has("s001", "absent")
    assert not cache.has("absent_strain", "g")


def test_gene_id_with_slash_is_sanitized(tmp_path: Path):
    """HDF5 forbids '/' in group names; cache should sanitize to '__'."""
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("s001", "gyrA/region", np.zeros(4, dtype=np.float32))
    assert cache.has("s001", "gyrA/region")


# ---- list operations ----


def test_list_strains_empty(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    assert cache.list_strains() == []


def test_list_strains_returns_sorted(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("zeta", "g", np.zeros(4, dtype=np.float32))
    cache.put("alpha", "g", np.zeros(4, dtype=np.float32))
    cache.put("beta", "g", np.zeros(4, dtype=np.float32))
    assert cache.list_strains() == ["alpha", "beta", "zeta"]


def test_list_genes_per_strain(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("s001", "gyrA", np.zeros(4, dtype=np.float32))
    cache.put("s001", "parC", np.zeros(4, dtype=np.float32))
    cache.put("s002", "blaCTX", np.zeros(4, dtype=np.float32))
    assert cache.list_genes("s001") == ["gyrA", "parC"]
    assert cache.list_genes("s002") == ["blaCTX"]
    assert cache.list_genes("missing") == []


# ---- bulk operations ----


def test_bulk_get_stacks_in_order(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("s001", "gyrA", np.array([1, 1, 1, 1], dtype=np.float32))
    cache.put("s002", "gyrA", np.array([2, 2, 2, 2], dtype=np.float32))
    cache.put("s003", "gyrA", np.array([3, 3, 3, 3], dtype=np.float32))

    out = cache.bulk_get([("s001", "gyrA"), ("s002", "gyrA"), ("s003", "gyrA")])
    assert out.shape == (3, 4)
    np.testing.assert_array_almost_equal(out[0], [1, 1, 1, 1])
    np.testing.assert_array_almost_equal(out[2], [3, 3, 3, 3])


def test_bulk_get_nan_fills_missing(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    cache.put("s001", "gyrA", np.array([1, 1, 1, 1], dtype=np.float32))
    out = cache.bulk_get([("s001", "gyrA"), ("s_missing", "gyrA")], fill_missing=True)
    np.testing.assert_array_almost_equal(out[0], [1, 1, 1, 1])
    assert np.all(np.isnan(out[1]))


def test_bulk_get_raises_when_missing_and_no_fill(tmp_path: Path):
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 4)
    with pytest.raises(EmbeddingCacheError, match="Missing"):
        cache.bulk_get([("absent", "g")], fill_missing=False)


# ---- populate end-to-end ----


def test_populate_writes_per_gene_per_strain(tmp_path: Path):
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)

    strain_sequences = {
        "s001": {"gyrA": "ATGCATGC", "parC": "GGGGAAAA"},
        "s002": {"gyrA": "TTTTCCCC"},
    }
    written = cache.populate(model, strain_sequences=strain_sequences)
    assert written == {"s001": 2, "s002": 1}
    assert cache.has("s001", "gyrA")
    assert cache.has("s001", "parC")
    assert cache.has("s002", "gyrA")


def test_populate_idempotent_skip_existing(tmp_path: Path):
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    seqs = {"s001": {"gyrA": "ATGC"}}

    first_pass = cache.populate(model, strain_sequences=seqs)
    second_pass = cache.populate(model, strain_sequences=seqs)
    assert first_pass["s001"] == 1
    assert second_pass["s001"] == 0  # skipped existing


def test_populate_progress_callback_fires_per_strain(tmp_path: Path):
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    seqs = {"s001": {"g": "ATGC"}, "s002": {"g": "GGGG"}}

    calls: list[tuple[str, int, int]] = []
    cache.populate(
        model,
        strain_sequences=seqs,
        progress_callback=lambda sid, done, total: calls.append((sid, done, total)),
    )
    assert len(calls) == 2
    assert ("s001", 1, 1) in calls
    assert ("s002", 1, 1) in calls


# ---- Wave 2.5 hardening C6: plan-contract signature + open-once HDF5 ----


def test_populate_plan_signature_extracts_cds_via_annotations(tmp_path: Path):
    """populate(strain_genomes, annotations) calls extract_cds_sequences internally."""
    import pandas as pd

    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)

    # Build a tiny FASTA + matching annotation table
    fasta = tmp_path / "s001.fna"
    fasta.write_text(">chr1\n" + "ATGC" * 25 + "\n", encoding="utf-8")

    annotations = pd.DataFrame(
        {
            "seqid": ["chr1"],
            "source": ["test"],
            "type": ["CDS"],
            "start": [1],
            "end": [12],
            "strand": ["+"],
            "gene_id": ["gene_A"],
            "locus_tag": [""],
            "product": [""],
        }
    )

    written = cache.populate(
        model,
        strain_genomes={"s001": fasta},
        annotations={"s001": annotations},
    )
    assert written["s001"] == 1
    assert cache.has("s001", "gene_A")


def test_populate_rejects_both_paths_provided(tmp_path: Path):
    """Cannot pass both strain_genomes AND strain_sequences."""
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    with pytest.raises(EmbeddingCacheError, match="not both"):
        cache.populate(
            model,
            strain_genomes={"s001": tmp_path / "x.fna"},
            strain_sequences={"s001": {"g": "ATGC"}},
        )


def test_populate_rejects_neither_path_provided(tmp_path: Path):
    """Must pass strain_genomes+annotations OR strain_sequences."""
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    with pytest.raises(EmbeddingCacheError, match="must provide"):
        cache.populate(model)


def test_populate_genomes_without_annotations_raises(tmp_path: Path):
    """strain_genomes requires matching annotations dict."""
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    with pytest.raises(EmbeddingCacheError, match="annotations"):
        cache.populate(model, strain_genomes={"s001": tmp_path / "x.fna"})


def test_populate_genomes_missing_annotation_for_strain_raises(tmp_path: Path):
    """annotations dict must include every strain in strain_genomes."""
    import pandas as pd

    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    fasta = tmp_path / "s001.fna"
    fasta.write_text(">chr1\nATGC\n", encoding="utf-8")
    with pytest.raises(EmbeddingCacheError, match="no annotation"):
        cache.populate(
            model,
            strain_genomes={"s001": fasta},
            annotations={"s002": pd.DataFrame()},  # wrong strain
        )
