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
    CompletenessReport,
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


# ---- pooling_strategy (Phase 2.5 hardening) ----


def test_cache_records_default_pooling_strategy(tmp_path: Path):
    """New caches default to 'single_seq_mean' for backward compat."""
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", 8)
    assert cache.metadata().pooling_strategy == "single_seq_mean"


def test_cache_records_custom_pooling_strategy(tmp_path: Path):
    cache = EmbeddingCache(
        tmp_path / "cache.h5",
        "mock",
        "v1",
        8,
        pooling_strategy="mask_aware_mean",
    )
    assert cache.metadata().pooling_strategy == "mask_aware_mean"


def test_cache_pooling_mismatch_raises(tmp_path: Path):
    """Reopening with different pooling_strategy → mismatch (prevents hybrid cache)."""
    cache_path = tmp_path / "cache.h5"
    EmbeddingCache(cache_path, "mock", "v1", 8, pooling_strategy="single_seq_mean")
    with pytest.raises(EmbeddingCacheVersionMismatch, match="pooling"):
        EmbeddingCache(cache_path, "mock", "v1", 8, pooling_strategy="mask_aware_mean")


def test_cache_legacy_file_without_pooling_attr_accepts_default(tmp_path: Path):
    """Caches written pre-2026-05-13 lack the pooling attr; default-match should pass."""
    import h5py
    cache_path = tmp_path / "legacy.h5"
    # Hand-craft an old-shape cache: attrs match defaults except no pooling_strategy
    with h5py.File(cache_path, "w") as f:
        f.attrs["model_name"] = "mock"
        f.attrs["model_version"] = "v1"
        f.attrs["embedding_dim"] = 8
        f.attrs["created_at"] = "2026-01-01T00:00:00"
        f.create_group("strains")
        # NOTE: no pooling_strategy attribute

    # Default constructor (pooling_strategy="single_seq_mean") should match a legacy
    # file without the attr — the reader defaults to "single_seq_mean".
    reopened = EmbeddingCache(cache_path, "mock", "v1", 8)
    assert reopened.metadata().pooling_strategy == "single_seq_mean"


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


def test_populate_flushes_hdf5_per_strain_for_usb_hiccup_resilience(tmp_path: Path, monkeypatch):
    """Per-strain f.flush() pins the HDF5 superblock + EOA forward to disk so a
    USB hiccup (Seagate Portable D: drive observed 2026-05-14) costs at most
    one strain's recompute, not the whole populate. Regression guard against
    silent EOA-truncation when __exit__ flush never runs.

    Strategy: monkeypatch h5py.File so we can record `flush()` invocations.
    Expect one flush per strain, in strain order.
    """
    import h5py
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=4, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=4)

    # Wrap h5py.File to record per-flush calls
    flush_count = {"n": 0}
    original_file_class = h5py.File

    class _FlushSpyingFile(original_file_class):
        def flush(self):
            flush_count["n"] += 1
            return super().flush()

    monkeypatch.setattr(cache, "_h5py", type("M", (), {"File": _FlushSpyingFile}))

    seqs = {"s001": {"g": "ATGC"}, "s002": {"g": "GGGG"}, "s003": {"g": "AAAA"}}
    cache.populate(model, strain_sequences=seqs)
    # 3 strains -> at least 3 flushes (one per strain). >=3 because the with-block
    # exit on the h5py File context manager may also implicitly flush.
    assert flush_count["n"] >= 3, f"Expected >=3 per-strain flushes, got {flush_count['n']}"


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


# ---- Phase 2.5 batching refactor (EMBED_BATCH_SIZE=4 chunking) ----


def test_populate_chunks_in_groups_of_four(tmp_path: Path):
    """populate() must chunk pending sequences in groups of EMBED_BATCH_SIZE=4.

    Guards the `for i in range(0, len(pending), EMBED_BATCH_SIZE)` loop:
    10 pending genes → embed_batch called 3 times with sizes [4, 4, 2].
    """
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)

    # Spy on embed_batch
    real_embed_batch = model.embed_batch
    observed_batch_sizes: list[int] = []

    def spy(seqs):
        observed_batch_sizes.append(len(seqs))
        return real_embed_batch(seqs)

    model.embed_batch = spy

    strain_sequences = {
        "s001": {f"g{i:02d}": f"ATGC{i:04d}" for i in range(10)},
    }
    written = cache.populate(model, strain_sequences=strain_sequences)
    assert written == {"s001": 10}
    assert observed_batch_sizes == [4, 4, 2]
    # All 10 datasets present
    assert len(cache.list_genes("s001")) == 10


def test_populate_skip_existing_only_embeds_pending(tmp_path: Path):
    """skip_existing=True must filter the pending list BEFORE batching, so
    already-cached genes are not re-embedded.
    """
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    # Pre-populate 2 genes
    cache.put("s001", "g00", np.zeros(8, dtype=np.float32))
    cache.put("s001", "g01", np.zeros(8, dtype=np.float32))

    observed_sequences: list[list[str]] = []
    real_embed_batch = model.embed_batch

    def spy(seqs):
        observed_sequences.append(list(seqs))
        return real_embed_batch(seqs)

    model.embed_batch = spy

    strain_sequences = {
        "s001": {f"g{i:02d}": f"ATGC{i:04d}" for i in range(5)},  # 5 genes; first 2 cached
    }
    written = cache.populate(model, strain_sequences=strain_sequences, skip_existing=True)
    assert written == {"s001": 3}  # only g02, g03, g04
    # One chunk of 3 (all under EMBED_BATCH_SIZE=4)
    assert len(observed_sequences) == 1
    assert len(observed_sequences[0]) == 3


def test_populate_rejects_bad_batch_shape(tmp_path: Path):
    """populate() must validate the model's embed_batch return shape is (B, D);
    a stale 1-D return (regression to pre-batching contract) must raise.
    """

    class BadBatchModel(MockFoundationModel):
        def embed_batch(self, sequences):
            # Wrong shape: returns 1-D instead of (B, D)
            return np.zeros(self.embedding_dim, dtype=np.float32)

    model = BadBatchModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    with pytest.raises(EmbeddingCacheError, match="embedding batch shape"):
        cache.populate(model, strain_sequences={"s001": {"g0": "ATGC"}})


def test_populate_rejects_wrong_embedding_dim_in_batch(tmp_path: Path):
    """populate() must reject batch outputs whose dim axis disagrees with the cache."""

    class WrongDimModel(MockFoundationModel):
        def embed_batch(self, sequences):
            # Right rank (2-D), wrong inner dim
            return np.zeros((len(sequences), 16), dtype=np.float32)

    model = WrongDimModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)
    with pytest.raises(EmbeddingCacheError, match="embedding batch shape"):
        cache.populate(model, strain_sequences={"s001": {"g0": "ATGC"}})


def test_populate_preserves_gene_to_embedding_correspondence(tmp_path: Path):
    """The (gene_id, embedding) zip after batching must preserve order — each
    gene must be stored against the embedding produced from its OWN sequence,
    not a neighbor's. Regression guard: if the zip ever drifted (e.g., from a
    shuffled batch), gene_X would silently get gene_Y's embedding.
    """
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=8)

    # 6 genes spans 2 chunks (4 + 2) — exercises both the full-chunk and
    # short-final-chunk paths.
    gene_to_seq = {f"g{i:02d}": f"GENE{i:04d}AAAA" for i in range(6)}
    cache.populate(model, strain_sequences={"s001": gene_to_seq})

    # For the mock, the per-sequence embedding is deterministic from the
    # sequence bytes alone — assert each cached embedding matches what the
    # mock produces for that gene's sequence.
    for gene_id, seq in gene_to_seq.items():
        cached = cache.get("s001", gene_id)
        expected = model._embed_window(seq)
        np.testing.assert_array_equal(cached, expected)


# ---- verify_complete (crash-truncated-strain integrity gate) ----


def _seed_cache(tmp_path: Path, payload: dict[str, dict[str, np.ndarray]]) -> EmbeddingCache:
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=4)
    for sid, genes in payload.items():
        for gid, arr in genes.items():
            cache.put(sid, gid, arr.astype(np.float32))
    return cache


def test_verify_complete_all_present(tmp_path: Path):
    cache = _seed_cache(
        tmp_path,
        {"s1": {"g1": np.ones(4), "g2": np.full(4, 2.0)}},
    )
    report = cache.verify_complete({"s1": {"g1", "g2"}})
    assert report.status == {"s1": "complete"}
    assert report.all_complete is True
    assert report.counts == {"complete": 1}
    assert report.missing_genes["s1"] == []
    assert report.corrupt_details["s1"] == []


def test_verify_complete_partial_strain_flagged(tmp_path: Path):
    # Half-flushed strain: 1 of 3 expected genes cached.
    cache = _seed_cache(tmp_path, {"s1": {"g1": np.ones(4)}})
    report = cache.verify_complete({"s1": {"g1", "g2", "g3"}})
    assert report.status == {"s1": "partial"}
    assert report.all_complete is False
    assert sorted(report.missing_genes["s1"]) == ["g2", "g3"]
    assert report.cached_n["s1"] == 1
    assert report.expected_n["s1"] == 3


def test_verify_complete_absent_strain(tmp_path: Path):
    # Strain expected by cohort but never written to cache.
    cache = _seed_cache(tmp_path, {"s1": {"g1": np.ones(4)}})
    report = cache.verify_complete({"s1": {"g1"}, "s2": {"gA", "gB"}})
    assert report.status["s1"] == "complete"
    assert report.status["s2"] == "absent"
    assert report.all_complete is False
    assert sorted(report.missing_genes["s2"]) == ["gA", "gB"]
    assert report.cached_n["s2"] == 0


def test_verify_complete_corrupt_non_finite_values(tmp_path: Path):
    nan_emb = np.array([1.0, np.nan, 3.0, 4.0])
    cache = _seed_cache(tmp_path, {"s1": {"g1": nan_emb, "g2": np.ones(4)}})
    report = cache.verify_complete({"s1": {"g1", "g2"}})
    assert report.status == {"s1": "corrupt"}
    assert report.all_complete is False
    assert any("non-finite" in d for d in report.corrupt_details["s1"])


def test_verify_complete_corrupt_wrong_shape(tmp_path: Path):
    # Bypass cache.put (which validates shape) to seed a corrupt dataset
    # directly via h5py — simulates a half-written embedding row.
    cache = EmbeddingCache(tmp_path / "cache.h5", "mock", "v1", embedding_dim=4)
    cache.put("s1", "g1", np.ones(4, dtype=np.float32))
    with h5py.File(cache.path, "a") as f:
        f.create_dataset("strains/s1/g2", data=np.ones(2, dtype=np.float32))
    report = cache.verify_complete({"s1": {"g1", "g2"}})
    assert report.status == {"s1": "corrupt"}
    assert any("shape=" in d for d in report.corrupt_details["s1"])


def test_verify_complete_stage1_consumer_must_refuse_on_partial(tmp_path: Path):
    # Regression pin for the half-flushed-strain landmine: Stage 1's loader
    # admits a strain on >=1 cached gene + mean-pools whatever's there. The
    # `all_complete` property is the consumer-side refuse signal.
    cache = _seed_cache(tmp_path, {"s1": {"g1": np.ones(4)}})
    report = cache.verify_complete({"s1": {"g1", "g2"}})
    assert report.all_complete is False  # consumer MUST bail


def test_verify_complete_empty_input_is_not_complete(tmp_path: Path):
    # Defensive: empty expected-genes dict shouldn't trick all_complete into True.
    cache = _seed_cache(tmp_path, {"s1": {"g1": np.ones(4)}})
    report = cache.verify_complete({})
    assert report.all_complete is False
    assert report.counts == {}


def test_verify_complete_returns_completeness_report_dataclass(tmp_path: Path):
    cache = _seed_cache(tmp_path, {"s1": {"g1": np.ones(4)}})
    report = cache.verify_complete({"s1": {"g1"}})
    assert isinstance(report, CompletenessReport)
