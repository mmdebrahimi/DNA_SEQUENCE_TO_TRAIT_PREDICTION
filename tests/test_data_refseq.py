"""Tests for Step 2 — NCBI RefSeq downloader."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from dna_decode.data import refseq


# ---- helpers ----


def _mock_ok_response(content: bytes = b"FAKE-GENOME-PACKAGE-ZIP-BYTES") -> MagicMock:
    resp = MagicMock(spec=requests.Response)
    resp.status_code = 200
    resp.ok = True
    resp.content = content
    return resp


def _mock_404_response() -> MagicMock:
    resp = MagicMock(spec=requests.Response)
    resp.status_code = 404
    resp.ok = False
    resp.content = b""
    resp.text = "Not Found"
    return resp


def _mock_500_response() -> MagicMock:
    resp = MagicMock(spec=requests.Response)
    resp.status_code = 500
    resp.ok = False
    resp.text = "Internal Server Error"
    return resp


# ---- canonical-accession + path helpers ----


def test_default_ecoli_k12_accession():
    assert refseq.default_ecoli_k12_accession() == "GCF_000005845.2"


def test_cache_dir_for_compose(tmp_path: Path):
    p = refseq.cache_dir_for("GCF_000005845.2", tmp_path)
    assert p == tmp_path / "GCF_000005845.2"


def test_path_helpers_compose(tmp_path: Path):
    assert refseq.fasta_path("GCF_X", tmp_path).name == "genome.fna"
    assert refseq.gff_path("GCF_X", tmp_path).name == "annotations.gff3"
    assert refseq.genbank_path("GCF_X", tmp_path).name == "annotations.gbk"


# ---- download_genome happy path ----


def test_download_genome_writes_package_and_sentinel(tmp_path: Path):
    with patch.object(refseq, "_http_get_with_retry", return_value=_mock_ok_response()):
        out = refseq.download_genome("GCF_000005845.2", tmp_path)

    assert out == tmp_path / "GCF_000005845.2"
    assert (out / "package.zip").exists()
    assert (out / ".complete").exists()
    assert refseq.is_cache_complete(out)


def test_download_genome_idempotent_short_circuits(tmp_path: Path):
    """Second call with complete cache short-circuits without hitting HTTP."""
    fetcher = MagicMock(return_value=_mock_ok_response())
    with patch.object(refseq, "_http_get_with_retry", fetcher):
        refseq.download_genome("GCF_X", tmp_path)
        assert fetcher.call_count == 1
        refseq.download_genome("GCF_X", tmp_path)  # should short-circuit
        assert fetcher.call_count == 1


def test_download_genome_force_redownloads(tmp_path: Path):
    """force=True bypasses short-circuit."""
    fetcher = MagicMock(return_value=_mock_ok_response())
    with patch.object(refseq, "_http_get_with_retry", fetcher):
        refseq.download_genome("GCF_X", tmp_path)
        refseq.download_genome("GCF_X", tmp_path, force=True)
        assert fetcher.call_count == 2


def test_download_genome_partial_cache_redownloads(tmp_path: Path):
    """Cache with no .complete sentinel triggers re-download (treated as partial)."""
    # Pre-create partial cache
    cache = tmp_path / "GCF_X"
    cache.mkdir()
    (cache / "package.zip").write_bytes(b"PARTIAL")
    assert not refseq.is_cache_complete(cache)

    fetcher = MagicMock(return_value=_mock_ok_response(b"FRESH"))
    with patch.object(refseq, "_http_get_with_retry", fetcher):
        refseq.download_genome("GCF_X", tmp_path)

    assert fetcher.call_count == 1
    assert (cache / "package.zip").read_bytes() == b"FRESH"
    assert refseq.is_cache_complete(cache)


# ---- error paths ----


def test_download_genome_404_raises(tmp_path: Path):
    with patch("requests.get", return_value=_mock_404_response()):
        with pytest.raises(refseq.RefSeqAccessionNotFound):
            refseq.download_genome("GCF_BOGUS", tmp_path)


def test_download_genome_5xx_after_retries_raises(tmp_path: Path):
    with patch("requests.get", return_value=_mock_500_response()):
        with pytest.raises(refseq.RefSeqDownloadError, match="failed after"):
            refseq.download_genome("GCF_X", tmp_path, max_retries=2)


def test_download_genome_timeout_after_retries_raises(tmp_path: Path):
    with patch("requests.get", side_effect=requests.Timeout("read timeout")):
        with pytest.raises(refseq.RefSeqDownloadError, match="failed after"):
            refseq.download_genome("GCF_X", tmp_path, max_retries=2)


# ---- list_cached ----


def test_list_cached_empty_when_no_cache(tmp_path: Path):
    assert refseq.list_cached(tmp_path / "does-not-exist") == []


def test_list_cached_returns_only_complete(tmp_path: Path):
    # Create two cache dirs: one complete, one partial
    (tmp_path / "GCF_A").mkdir()
    (tmp_path / "GCF_A" / ".complete").touch()
    (tmp_path / "GCF_B").mkdir()
    (tmp_path / "GCF_B" / "package.zip").touch()  # partial — no sentinel

    cached = refseq.list_cached(tmp_path)
    assert cached == ["GCF_A"]


def test_url_builder_includes_required_annotations():
    url = refseq._build_dataset_url("GCF_000005845.2")
    assert "GCF_000005845.2" in url
    assert "GENOME_FASTA" in url
    assert "GENOME_GFF" in url
    assert "GENOME_GBFF" in url
