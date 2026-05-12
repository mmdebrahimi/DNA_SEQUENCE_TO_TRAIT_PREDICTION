"""Step 2 — NCBI RefSeq downloader with on-disk caching.

Fetches FASTA + GFF3 + GenBank for an accession and caches under
`<cache_dir>/<accession>/`. Idempotent: short-circuits on existing cache
unless `force=True`. Atomic via `.complete` sentinel.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Iterable

import requests

DEFAULT_TIMEOUT_SECONDS = 30
DEFAULT_MAX_RETRIES = 3
DEFAULT_BACKOFF_SECONDS = 2
ECOLI_K12_MG1655_ACCESSION = "GCF_000005845.2"

# NCBI Datasets API for genome packages (FASTA + GFF + GenBank)
# Real endpoint format:
#   https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<accession>/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=GENOME_GBFF
# Returns a ZIP. For Phase 1 the wrapper exposes per-file accessors.

NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"


class RefSeqDownloadError(Exception):
    """Network / 5xx / timeout failure after retries."""


class RefSeqAccessionNotFound(Exception):
    """NCBI returned 404 for the accession."""


def default_ecoli_k12_accession() -> str:
    """Canonical Phase 1 reference: E. coli K-12 substr. MG1655."""
    return ECOLI_K12_MG1655_ACCESSION


def cache_dir_for(accession: str, cache_root: Path) -> Path:
    """Return the per-accession cache directory."""
    return Path(cache_root) / accession


def is_cache_complete(cache_path: Path) -> bool:
    """Check if the per-accession cache has the .complete sentinel."""
    return (cache_path / ".complete").exists()


def list_cached(cache_root: Path | str) -> list[str]:
    """Enumerate accessions with a complete cache under cache_root."""
    root = Path(cache_root)
    if not root.exists():
        return []
    return sorted(
        d.name
        for d in root.iterdir()
        if d.is_dir() and is_cache_complete(d)
    )


def _http_get_with_retry(
    url: str,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    max_retries: int = DEFAULT_MAX_RETRIES,
) -> requests.Response:
    """GET with exponential backoff retry on 5xx / network errors.

    404 raises immediately (not retried).
    """
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=timeout_seconds, stream=True)
        except (requests.Timeout, requests.ConnectionError) as e:
            last_err = e
            continue

        if resp.status_code == 404:
            raise RefSeqAccessionNotFound(f"NCBI 404 for {url}")
        if resp.status_code >= 500:
            last_err = RefSeqDownloadError(f"NCBI 5xx ({resp.status_code}) for {url}")
            continue
        if not resp.ok:
            raise RefSeqDownloadError(
                f"NCBI HTTP {resp.status_code} for {url}: {resp.text[:200]}"
            )
        return resp

    raise RefSeqDownloadError(
        f"NCBI download failed after {max_retries} attempts: {last_err}"
    )


def _build_dataset_url(accession: str) -> str:
    """Compose the NCBI Datasets v2 genome-package URL for an accession."""
    base = f"{NCBI_DATASETS_BASE}/genome/accession/{accession}/download"
    params = (
        "include_annotation_type=GENOME_FASTA"
        "&include_annotation_type=GENOME_GFF"
        "&include_annotation_type=GENOME_GBFF"
    )
    return f"{base}?{params}"


def _write_atomic(target: Path, content: bytes | str) -> None:
    """Write file content via tmp → rename for atomicity."""
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".tmp")
    mode = "wb" if isinstance(content, bytes) else "w"
    encoding = None if isinstance(content, bytes) else "utf-8"
    with open(tmp, mode, encoding=encoding) as f:
        f.write(content)
    tmp.replace(target)


def download_genome(
    accession: str,
    cache_root: Path | str,
    force: bool = False,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    max_retries: int = DEFAULT_MAX_RETRIES,
) -> Path:
    """Download FASTA + GFF3 + GenBank for an accession; cache idempotently.

    Returns the per-accession cache directory. Subsequent calls short-circuit
    unless `force=True` or the cache lacks a `.complete` sentinel (interpreted
    as partial cache → re-download).
    """
    cache_path = cache_dir_for(accession, Path(cache_root))

    if not force and is_cache_complete(cache_path):
        return cache_path

    # Partial cache (no sentinel) → wipe + redownload
    if cache_path.exists() and not is_cache_complete(cache_path):
        shutil.rmtree(cache_path)
    cache_path.mkdir(parents=True, exist_ok=True)

    url = _build_dataset_url(accession)
    resp = _http_get_with_retry(url, timeout_seconds=timeout_seconds, max_retries=max_retries)

    # NCBI returns a ZIP archive containing fna/gff/gbff. For the v1 scaffold we
    # write the raw response body and defer the unzip to the consumer
    # (annotations.py in Step 3 will handle GFF parsing; refseq stores raw files).
    # Real implementation hook: real ingestion writes 3 separate files from the ZIP.
    archive_path = cache_path / "package.zip"
    _write_atomic(archive_path, resp.content)

    # Sentinel last → marks atomicity-complete
    _write_atomic(cache_path / ".complete", "")

    return cache_path


def fasta_path(accession: str, cache_root: Path | str) -> Path:
    """Path to the FASTA file inside the per-accession cache."""
    return cache_dir_for(accession, Path(cache_root)) / "genome.fna"


def gff_path(accession: str, cache_root: Path | str) -> Path:
    """Path to the GFF3 file inside the per-accession cache."""
    return cache_dir_for(accession, Path(cache_root)) / "annotations.gff3"


def genbank_path(accession: str, cache_root: Path | str) -> Path:
    """Path to the GenBank file inside the per-accession cache."""
    return cache_dir_for(accession, Path(cache_root)) / "annotations.gbk"
