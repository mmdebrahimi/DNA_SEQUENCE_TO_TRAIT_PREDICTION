"""Step 2 — NCBI RefSeq downloader with on-disk caching.

Fetches FASTA + GFF3 + GenBank for an accession and caches under
`<cache_dir>/<accession>/`. Idempotent: short-circuits on existing cache
unless `force=True`. Atomic via `.complete` sentinel.

Wave 1.5 hardening: after downloading the NCBI Datasets ZIP, unpacks it
into `genome.fna`, `annotations.gff3`, `annotations.gbk` so the path
helpers point at real consumable files (not just the raw archive).
"""
from __future__ import annotations

import shutil
import zipfile
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

    Writes the raw NCBI Datasets ZIP to `package.zip`, then unpacks the
    expected genomic files into `genome.fna`, `annotations.gff3`, and
    `annotations.gbk`. Subsequent calls short-circuit unless `force=True` or
    the cache lacks a `.complete` sentinel (treated as partial → re-download).
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

    archive_path = cache_path / "package.zip"
    _write_atomic(archive_path, resp.content)

    # Wave 1.5: unpack the ZIP into the 3 canonical files
    _unpack_ncbi_datasets_zip(archive_path, cache_path)

    # Sentinel last → marks atomicity-complete
    _write_atomic(cache_path / ".complete", "")

    return cache_path


def _unpack_ncbi_datasets_zip(archive_path: Path, cache_path: Path) -> None:
    """Extract genome FASTA + GFF + GBFF from an NCBI Datasets ZIP archive.

    NCBI Datasets ZIPs follow a known structure:
        ncbi_dataset/data/<accession>/<accession>_genomic.fna
        ncbi_dataset/data/<accession>/genomic.gff
        ncbi_dataset/data/<accession>/genomic.gbff
    (plus a `README.md`, `dataset_catalog.json`, and other metadata.)

    Writes to:
        cache_path/genome.fna
        cache_path/annotations.gff3
        cache_path/annotations.gbk

    Missing-file behavior: each target file is best-effort. The annotation
    parser (Step 3) will fail gracefully if a specific format is absent — for
    instance, some assemblies ship without GBFF.
    """
    if not zipfile.is_zipfile(archive_path):
        raise RefSeqDownloadError(
            f"Downloaded archive at {archive_path} is not a valid ZIP. "
            "NCBI may have returned an error page or partial file."
        )

    target_paths = {
        ".fna": cache_path / "genome.fna",
        ".gff": cache_path / "annotations.gff3",
        ".gbff": cache_path / "annotations.gbk",
    }
    extracted: dict[str, bool] = {".fna": False, ".gff": False, ".gbff": False}

    with zipfile.ZipFile(archive_path) as zf:
        for member in zf.namelist():
            lower = member.lower()
            for suffix, target in target_paths.items():
                # Match the first occurrence per suffix; NCBI sometimes ships
                # both protein.faa and genomic.fna — we want genomic only.
                if extracted[suffix]:
                    continue
                if not lower.endswith(suffix):
                    continue
                if suffix == ".fna" and "protein" in lower:
                    # Skip protein FASTA; we want the nucleotide genome
                    continue
                with zf.open(member) as src:
                    content = src.read()
                target.write_bytes(content)
                extracted[suffix] = True

    if not extracted[".fna"]:
        # No FASTA at all → the cache is unusable downstream
        raise RefSeqDownloadError(
            f"NCBI archive at {archive_path} contained no genomic .fna file"
        )


def fasta_path(accession: str, cache_root: Path | str) -> Path:
    """Path to the FASTA file inside the per-accession cache."""
    return cache_dir_for(accession, Path(cache_root)) / "genome.fna"


def gff_path(accession: str, cache_root: Path | str) -> Path:
    """Path to the GFF3 file inside the per-accession cache."""
    return cache_dir_for(accession, Path(cache_root)) / "annotations.gff3"


def genbank_path(accession: str, cache_root: Path | str) -> Path:
    """Path to the GenBank file inside the per-accession cache."""
    return cache_dir_for(accession, Path(cache_root)) / "annotations.gbk"
