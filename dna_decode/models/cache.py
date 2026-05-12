"""Step 8 — HDF5-backed embedding cache.

Pre-compute foundation-model embeddings once per (strain × gene); cache to
disk; reuse across all CV folds + leaderboard runs. Without this, every
Wave 3 / Wave 6 run would re-invoke the foundation model on millions of
sequences.

HDF5 layout:
    /                              file-level attrs: model_name, model_version,
                                                     embedding_dim, created_at
    /strains/<strain_id>/<gene_id> dataset: 1-D float32 array (embedding_dim,)

Version-mismatch refusal: opening a cache with a model_version different from
the wrapping handle's expected version raises EmbeddingCacheVersionMismatch
rather than silently overwriting. v0.2 may add a migration path.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


class EmbeddingCacheError(Exception):
    """Generic cache-handling failure."""


class EmbeddingCacheVersionMismatch(EmbeddingCacheError):
    """Cache exists with a different model_version than expected."""


@dataclass(frozen=True)
class CacheMetadata:
    """File-level metadata stored as HDF5 root attributes."""

    model_name: str
    model_version: str
    embedding_dim: int
    created_at: str  # ISO-8601 string


class EmbeddingCache:
    """HDF5 wrapper for (strain, gene) -> embedding lookups.

    Construction is non-destructive: if `path` exists, opens for append; else
    creates a new HDF5 file with the given metadata. Version mismatch raises
    instead of overwriting.
    """

    def __init__(
        self,
        path: Path | str,
        model_name: str,
        model_version: str,
        embedding_dim: int,
    ):
        # Import h5py lazily so this module imports cleanly without h5py installed
        # (full pytest run + actual cache use requires `uv sync`).
        try:
            import h5py
        except ImportError as e:
            raise EmbeddingCacheError(
                "h5py not installed; run `uv sync` to install Phase 1 deps"
            ) from e

        self._h5py = h5py
        self.path = Path(path)
        self.model_name = model_name
        self.model_version = model_version
        self.embedding_dim = embedding_dim

        if self.path.exists():
            self._validate_existing_metadata()
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._initialize_new_file()

    # ---- file lifecycle ----

    def _initialize_new_file(self) -> None:
        with self._h5py.File(self.path, "w") as f:
            f.attrs["model_name"] = self.model_name
            f.attrs["model_version"] = self.model_version
            f.attrs["embedding_dim"] = self.embedding_dim
            f.attrs["created_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            f.create_group("strains")

    def _validate_existing_metadata(self) -> None:
        with self._h5py.File(self.path, "r") as f:
            existing_name = f.attrs.get("model_name", "")
            existing_version = f.attrs.get("model_version", "")
            existing_dim = int(f.attrs.get("embedding_dim", -1))
            if (
                existing_name != self.model_name
                or existing_version != self.model_version
                or existing_dim != self.embedding_dim
            ):
                raise EmbeddingCacheVersionMismatch(
                    f"Cache at {self.path} has model_name={existing_name!r}, "
                    f"model_version={existing_version!r}, embedding_dim={existing_dim}; "
                    f"expected name={self.model_name!r}, "
                    f"version={self.model_version!r}, dim={self.embedding_dim}. "
                    f"Delete the cache + re-populate, OR use a migration tool."
                )

    def metadata(self) -> CacheMetadata:
        with self._h5py.File(self.path, "r") as f:
            return CacheMetadata(
                model_name=str(f.attrs.get("model_name", "")),
                model_version=str(f.attrs.get("model_version", "")),
                embedding_dim=int(f.attrs.get("embedding_dim", -1)),
                created_at=str(f.attrs.get("created_at", "")),
            )

    # ---- per-entry operations ----

    @staticmethod
    def _dataset_path(strain_id: str, gene_id: str) -> str:
        # HDF5 group names must not contain '/'. Sanitize gene_id by replacing
        # any '/' with '__'. strain_id sanitization left to caller (strain IDs
        # in NCBI accession form don't have slashes).
        sanitized_gene = gene_id.replace("/", "__")
        return f"strains/{strain_id}/{sanitized_gene}"

    def put(self, strain_id: str, gene_id: str, embedding: np.ndarray) -> None:
        if embedding.ndim != 1 or embedding.shape[0] != self.embedding_dim:
            raise EmbeddingCacheError(
                f"embedding shape {embedding.shape} mismatch with cache "
                f"embedding_dim={self.embedding_dim}"
            )
        ds_path = self._dataset_path(strain_id, gene_id)
        with self._h5py.File(self.path, "a") as f:
            if ds_path in f:
                del f[ds_path]
            f.create_dataset(ds_path, data=embedding.astype(np.float32))

    def get(self, strain_id: str, gene_id: str) -> np.ndarray | None:
        ds_path = self._dataset_path(strain_id, gene_id)
        with self._h5py.File(self.path, "r") as f:
            if ds_path not in f:
                return None
            return f[ds_path][()].astype(np.float32)

    def has(self, strain_id: str, gene_id: str) -> bool:
        ds_path = self._dataset_path(strain_id, gene_id)
        with self._h5py.File(self.path, "r") as f:
            return ds_path in f

    def list_strains(self) -> list[str]:
        with self._h5py.File(self.path, "r") as f:
            if "strains" not in f:
                return []
            return sorted(f["strains"].keys())

    def list_genes(self, strain_id: str) -> list[str]:
        with self._h5py.File(self.path, "r") as f:
            grp_path = f"strains/{strain_id}"
            if grp_path not in f:
                return []
            return sorted(f[grp_path].keys())

    # ---- bulk operations ----

    def bulk_get(
        self,
        pairs: list[tuple[str, str]],
        fill_missing: bool = True,
    ) -> np.ndarray:
        """Stack embeddings for `pairs` into a (len(pairs), embedding_dim) array.

        Missing pairs are NaN-filled when `fill_missing=True`; otherwise raise.
        """
        n = len(pairs)
        out = np.full((n, self.embedding_dim), np.nan, dtype=np.float32)
        with self._h5py.File(self.path, "r") as f:
            for i, (strain_id, gene_id) in enumerate(pairs):
                ds_path = self._dataset_path(strain_id, gene_id)
                if ds_path in f:
                    out[i] = f[ds_path][()].astype(np.float32)
                elif not fill_missing:
                    raise EmbeddingCacheError(
                        f"Missing embedding for {strain_id}/{gene_id} "
                        f"(use fill_missing=True to NaN-fill)"
                    )
        return out

    # ---- end-to-end population ----

    def populate(
        self,
        model,  # FoundationModel (avoid hard import to keep dep light)
        strain_sequences: dict[str, dict[str, str]],
        skip_existing: bool = True,
        progress_callback=None,
    ) -> dict[str, int]:
        """Compute + cache embeddings for every (strain, gene) sequence pair.

        Args:
            model: FoundationModel instance from Step 7.
            strain_sequences: mapping strain_id -> { gene_id -> sequence }.
            skip_existing: if True, skip pairs already in the cache (idempotent).
            progress_callback: optional callable(strain_id, n_done, n_total).

        Returns:
            Per-strain dict of how many new entries were written.
        """
        written_per_strain: dict[str, int] = {}
        for strain_id, gene_map in strain_sequences.items():
            count = 0
            for gene_id, sequence in gene_map.items():
                if skip_existing and self.has(strain_id, gene_id):
                    continue
                # embed_batch returns (len(sequences), embedding_dim); we want
                # the single-row mean-pooled embedding
                emb = model.embed_batch([sequence])[0]
                self.put(strain_id, gene_id, emb)
                count += 1
            written_per_strain[strain_id] = count
            if progress_callback:
                progress_callback(strain_id, count, len(gene_map))
        return written_per_strain
