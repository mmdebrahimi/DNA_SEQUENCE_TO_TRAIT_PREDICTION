"""Step 8 — HDF5-backed embedding cache.

Pre-compute foundation-model embeddings once per (strain × gene); cache to
disk; reuse across all CV folds + leaderboard runs. Without this, every
Wave 3 / Wave 6 run would re-invoke the foundation model on millions of
sequences.

HDF5 layout:
    /                              file-level attrs: model_name, model_version,
                                                     embedding_dim, created_at,
                                                     pooling_strategy
    /strains/<strain_id>/<gene_id> dataset: 1-D float32 array (embedding_dim,)

Version-mismatch refusal: opening a cache with a model_version different from
the wrapping handle's expected version raises EmbeddingCacheVersionMismatch
rather than silently overwriting. v0.2 may add a migration path.

`pooling_strategy` (added 2026-05-13 per /brainstorm D-2): captures HOW the
foundation model reduces token-level outputs to a per-window embedding. Today
the single-sequence inference path uses "single_seq_mean" (mean over every
token in the tokenizer output, including special tokens). A future batched
implementation must use mask-aware mean pooling and tag itself "mask_aware_mean"
to prevent hybrid caches (some genes embedded under strategy_A, others under
strategy_B, with no way to tell which). Strict-match on reopen.
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


DEFAULT_POOLING_STRATEGY = "single_seq_mean"
"""Tag for the current foundation-model embedding path: tokenize one sequence at
a time + mean over every output token position (incl. special tokens). All
Phase 1 + early Phase 2 callers inherit this default. Future batched / mask-
aware pooling must use a different tag (e.g., "mask_aware_mean") so the cache's
strict-match validation prevents hybrid embeddings on rebuilt caches."""


@dataclass(frozen=True)
class CacheMetadata:
    """File-level metadata stored as HDF5 root attributes."""

    model_name: str
    model_version: str
    embedding_dim: int
    created_at: str  # ISO-8601 string
    pooling_strategy: str = DEFAULT_POOLING_STRATEGY


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
        pooling_strategy: str = DEFAULT_POOLING_STRATEGY,
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
        self.pooling_strategy = pooling_strategy

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
            f.attrs["pooling_strategy"] = self.pooling_strategy
            f.create_group("strains")

    def _validate_existing_metadata(self) -> None:
        with self._h5py.File(self.path, "r") as f:
            existing_name = f.attrs.get("model_name", "")
            existing_version = f.attrs.get("model_version", "")
            existing_dim = int(f.attrs.get("embedding_dim", -1))
            # Legacy caches written before 2026-05-13 lack this attr; default
            # to the historical "single_seq_mean" value rather than raising.
            existing_pooling = str(
                f.attrs.get("pooling_strategy", DEFAULT_POOLING_STRATEGY)
            )
            if (
                existing_name != self.model_name
                or existing_version != self.model_version
                or existing_dim != self.embedding_dim
                or existing_pooling != self.pooling_strategy
            ):
                raise EmbeddingCacheVersionMismatch(
                    f"Cache at {self.path} has model_name={existing_name!r}, "
                    f"model_version={existing_version!r}, embedding_dim={existing_dim}, "
                    f"pooling_strategy={existing_pooling!r}; "
                    f"expected name={self.model_name!r}, "
                    f"version={self.model_version!r}, dim={self.embedding_dim}, "
                    f"pooling={self.pooling_strategy!r}. "
                    f"Delete the cache + re-populate, OR use a migration tool."
                )

    def metadata(self) -> CacheMetadata:
        with self._h5py.File(self.path, "r") as f:
            return CacheMetadata(
                model_name=str(f.attrs.get("model_name", "")),
                model_version=str(f.attrs.get("model_version", "")),
                embedding_dim=int(f.attrs.get("embedding_dim", -1)),
                created_at=str(f.attrs.get("created_at", "")),
                pooling_strategy=str(
                    f.attrs.get("pooling_strategy", DEFAULT_POOLING_STRATEGY)
                ),
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

    # ---- end-to-end population (plan contract signature) ----

    def populate(
        self,
        model,  # FoundationModel from Step 7 (avoid hard import)
        strain_genomes: dict[str, "Path"] | None = None,
        annotations: dict[str, "pd.DataFrame"] | None = None,
        strain_sequences: dict[str, dict[str, str]] | None = None,
        skip_existing: bool = True,
        progress_callback=None,
    ) -> dict[str, int]:
        """Compute + cache embeddings for every (strain, gene) sequence pair.

        Wave 2.5 hardening C6 — restored plan contract: takes either
        `strain_genomes` + `annotations` (preferred, runs Step 3 CDS extraction
        internally) OR `strain_sequences` (pre-extracted, back-compat). Opens
        the HDF5 file ONCE for the entire batch — was 4.5M opens per Phase 1
        run before, now 1 open per populate() call.

        Args:
            model: FoundationModel instance from Step 7.
            strain_genomes: mapping strain_id -> FASTA path. When provided,
                `annotations` must also be provided. Calls
                `dna_decode.data.annotations.extract_cds_sequences` per strain.
            annotations: mapping strain_id -> AnnotationTable (from Step 3
                `parse_gff3` / `parse_genbank`).
            strain_sequences: legacy entry — mapping strain_id -> { gene_id -> sequence }.
                Mutually exclusive with `strain_genomes`.
            skip_existing: if True, skip pairs already cached (idempotent).
            progress_callback: optional callable(strain_id, n_written, n_total).

        Returns:
            Per-strain dict of how many new entries were written.
        """
        if strain_genomes is not None and strain_sequences is not None:
            raise EmbeddingCacheError(
                "populate(): pass either strain_genomes+annotations OR strain_sequences, not both"
            )
        if strain_genomes is None and strain_sequences is None:
            raise EmbeddingCacheError(
                "populate(): must provide strain_genomes+annotations OR strain_sequences"
            )

        if strain_genomes is not None:
            # Plan-contract path: extract CDS sequences from FASTA + annotations
            if annotations is None:
                raise EmbeddingCacheError(
                    "populate(): strain_genomes requires matching annotations dict"
                )
            from dna_decode.data.annotations import extract_cds_sequences

            strain_sequences = {}
            for strain_id, genome_path in strain_genomes.items():
                ann = annotations.get(strain_id)
                if ann is None:
                    raise EmbeddingCacheError(
                        f"populate(): no annotation for strain {strain_id!r}"
                    )
                strain_sequences[strain_id] = extract_cds_sequences(genome_path, ann)

        # Open HDF5 ONCE for the full populate() batch (Wave 2.5 C6 fix)
        written_per_strain: dict[str, int] = {}
        with self._h5py.File(self.path, "a") as f:
            for strain_id, gene_map in strain_sequences.items():
                count = 0
                for gene_id, sequence in gene_map.items():
                    ds_path = self._dataset_path(strain_id, gene_id)
                    if skip_existing and ds_path in f:
                        continue
                    emb = model.embed_batch([sequence])[0]
                    if emb.ndim != 1 or emb.shape[0] != self.embedding_dim:
                        raise EmbeddingCacheError(
                            f"embedding shape {emb.shape} mismatch with cache "
                            f"embedding_dim={self.embedding_dim}"
                        )
                    if ds_path in f:
                        del f[ds_path]
                    f.create_dataset(ds_path, data=emb.astype(np.float32))
                    count += 1
                written_per_strain[strain_id] = count
                if progress_callback:
                    progress_callback(strain_id, count, len(gene_map))
        return written_per_strain
