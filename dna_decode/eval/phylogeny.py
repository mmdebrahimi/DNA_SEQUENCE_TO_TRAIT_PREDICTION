"""Step 10 — Mash/ANI distance clustering for phylogenetic CV.

Higher-resolution than MLST. Phase 1 primary phylogeny control per
post-tech-plan brainstorm M2.

External binary dependency: `mash` CLI (https://github.com/marbl/Mash).
Tests mock subprocess; real-data runs require system-installed binary.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np


DEFAULT_MASH_BINARY = "mash"
DEFAULT_ANI_THRESHOLD = 0.02  # ~98% genome identity at sub-species level


class MashError(Exception):
    """Mash CLI invocation failure."""


@dataclass
class DistanceMatrix:
    """Pairwise distance matrix indexed by strain_id."""

    strain_ids: list[str]
    matrix: np.ndarray  # (n, n), symmetric, diagonal=0

    def __post_init__(self) -> None:
        n = len(self.strain_ids)
        if self.matrix.shape != (n, n):
            raise ValueError(
                f"matrix shape {self.matrix.shape} != ({n}, {n})"
            )


def _run_mash(args: list[str], mash_binary: str = DEFAULT_MASH_BINARY) -> str:
    """Invoke `mash` CLI; return stdout."""
    cmd = [mash_binary] + args
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=False, encoding="utf-8"
        )
    except FileNotFoundError as e:
        raise MashError(
            f"`{mash_binary}` binary not found on PATH. Install Mash: "
            "https://github.com/marbl/Mash/releases (Windows) or "
            "`apt install mash` (Linux/WSL)."
        ) from e
    if result.returncode != 0:
        raise MashError(
            f"mash {' '.join(args)} failed (exit {result.returncode}): {result.stderr.strip()[:200]}"
        )
    return result.stdout


def compute_mash_distances(
    strain_genomes: dict[str, Path],
    mash_binary: str = DEFAULT_MASH_BINARY,
) -> DistanceMatrix:
    """Run `mash sketch` + `mash dist` to produce pairwise ANI distances.

    Args:
        strain_genomes: mapping strain_id -> path to genome FASTA.
        mash_binary: name or path of the mash executable.

    Returns:
        DistanceMatrix where matrix[i, j] is the Mash distance (~ANI mismatch)
        between strain_ids[i] and strain_ids[j].
    """
    strain_ids = list(strain_genomes.keys())
    n = len(strain_ids)
    if n == 0:
        return DistanceMatrix(strain_ids=[], matrix=np.zeros((0, 0)))

    matrix = np.zeros((n, n), dtype=np.float32)
    for i, sid_i in enumerate(strain_ids):
        for j, sid_j in enumerate(strain_ids):
            if i >= j:
                continue
            stdout = _run_mash(
                ["dist", str(strain_genomes[sid_i]), str(strain_genomes[sid_j])],
                mash_binary=mash_binary,
            )
            distance = _parse_mash_dist_line(stdout)
            matrix[i, j] = distance
            matrix[j, i] = distance

    return DistanceMatrix(strain_ids=strain_ids, matrix=matrix)


def _parse_mash_dist_line(stdout: str) -> float:
    """Parse `mash dist` output: ref  query  distance  pvalue  shared/total."""
    if not stdout.strip():
        raise MashError("empty mash dist stdout")
    line = stdout.strip().splitlines()[0]
    fields = line.split("\t")
    if len(fields) < 3:
        raise MashError(f"unexpected mash dist line: {line!r}")
    try:
        return float(fields[2])
    except ValueError as e:
        raise MashError(f"non-numeric distance in: {line!r}") from e


def cluster_by_ani(
    distance_matrix: DistanceMatrix,
    threshold: float = DEFAULT_ANI_THRESHOLD,
) -> dict[str, int]:
    """Hierarchical clustering on the distance matrix.

    Returns strain_id -> cluster_id mapping. Cluster IDs start at 0.

    Implementation: simple union-find on pairs with distance < threshold.
    Avoids the scipy dependency for Phase 1; produces the same single-linkage
    clusters at the chosen threshold.
    """
    n = len(distance_matrix.strain_ids)
    if n == 0:
        return {}

    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix.matrix[i, j] < threshold:
                union(i, j)

    # Compact cluster IDs to 0..k-1
    root_to_id: dict[int, int] = {}
    out: dict[str, int] = {}
    next_id = 0
    for i, sid in enumerate(distance_matrix.strain_ids):
        root = find(i)
        if root not in root_to_id:
            root_to_id[root] = next_id
            next_id += 1
        out[sid] = root_to_id[root]
    return out
