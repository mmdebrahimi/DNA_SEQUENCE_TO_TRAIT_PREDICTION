"""Step 10 — Mash/ANI distance clustering for phylogenetic CV.

Higher-resolution than MLST. Phase 1 primary phylogeny control per
post-tech-plan brainstorm M2.

External binary dependency: Mash CLI. Two execution paths:
    (a) Native binary on PATH (Linux / WSL2 / Windows binary release)
    (b) Containerized via `tools/docker_runner.py` (Stage 2 preferred on
        Windows hosts that lack a native build)

Batched-call discipline (Stage 2 prep, 2026-05-15): the all-pairs distance
matrix is computed in 2 mash invocations — `mash sketch -o sketch.msh <all
fastas...>` then `mash dist sketch.msh sketch.msh` — instead of the prior
N*(N-1)/2 invocations. At N=147 that is 2 calls vs 10,731. With docker
spin-up overhead the prior pattern would cost hours of preflight time per
LOMO-clade-out CV pass; the batched call is seconds.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


DEFAULT_MASH_BINARY = "mash"
DEFAULT_MASH_DOCKER_IMAGE = "quay.io/biocontainers/mash:2.3--hb105d93_10"
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


def _run_mash_native(args: list[str], mash_binary: str = DEFAULT_MASH_BINARY) -> str:
    """Invoke a native `mash` binary; return stdout."""
    cmd = [mash_binary] + args
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=False, encoding="utf-8"
        )
    except FileNotFoundError as e:
        raise MashError(
            f"`{mash_binary}` binary not found on PATH. Install Mash: "
            "https://github.com/marbl/Mash/releases (Windows) or "
            "`apt install mash` (Linux/WSL), or pass use_docker=True."
        ) from e
    if result.returncode != 0:
        raise MashError(
            f"mash {' '.join(args)} failed (exit {result.returncode}): "
            f"{result.stderr.strip()[:200]}"
        )
    return result.stdout


def _run_mash_docker(
    args: list[str],
    mount_host_dir: Path,
    image: str = DEFAULT_MASH_DOCKER_IMAGE,
) -> str:
    """Invoke containerized mash with the working dir mounted at /data."""
    from tools.docker_runner import DockerRunnerError, run as docker_run

    try:
        proc = docker_run(
            image,
            ["mash"] + args,
            mounts={str(mount_host_dir): "/data"},
            capture_output=True,
        )
    except DockerRunnerError as e:
        raise MashError(f"containerized mash failed: {e}") from e
    return proc.stdout


def compute_mash_distances(
    strain_genomes: dict[str, Path],
    *,
    mash_binary: str = DEFAULT_MASH_BINARY,
    use_docker: bool = False,
    docker_image: str = DEFAULT_MASH_DOCKER_IMAGE,
) -> DistanceMatrix:
    """Run `mash sketch` + `mash dist` to produce pairwise ANI distances.

    Uses ONE `mash sketch` call across all genomes, then ONE `mash dist`
    call on the sketch against itself, producing the full N*N distance
    matrix. At N=147 this is 2 invocations instead of 10,731.

    Args:
        strain_genomes: mapping strain_id -> path to genome FASTA.
        mash_binary: name or path of the mash executable (native path).
        use_docker: route through `tools.docker_runner.run` using
            `docker_image`. Required on Windows hosts without a native mash.
        docker_image: pinned biocontainers image tag.

    Returns:
        DistanceMatrix where matrix[i, j] is the Mash distance between
        strain_ids[i] and strain_ids[j].
    """
    strain_ids = list(strain_genomes.keys())
    n = len(strain_ids)
    if n == 0:
        return DistanceMatrix(strain_ids=[], matrix=np.zeros((0, 0)))

    with tempfile.TemporaryDirectory(prefix="mash_batch_") as tmp:
        work = Path(tmp)
        sid_to_arg: dict[str, str] = {}
        for sid, fna in strain_genomes.items():
            staged_name = f"{sid}.fna"
            staged = work / staged_name
            shutil.copy(fna, staged)
            sid_to_arg[sid] = staged_name if use_docker else str(staged)

        if use_docker:
            sketch_args = ["sketch", "-o", "/data/sketch"] + [
                f"/data/{sid_to_arg[sid]}" for sid in strain_ids
            ]
            _run_mash_docker(sketch_args, work, image=docker_image)

            dist_stdout = _run_mash_docker(
                ["dist", "/data/sketch.msh", "/data/sketch.msh"],
                work,
                image=docker_image,
            )
        else:
            sketch_path = work / "sketch"
            sketch_args = ["sketch", "-o", str(sketch_path)] + [
                sid_to_arg[sid] for sid in strain_ids
            ]
            _run_mash_native(sketch_args, mash_binary=mash_binary)

            msh = str(sketch_path) + ".msh"
            dist_stdout = _run_mash_native(
                ["dist", msh, msh], mash_binary=mash_binary
            )

    return _parse_all_pairs_dist(dist_stdout, strain_ids, sid_to_arg)


def _parse_all_pairs_dist(
    stdout: str,
    strain_ids: list[str],
    sid_to_arg: dict[str, str],
) -> DistanceMatrix:
    """Parse mash dist all-pairs output into a (n, n) DistanceMatrix.

    `mash dist sketch.msh sketch.msh` writes one line per (ref, query)
    pair:
        <ref-path>\t<query-path>\t<distance>\t<pvalue>\t<shared/total>
    Order is ref-major. Self-pairs appear with distance 0.
    """
    if not stdout.strip():
        raise MashError("empty mash dist stdout (expected all-pairs output)")

    n = len(strain_ids)
    matrix = np.zeros((n, n), dtype=np.float32)

    # Map any staged-path tail back to the strain_id index.
    name_to_idx: dict[str, int] = {}
    for i, sid in enumerate(strain_ids):
        # strain_ids[i]'s staged path ends with sid_to_arg[sid]; both /data/<x>.fna
        # (docker) and absolute (native) tails end with "<x>.fna".
        name_to_idx[Path(sid_to_arg[sid]).name] = i

    for line in stdout.strip().splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            raise MashError(f"unexpected mash dist line: {line!r}")
        ref_name = Path(fields[0]).name
        query_name = Path(fields[1]).name
        try:
            d = float(fields[2])
        except ValueError as e:
            raise MashError(f"non-numeric distance in: {line!r}") from e
        if ref_name not in name_to_idx or query_name not in name_to_idx:
            raise MashError(
                f"mash dist row references unknown staged file: ref={ref_name!r} "
                f"query={query_name!r}"
            )
        i = name_to_idx[ref_name]
        j = name_to_idx[query_name]
        matrix[i, j] = d

    return DistanceMatrix(strain_ids=strain_ids, matrix=matrix)


def _parse_mash_dist_line(stdout: str) -> float:
    """Parse a single `mash dist` output line: ref query distance pvalue shared/total.

    Retained for backwards compatibility with prior single-pair callers
    + existing tests. New code should use `compute_mash_distances` which
    parses the all-pairs output internally.
    """
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
