"""Mash invocation for the identification router -- the only part that touches Docker or the disk.

Split out from the CLI deliberately: Steps 4 (threshold analysis) and 6 (validation) both need to
query the reference, and neither should go through an argparse surface to do it.

The boundary is `tools.docker_runner.run`, which spawns no shell and is therefore immune to the MSYS
path-conversion trap that silently turns `-v /data` into `C:/Program Files/Git/data`.

ONE INVOCATION SHAPE, verified live against the pinned image:
    mash dist <reference.msh> <query.fna>   ->   [reference-ID, query-ID, distance, p-value, shared-hashes]

`mash screen` is NOT used. Its own `-h` documents it for read/contig MIXTURES; the query here is one
assembled genome, and `dist` is the right subcommand for that. FASTQ input is a named non-goal.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from dna_decode.identify.core import AbstainReason, Hit, parse_mash_dist

MASH_IMAGE = "quay.io/biocontainers/mash:2.3--hb105d93_10"
HOST_DRIVE = Path("D:/")
CONTAINER_MOUNT = "/d"
DEFAULT_REFERENCE = Path("D:/dna_decode_cache/identify/reference.msh")


class ReferenceUnavailable(RuntimeError):
    """The reference sketch or the Mash runner could not be used.

    A SEPARATE failure class from "no organism matched" on purpose. If these collapse, a wedged Docker
    Desktop or a missing sketch reads as a biological result -- and on this host a corrupted Docker D:
    mount is a KNOWN failure mode (`wsl --shutdown` is the recorded fix).
    """


@dataclass(frozen=True)
class QueryResult:
    hits: tuple[Hit, ...]
    stderr: str = ""


def _to_container_path(host: Path) -> str:
    """Rewrite a host path under D: for the container. Raises if it is not under the mount.

    Fail-closed: a query genome outside the mounted drive does not exist inside the container, so mash
    would report no hits and the router would abstain NO_HIT -- an operational fault wearing the
    costume of "not a supported organism".
    """
    r = host.resolve()
    try:
        rel = r.relative_to(HOST_DRIVE.resolve())
    except ValueError as exc:
        raise ReferenceUnavailable(
            f"{r} is not under the mounted drive {HOST_DRIVE}; it would be invisible inside the "
            f"container and the router would abstain as if the organism were unsupported") from exc
    return CONTAINER_MOUNT + "/" + rel.as_posix()


def query(genome_fasta: Path, reference: Path = DEFAULT_REFERENCE, *,
          timeout: float = 900.0) -> QueryResult:
    """Run one `mash dist` of the reference against this genome and parse the hits.

    Raises `ReferenceUnavailable` (never returns an empty result) when the sketch is missing, the
    query is missing, or the runner fails. The caller maps that to
    `AbstainReason.REFERENCE_UNAVAILABLE`, which is what keeps an infrastructure fault from being
    reported as a clean biological abstention.
    """
    from tools.docker_runner import DockerRunnerError
    from tools.docker_runner import run as docker_run

    ref, q = Path(reference), Path(genome_fasta)
    if not ref.exists():
        raise ReferenceUnavailable(
            f"reference sketch not found at {ref}. Build it with "
            f"`uv run python scripts/build_identify_reference.py`.")
    if not q.exists():
        raise ReferenceUnavailable(f"query genome not found: {q}")

    args = ["mash", "dist", _to_container_path(ref), _to_container_path(q)]
    try:
        res = docker_run(MASH_IMAGE, args, mounts={str(HOST_DRIVE): CONTAINER_MOUNT},
                         capture_output=True, check=True, timeout=timeout)
    except DockerRunnerError as exc:
        raise ReferenceUnavailable(f"mash dist failed: {exc}") from exc

    hits = parse_mash_dist(res.stdout or "")
    if not hits:
        # mash exited 0 but produced nothing parseable. Exit 0 is not evidence -- this run already
        # saw a sketch build reported as "exit code 0" while the script had exited 1.
        raise ReferenceUnavailable(
            f"mash dist returned no parseable rows against {ref.name}; treat exit 0 as unverified. "
            f"stderr: {(res.stderr or '')[:400]!r}")
    return QueryResult(tuple(hits), res.stderr or "")


def abstain_reason_for_unavailable() -> AbstainReason:
    """The single mapping from an infrastructure fault to its abstain reason."""
    return AbstainReason.REFERENCE_UNAVAILABLE
