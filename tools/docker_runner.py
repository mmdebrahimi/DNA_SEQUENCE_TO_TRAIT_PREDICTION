"""Stage 2 bioinformatics-tool runner via Docker Desktop.

Single Python entry point for Mash / Bakta / AMRFinderPlus invocation on a
Windows host with no native binaries. Replaces the `.sh` wrappers proposed
in an earlier draft of the install plan (subprocess.run does not launch
.sh on Windows without an explicit bash.exe shim, and `-v C:\\...` Windows
paths don't translate inside Linux containers).

Caller passes a `mounts` dict of host paths → container paths; the runner
normalizes Windows paths to forward-slash form (Docker Desktop on Windows
accepts `C:/Users/...:/container/path` natively) and builds the
`docker run --rm -v <host>:<container> <image> <args>` invocation.

Inside-container arg paths must already reference container-side mount
points (e.g., `--database /db/latest` after `mounts={"C:/amrfinder_db": "/db"}`).
The runner does NOT auto-translate inside-container args.
"""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


class DockerRunnerError(RuntimeError):
    """Docker invocation failed or returned a non-zero exit code."""


def _normalize_host_path(host: str) -> str:
    return str(Path(host).resolve()).replace("\\", "/")


def run(
    image: str,
    args: list[str],
    *,
    mounts: dict[str, str] | None = None,
    env: dict[str, str] | None = None,
    capture_output: bool = True,
    check: bool = True,
    timeout: float | None = None,
) -> subprocess.CompletedProcess:
    """Run `docker run --rm <image> <args>` with the given bind mounts.

    Args:
        image: Fully-qualified Docker image tag (e.g.,
            ``quay.io/biocontainers/mash:2.3--he348c14_4``).
        args: CLI arguments passed inside the container, after the image.
        mounts: ``host_path -> container_path`` bind-mount dict. Host paths
            are normalized to forward-slash form. Container paths must be
            Linux-style absolute paths. Append ``:ro`` to the container
            value for read-only mounts (e.g., ``{"D:/refseq": "/refseq:ro"}``
            -> ``-v D:/refseq:/refseq:ro``).
        env: Environment variables to inject via ``-e KEY=VALUE``.
        capture_output: Capture stdout+stderr (True) or stream to console.
        check: If True, raise DockerRunnerError on non-zero exit code.
        timeout: Max wall-clock seconds before terminating the container.
            None (default) means no timeout — caller takes responsibility
            for hung-daemon / paused-Docker-Desktop / stuck-pull cases.

    Returns:
        subprocess.CompletedProcess (stdout/stderr decoded as text).

    Raises:
        DockerRunnerError: when docker exits non-zero and check=True, or
            the `docker` binary is not on PATH, or the call exceeds
            `timeout` seconds.
    """
    cmd: list[str] = ["docker", "run", "--rm"]

    if mounts:
        for host, container in mounts.items():
            cmd.extend(["-v", f"{_normalize_host_path(host)}:{container}"])

    if env:
        for k, v in env.items():
            cmd.extend(["-e", f"{k}={v}"])

    cmd.append(image)
    cmd.extend(args)

    log.debug("docker_runner: %s", " ".join(cmd))

    try:
        proc = subprocess.run(
            cmd,
            capture_output=capture_output,
            text=True,
            check=False,
            timeout=timeout,
        )
    except FileNotFoundError as e:
        raise DockerRunnerError(
            f"`docker` not found on PATH — is Docker Desktop installed and running?\n"
            f"cmd: {' '.join(cmd)}"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise DockerRunnerError(
            f"docker invocation exceeded timeout={timeout}s\n"
            f"cmd: {' '.join(cmd)}"
        ) from e

    if check and proc.returncode != 0:
        raise DockerRunnerError(
            f"docker exited {proc.returncode}\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stderr: {proc.stderr}"
        )

    return proc
