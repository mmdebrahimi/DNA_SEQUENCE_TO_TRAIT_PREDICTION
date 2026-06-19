"""Thin Bakta (db-light) annotation runner (Step 1).

Adapts ``scripts/pathotype_laptop_pipeline.bakta_annotate`` to a generic
single-FASTA signature for the genome-map package. Pinned image, db-light,
D:-cache-friendly, idempotent (skip-existing). The Bakta image entrypoint is
already ``bakta`` so the args must NOT repeat it.

Bakta annotation was deferred / never smoke-tested on this host (CPU-heavy;
Docker-mount-corruption history) — the Step-2 tool-surface manifest is where
it is first proven. A wedge surfaces as ``BAKTA_ANNOTATION_BLOCKED`` upstream,
never a fake GFF.
"""
from __future__ import annotations

import os
from pathlib import Path

from tools.docker_runner import run as docker_run

# Pinned Bakta image (matches CLAUDE.md Stage 2 toolchain + pathotype pipeline).
BAKTA_IMAGE = "oschwengers/bakta:v1.11.4"

# Machine-agnostic default DB parent dir (contains the ``db-light`` subdir).
# Resolves on this laptop AND the workhorse; override with DNA_DECODE_BAKTA_DB.
BAKTA_DB = str(
    Path(os.environ.get("DNA_DECODE_BAKTA_DB", str(Path.home() / "dna_decode_stage2" / "bakta_db")))
)


def build_bakta_args(fasta_name: str, prefix: str, threads: int) -> list[str]:
    """Build the container-side Bakta CLI args (no image, no `bakta` prefix).

    Pure + testable: the mounts put the FASTA at ``/data/<fasta_name>``, the DB
    parent at ``/db`` (so ``/db/db-light``), and the output dir at ``/out``.
    """
    return [
        "--db", "/db/db-light",
        "--output", "/out",
        "--prefix", prefix,
        "--skip-plot",
        "--force",
        "--threads", str(threads),
        f"/data/{fasta_name}",
    ]


def run_bakta(
    fasta: Path | str,
    out_dir: Path | str,
    *,
    prefix: str | None = None,
    bakta_db: str | None = None,
    threads: int = 4,
    timeout: float = 2400,
) -> Path:
    """Annotate one genome FASTA with Bakta (Docker). Returns the GFF3 path.

    Idempotent: returns the existing ``<out_dir>/<prefix>.gff3`` without
    re-running if present. ``prefix`` defaults to the FASTA stem.
    """
    fasta = Path(fasta)
    out_dir = Path(out_dir)
    prefix = prefix or fasta.stem
    bakta_db = bakta_db or BAKTA_DB
    gff = out_dir / f"{prefix}.gff3"
    if gff.exists():
        return gff
    out_dir.mkdir(parents=True, exist_ok=True)
    docker_run(
        image=BAKTA_IMAGE,
        # image entrypoint is already `bakta` -> do NOT prepend "bakta".
        args=build_bakta_args(fasta.name, prefix, threads),
        mounts={
            str(fasta.parent): "/data:ro",
            bakta_db: "/db:ro",
            str(out_dir): "/out",
        },
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    if not gff.exists():
        raise RuntimeError(f"Bakta produced no GFF3 for {fasta.name} (prefix={prefix})")
    return gff
