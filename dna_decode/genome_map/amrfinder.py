"""Thin AMRFinderPlus runner (Step 1, brainstorm catch C1).

AMRFinder must be RUN — the overlay (Step 4) parses the raw ``main.tsv`` for
coords/protein-ids, so a precomputed cell verdict is not enough. This reuses
``scripts.drug_mechanism_audit._run_amrfinder`` (pinned image + cross-drive DB
symlink handling) for the organism path, and adds a no-``-O`` path for a
genuinely unfamiliar genome (the 3rd spike prototype) where organism-specific
point-mutation detection is neither available nor expected.

The organism is an EXPLICIT input — v1 does NOT auto-detect it from Bakta
(an explicit ``-O`` is what makes gyrA/parC QRDR calls organism-correct).
"""
from __future__ import annotations

from pathlib import Path

from tools.docker_runner import run as docker_run


def run_amrfinder(
    fasta: Path | str,
    out_dir: Path | str,
    organism: str | None = "Escherichia",
    *,
    timeout_sec: float = 600,
) -> tuple[Path, Path]:
    """Run AMRFinderPlus on one FASTA -> (main_tsv, mutations_tsv).

    ``organism`` selects AMRFinder's ``-O`` (organism-specific QRDR / point
    mutations). Pass ``None`` to run WITHOUT ``-O`` (a generic scan for an
    organism AMRFinder doesn't curate — the unfamiliar-genome prototype);
    acquired-gene determinants are still detected, only organism-specific
    point mutations are skipped.
    """
    fasta = Path(fasta)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Lazy import: scripts.drug_mechanism_audit pulls pandas + cohort deps; keep
    # the package import light and reuse its pinned image + DB-resolution logic.
    import scripts.drug_mechanism_audit as dma

    if organism:
        # Reuse the validated organism path verbatim (DB symlink-resolution etc.).
        return dma._run_amrfinder(fasta, out_dir, timeout_sec=timeout_sec, organism=organism)

    # No-organism path: mirror _run_amrfinder's DB resolution, drop -O.
    main_out = out_dir / "main.tsv"
    mut_out = out_dir / "mutations.tsv"
    db_latest = Path(dma.AMRFINDER_DB) / "latest"
    real_db = db_latest.resolve() if db_latest.exists() else db_latest
    docker_run(
        dma.AMRFINDER_IMAGE,
        [
            "amrfinder",
            "-n", f"/in/{fasta.name}",
            "--database", "/db/latest",
            "--mutation_all", "/out/mutations.tsv",
            "-o", "/out/main.tsv",
        ],
        mounts={
            str(fasta.parent): "/in:ro",
            str(real_db): "/db/latest:ro",
            str(out_dir): "/out",
        },
        timeout=timeout_sec,
    )
    return main_out, mut_out
