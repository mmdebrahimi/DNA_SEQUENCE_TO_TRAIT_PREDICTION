"""Tests for scripts/populate_cache.py — standalone embedding-cache driver."""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest


def _write_minimal_fasta(path: Path, gene_seqs: dict[str, str]) -> None:
    """Write a multi-record FASTA with one record per gene_id."""
    with open(path, "w", encoding="utf-8") as f:
        for gid, seq in gene_seqs.items():
            f.write(f">{gid}\n{seq}\n")


def _write_minimal_gff3(path: Path, gene_records: list[tuple[str, str, int, int]]) -> None:
    """Write a tiny GFF3 with CDS rows.

    gene_records: list of (gene_id, locus_tag, start, end). Assumes a single
    contig 'chrom1'.
    """
    lines = ["##gff-version 3"]
    for gid, locus, start, end in gene_records:
        attrs = f"ID={gid};locus_tag={locus};gene={gid}"
        lines.append(
            f"chrom1\tphase2-test\tCDS\t{start}\t{end}\t.\t+\t0\t{attrs}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def fake_refseq_cache(tmp_path: Path) -> Path:
    """Build a minimal RefSeq cache with two accessions, each with fna + gff3."""
    root = tmp_path / "refseq_cache"

    for acc, contig_seq in [
        ("GCF_AAA.1", "ATGAAAGGGCCCTAA" * 20),
        ("GCF_BBB.1", "ATGCCCAAAGGGTAA" * 20),
    ]:
        acc_dir = root / acc
        acc_dir.mkdir(parents=True)
        # Write a single-contig FASTA covering the gene positions used in the gff3
        (acc_dir / "genome.fna").write_text(
            f">chrom1\n{contig_seq}\n", encoding="utf-8"
        )
        _write_minimal_gff3(
            acc_dir / "annotations.gff3",
            [("g1", "TAG_001", 1, 30), ("g2", "TAG_002", 31, 60)],
        )
        (acc_dir / ".complete").touch()

    return root


@pytest.fixture
def fake_cohort_parquet(tmp_path: Path) -> Path:
    """Build a tiny cohort.parquet with two strains pointing at the fake refseq accessions."""
    from dna_decode.data.cohort import (
        CandidateStrain,
        StrainCohort,
        save_cohort,
    )

    strains = [
        CandidateStrain(
            strain_id="synth_001",
            assembly_accession="GCF_AAA.1",
            mlst="ST1",
            contig_count=10,
            n50=200_000,
            ast_labels={"ciprofloxacin": 1},
        ),
        CandidateStrain(
            strain_id="synth_002",
            assembly_accession="GCF_BBB.1",
            mlst="ST2",
            contig_count=10,
            n50=200_000,
            ast_labels={"ciprofloxacin": 0},
        ),
    ]
    cohort = StrainCohort(
        strains=strains,
        per_drug_strain_ids={"ciprofloxacin": ["synth_001", "synth_002"]},
        three_drug_intersection=["synth_001", "synth_002"],
    )
    parquet_path = tmp_path / "cohort.parquet"
    save_cohort(cohort, parquet_path)
    return parquet_path


# ---- build_parser ----


def test_build_parser_has_required_args():
    from scripts.populate_cache import build_parser

    parser = build_parser()
    args = parser.parse_args(
        [
            "--cohort", "c.parquet",
            "--cache", "cache.h5",
            "--refseq-cache", "rc",
            "--allow-mock",
            "--model", "mock",
        ]
    )
    assert args.model == "mock"
    assert args.allow_mock is True


def test_parser_requires_cohort_cache_refseq():
    """Argparse enforces the three required flags."""
    from scripts.populate_cache import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])  # all three required args missing


# ---- main: gate behavior ----


def test_mock_without_allow_mock_exits_2(tmp_path: Path):
    """--model mock without --allow-mock fails fast."""
    from scripts.populate_cache import main

    exit_code = main(
        [
            "--cohort", str(tmp_path / "fake.parquet"),
            "--cache", str(tmp_path / "fake.h5"),
            "--refseq-cache", str(tmp_path),
            "--model", "mock",
        ]
    )
    assert exit_code == 2


def test_missing_cohort_exits_2(tmp_path: Path):
    from scripts.populate_cache import main

    exit_code = main(
        [
            "--cohort", str(tmp_path / "missing.parquet"),
            "--cache", str(tmp_path / "cache.h5"),
            "--refseq-cache", str(tmp_path),
            "--model", "mock",
            "--allow-mock",
        ]
    )
    assert exit_code == 2


def test_missing_refseq_cache_exits_2(tmp_path: Path, fake_cohort_parquet: Path):
    from scripts.populate_cache import main

    exit_code = main(
        [
            "--cohort", str(fake_cohort_parquet),
            "--cache", str(tmp_path / "cache.h5"),
            "--refseq-cache", str(tmp_path / "missing_refseq"),
            "--model", "mock",
            "--allow-mock",
        ]
    )
    assert exit_code == 2


# ---- resolve_strain_assets ----


def test_resolve_strain_assets_happy_path(
    fake_cohort_parquet: Path, fake_refseq_cache: Path
):
    from dna_decode.data.cohort import load_cohort
    from scripts.populate_cache import resolve_strain_assets

    cohort = load_cohort(fake_cohort_parquet)
    genomes, anns, skipped = resolve_strain_assets(cohort, fake_refseq_cache)

    assert set(genomes.keys()) == {"synth_001", "synth_002"}
    assert set(anns.keys()) == {"synth_001", "synth_002"}
    assert skipped == []


def test_resolve_strain_assets_skips_missing_fna(
    fake_cohort_parquet: Path, fake_refseq_cache: Path
):
    """Remove one strain's genome.fna → that strain is skipped, the other resolved."""
    (fake_refseq_cache / "GCF_AAA.1" / "genome.fna").unlink()

    from dna_decode.data.cohort import load_cohort
    from scripts.populate_cache import resolve_strain_assets

    cohort = load_cohort(fake_cohort_parquet)
    genomes, anns, skipped = resolve_strain_assets(cohort, fake_refseq_cache)

    assert "synth_001" not in genomes
    assert "synth_002" in genomes
    assert len(skipped) == 1
    skipped_ids = {sid for sid, _ in skipped}
    assert "synth_001" in skipped_ids


def test_resolve_strain_assets_skips_missing_accession(tmp_path: Path, fake_refseq_cache: Path):
    """A strain with empty assembly_accession is skipped with a clear reason."""
    from dna_decode.data.cohort import (
        CandidateStrain,
        StrainCohort,
        save_cohort,
        load_cohort,
    )
    from scripts.populate_cache import resolve_strain_assets

    strains = [
        CandidateStrain(
            strain_id="no_acc",
            assembly_accession="",
            mlst="ST_X",
            contig_count=10,
            n50=200_000,
            ast_labels={"ciprofloxacin": 1},
        ),
    ]
    cohort = StrainCohort(
        strains=strains,
        per_drug_strain_ids={"ciprofloxacin": ["no_acc"]},
        three_drug_intersection=["no_acc"],
    )
    path = tmp_path / "no_acc_cohort.parquet"
    save_cohort(cohort, path)

    genomes, anns, skipped = resolve_strain_assets(load_cohort(path), fake_refseq_cache)

    assert genomes == {}
    assert anns == {}
    assert len(skipped) == 1
    assert "missing assembly_accession" in skipped[0][1]


# ---- end-to-end mock populate ----


def test_main_mock_populates_cache(
    tmp_path: Path,
    fake_cohort_parquet: Path,
    fake_refseq_cache: Path,
):
    """Full end-to-end: mock model + fake cohort + fake refseq → HDF5 with embeddings."""
    from scripts.populate_cache import main

    cache_path = tmp_path / "populated.h5"
    exit_code = main(
        [
            "--cohort", str(fake_cohort_parquet),
            "--cache", str(cache_path),
            "--refseq-cache", str(fake_refseq_cache),
            "--model", "mock",
            "--allow-mock",
        ]
    )
    assert exit_code == 0
    assert cache_path.exists()

    with h5py.File(cache_path, "r") as f:
        # MockFoundationModel default embedding_dim = 128
        # Each strain has 2 CDS in the fake gff3 → 2 datasets per strain
        for sid in ("synth_001", "synth_002"):
            grp = f[f"strains/{sid}"]
            assert len(grp) == 2  # g1, g2
            for gid in grp:
                arr = np.array(grp[gid])
                assert arr.dtype == np.float32
                assert arr.shape == (128,)
