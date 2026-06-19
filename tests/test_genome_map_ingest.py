"""Step 1 tests — shared ##FASTA-safe GFF loader + Bakta + AMRFinder runner args.

No Docker / no live tools: the loader is tested on synthetic GFF text, and the
runner wrappers are tested for correct docker-arg construction with the docker
primitive mocked.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.genome_map import ingest


# A minimal valid Bakta-style GFF3 with the embedded ##FASTA block at the end.
_GFF_WITH_FASTA = """\
##gff-version 3
contig_1\tBakta\tgene\t100\t900\t.\t+\t.\tID=gene-1;gene=gyrA
contig_1\tBakta\tCDS\t100\t900\t.\t+\t0\tID=cds-1;Parent=gene-1;locus_tag=TAG_001;product=DNA gyrase subunit A
contig_1\tBakta\tCDS\t1000\t1500\t.\t-\t0\tID=cds-2;locus_tag=TAG_002;product=hypothetical protein
##FASTA
>contig_1
ACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGT
"""

_GFF_PLAIN = """\
##gff-version 3
contig_1\tRefSeq\tCDS\t100\t900\t.\t+\t0\tID=cds-1;gene=parC;locus_tag=b0001;product=topoisomerase IV subunit A
"""


# ---- strip_fasta_block ----


def test_strip_fasta_block_removes_sequence():
    out = ingest.strip_fasta_block(_GFF_WITH_FASTA)
    assert "##FASTA" not in out
    assert ">contig_1" not in out
    assert "DNA gyrase subunit A" in out


def test_strip_fasta_block_idempotent_on_plain():
    assert ingest.strip_fasta_block(_GFF_PLAIN) == _GFF_PLAIN


# ---- load_genome_gff ----


def test_load_genome_gff_parses_fasta_bearing(tmp_path: Path):
    gff = tmp_path / "g.gff3"
    gff.write_text(_GFF_WITH_FASTA, encoding="utf-8")
    table = ingest.load_genome_gff(gff)
    # 3 feature rows (gene + 2 CDS); the FASTA lines must NOT appear as rows.
    assert len(table) == 3
    products = set(table["product"])
    assert "DNA gyrase subunit A" in products
    assert "hypothetical protein" in products
    # the stripped temp file was written next to the source
    assert (tmp_path / "g.nofasta.gff3").exists()


def test_load_genome_gff_propagates_gene_symbol_to_cds(tmp_path: Path):
    # Bakta puts gene= on the parent gene row; the CDS inherits via Parent=.
    gff = tmp_path / "g.gff3"
    gff.write_text(_GFF_WITH_FASTA, encoding="utf-8")
    table = ingest.load_genome_gff(gff)
    cds1 = table[table["locus_tag"] == "TAG_001"].iloc[0]
    assert cds1["gene_symbol"] == "gyrA"


def test_load_genome_gff_plain_unaffected(tmp_path: Path):
    gff = tmp_path / "plain.gff3"
    gff.write_text(_GFF_PLAIN, encoding="utf-8")
    table = ingest.load_genome_gff(gff)
    assert len(table) == 1
    assert table.iloc[0]["gene_symbol"] == "parC"
    # no ##FASTA -> no temp file written
    assert not (tmp_path / "plain.nofasta.gff3").exists()


def test_load_genome_gff_raises_on_malformed(tmp_path: Path):
    from dna_decode.data.annotations import AnnotationParseError

    gff = tmp_path / "bad.gff3"
    gff.write_text("contig_1\tBakta\tCDS\tonly\tthree\n", encoding="utf-8")
    with pytest.raises(AnnotationParseError):
        ingest.load_genome_gff(gff)


# ---- annotate.run_bakta arg construction ----


def test_build_bakta_args_db_light_and_prefix():
    from dna_decode.genome_map import annotate

    args = annotate.build_bakta_args("genome.fna", "ACC1", threads=8)
    assert "/db/db-light" in args
    assert "ACC1" in args
    assert "/data/genome.fna" in args
    assert "--skip-plot" in args
    # entrypoint is already `bakta` — must not be repeated as an arg
    assert args[0] != "bakta"


def test_run_bakta_skips_when_gff_exists(tmp_path: Path, monkeypatch):
    from dna_decode.genome_map import annotate

    out = tmp_path / "out"
    out.mkdir()
    (out / "ACC1.gff3").write_text("##gff-version 3\n", encoding="utf-8")

    called = {"ran": False}

    def _boom(*a, **k):
        called["ran"] = True
        raise AssertionError("docker should not run when GFF exists")

    monkeypatch.setattr(annotate, "docker_run", _boom)
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    gff = annotate.run_bakta(fasta, out, prefix="ACC1")
    assert gff == out / "ACC1.gff3"
    assert called["ran"] is False


def test_run_bakta_builds_mounts_and_runs(tmp_path: Path, monkeypatch):
    from dna_decode.genome_map import annotate

    captured = {}

    def _fake(*, image, args, mounts, capture_output, check, timeout):
        captured["image"] = image
        captured["args"] = args
        captured["mounts"] = mounts
        # simulate Bakta producing the GFF
        out_container = [v for v in mounts.values() if v == "/out"]
        assert out_container
        gff = Path([h for h, c in mounts.items() if c == "/out"][0]) / "ACC1.gff3"
        gff.write_text("##gff-version 3\n", encoding="utf-8")
        return None

    monkeypatch.setattr(annotate, "docker_run", _fake)
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    out = tmp_path / "out"
    gff = annotate.run_bakta(fasta, out, prefix="ACC1")
    assert gff.exists()
    assert captured["image"] == annotate.BAKTA_IMAGE
    # the FASTA's parent is mounted read-only at /data
    assert any(c == "/data:ro" for c in captured["mounts"].values())


# ---- amrfinder.run_amrfinder ----


def test_run_amrfinder_organism_delegates(tmp_path: Path, monkeypatch):
    from dna_decode.genome_map import amrfinder
    import scripts.drug_mechanism_audit as dma

    calls = {}

    def _fake_runner(fasta, out_dir, timeout_sec=600, organism="Escherichia"):
        calls["organism"] = organism
        return out_dir / "main.tsv", out_dir / "mutations.tsv"

    monkeypatch.setattr(dma, "_run_amrfinder", _fake_runner)
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    main, mut = amrfinder.run_amrfinder(fasta, tmp_path / "out", organism="Klebsiella_pneumoniae")
    assert calls["organism"] == "Klebsiella_pneumoniae"
    assert main.name == "main.tsv"


def test_run_amrfinder_no_organism_drops_flag(tmp_path: Path, monkeypatch):
    from dna_decode.genome_map import amrfinder

    captured = {}

    def _fake_docker(image, args, *, mounts, timeout):
        captured["args"] = args
        captured["image"] = image
        return None

    monkeypatch.setattr(amrfinder, "docker_run", _fake_docker)
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    main, mut = amrfinder.run_amrfinder(fasta, tmp_path / "out", organism=None)
    assert "-O" not in captured["args"]
    assert "amrfinder" in captured["args"]
    assert main.name == "main.tsv"
