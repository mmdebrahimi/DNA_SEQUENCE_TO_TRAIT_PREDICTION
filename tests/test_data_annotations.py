"""Tests for Step 3 — Genome annotation parser."""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.data import annotations as ann


SAMPLE_GFF3 = """\
##gff-version 3
seq1\tBakta\tCDS\t1\t99\t.\t+\t0\tID=g1;locus_tag=TAG_001;product=hypothetical protein
seq1\tBakta\tCDS\t200\t299\t.\t-\t0\tID=g2;locus_tag=TAG_002;product=DNA gyrase subunit A
seq1\tBakta\tgene\t1\t99\t.\t+\t.\tID=gene1;gene=gyrA
"""

SAMPLE_FASTA = """\
>seq1
ATGCGTAAACCCTTTGGGCATCATGAATTCGCGCGAGGGCCCAAATTTCCCGGGAGGCGCTAGAGCAGTAGGCATATTGCTCAGCATCGCATGATAAAGTC\
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""


@pytest.fixture
def gff3_file(tmp_path: Path) -> Path:
    p = tmp_path / "sample.gff3"
    p.write_text(SAMPLE_GFF3, encoding="utf-8")
    return p


@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "sample.fna"
    p.write_text(SAMPLE_FASTA, encoding="utf-8")
    return p


# ---- attribute parsing ----


def test_parse_gff3_attrs_basic():
    parsed = ann._parse_gff3_attrs("ID=g1;locus_tag=TAG_001;product=hypothetical")
    assert parsed == {"ID": "g1", "locus_tag": "TAG_001", "product": "hypothetical"}


def test_parse_gff3_attrs_empty():
    assert ann._parse_gff3_attrs("") == {}


def test_parse_gff3_attrs_skips_malformed():
    parsed = ann._parse_gff3_attrs("ID=g1;invalid;locus_tag=TAG")
    assert parsed == {"ID": "g1", "locus_tag": "TAG"}


# ---- GFF3 parsing ----


def test_parse_gff3_produces_stable_columns(gff3_file: Path):
    df = ann.parse_gff3(gff3_file)
    assert list(df.columns) == list(ann.ANNOTATION_COLUMNS)


def test_parse_gff3_row_count_matches_data_lines(gff3_file: Path):
    df = ann.parse_gff3(gff3_file)
    assert len(df) == 3


def test_parse_gff3_extracts_gene_id_and_locus_tag(gff3_file: Path):
    df = ann.parse_gff3(gff3_file)
    cds = df[df["type"] == "CDS"]
    assert "g1" in cds["gene_id"].values
    assert "TAG_001" in cds["locus_tag"].values


def test_parse_gff3_preserves_strand(gff3_file: Path):
    df = ann.parse_gff3(gff3_file)
    strands = set(df[df["type"] == "CDS"]["strand"])
    assert strands == {"+", "-"}


def test_parse_gff3_raises_on_wrong_field_count(tmp_path: Path):
    bad = tmp_path / "bad.gff3"
    bad.write_text("seq1\tBakta\tCDS\t1\t99\t.\t+\n", encoding="utf-8")  # only 7 fields
    with pytest.raises(ann.AnnotationParseError, match="9"):
        ann.parse_gff3(bad)


def test_parse_gff3_raises_on_non_integer_start(tmp_path: Path):
    bad = tmp_path / "bad.gff3"
    bad.write_text("seq1\tBakta\tCDS\tNOT_AN_INT\t99\t.\t+\t0\tID=g\n", encoding="utf-8")
    with pytest.raises(ann.AnnotationParseError, match="non-integer"):
        ann.parse_gff3(bad)


def test_parse_gff3_empty_file(tmp_path: Path):
    empty = tmp_path / "empty.gff3"
    empty.write_text("##gff-version 3\n", encoding="utf-8")  # header only
    df = ann.parse_gff3(empty)
    assert len(df) == 0
    assert list(df.columns) == list(ann.ANNOTATION_COLUMNS)


# ---- revcomp ----


def test_revcomp_basic():
    assert ann._revcomp("ATGC") == "GCAT"


def test_revcomp_mixed_case_normalizes_handled():
    # Lowercase tolerated; downstream extract_cds_sequences uppercases via _load_genome_dict
    assert ann._revcomp("atgc") == "gcat"


def test_revcomp_with_n_preserved():
    assert ann._revcomp("ATNC") == "GNAT"


# ---- CDS extraction ----


def test_extract_cds_sequences_plus_strand(gff3_file: Path, fasta_file: Path):
    df = ann.parse_gff3(gff3_file)
    seqs = ann.extract_cds_sequences(fasta_file, df)
    # g1 is CDS on seq1 from 1-99 on + strand → first 99 bases
    assert "g1" in seqs
    assert len(seqs["g1"]) == 99


def test_extract_cds_sequences_minus_strand_revcomps(gff3_file: Path, fasta_file: Path):
    df = ann.parse_gff3(gff3_file)
    seqs = ann.extract_cds_sequences(fasta_file, df)
    # g2 is CDS on seq1 from 200-299 on - strand → revcomp of bases 200-299
    assert "g2" in seqs
    assert len(seqs["g2"]) == 100


def test_extract_cds_sequences_skips_non_cds(gff3_file: Path, fasta_file: Path):
    df = ann.parse_gff3(gff3_file)
    seqs = ann.extract_cds_sequences(fasta_file, df)
    # 'gene' rows are not CDS → not extracted
    assert "gene1" not in seqs


# ---- intergenic extraction ----


def test_extract_intergenic_regions_finds_gap(gff3_file: Path, fasta_file: Path):
    df = ann.parse_gff3(gff3_file)
    # CDS g1 ends at 99, g2 starts at 200 → gap of 100 bp
    intergenic = ann.extract_intergenic_regions(fasta_file, df, min_length=30)
    assert "g1__g2" in intergenic
    assert len(intergenic["g1__g2"]) == 100


def test_extract_intergenic_regions_respects_min_length(gff3_file: Path, fasta_file: Path):
    df = ann.parse_gff3(gff3_file)
    intergenic = ann.extract_intergenic_regions(fasta_file, df, min_length=200)
    assert intergenic == {}
