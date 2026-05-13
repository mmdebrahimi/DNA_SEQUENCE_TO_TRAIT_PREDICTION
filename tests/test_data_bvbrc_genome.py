"""Tests for dna_decode/data/bvbrc_genome.py — BV-BRC genome-metadata adapter."""
from __future__ import annotations

import warnings
from pathlib import Path

import pytest

from dna_decode.data.bvbrc_genome import (
    BvBrcGenomeError,
    GENOME_METADATA_KEYS,
    load_bvbrc_genome_metadata,
)


# Real BV-BRC Genomes export header (Title Case + spaces). Subset of the full
# column list — only the columns the adapter consumes are required.
REAL_BVBRC_GENOME_CSV = (
    "Genome ID,Genome Name,Species,MLST,Assembly Accession,Contigs,"
    "Contig N50,Size,CheckM Completeness,Collection Year,Isolation Country,Plasmids\n"
    '"562.5691","Escherichia coli CVM N36113PS","Escherichia coli","MLST.ecoli.131",'
    '"GCF_001234567.1","48","85234","4900000","99.5","2014","USA","1"\n'
    '"562.144245","Escherichia coli ERR7221502","Escherichia coli","MLST.ecoli.10",'
    '"GCF_009876543.2","12","250000","5100000","99.8","2016","Japan","0"\n'
    '"562.1001","Escherichia coli foo","Escherichia coli","MLST.ecoli.69",'
    '"GCF_000111222.1","200","60000","4750000","98.2","2018","United Kingdom","2"\n'
    '"1280.3350","Staphylococcus aureus NRS100","Staphylococcus aureus","MLST.saureus.250",'
    '"GCF_000626615.1","1","2823087","2823087","99.5","2014","USA","0"\n'
)


@pytest.fixture
def bvbrc_genome_csv(tmp_path: Path) -> Path:
    p = tmp_path / "BVBRC_genome.csv"
    p.write_text(REAL_BVBRC_GENOME_CSV, encoding="utf-8")
    return p


# ---- happy path ----


def test_load_returns_dict_keyed_by_strain_id(bvbrc_genome_csv: Path):
    meta = load_bvbrc_genome_metadata(bvbrc_genome_csv)
    # 3 E. coli rows retained; Staph filtered by organism
    assert set(meta.keys()) == {"562.5691", "562.144245", "562.1001"}
    assert "1280.3350" not in meta  # Staph excluded


def test_load_populates_all_metadata_keys(bvbrc_genome_csv: Path):
    meta = load_bvbrc_genome_metadata(bvbrc_genome_csv)
    row = meta["562.5691"]
    assert set(row.keys()) == set(GENOME_METADATA_KEYS)


def test_load_parses_int_fields_correctly(bvbrc_genome_csv: Path):
    meta = load_bvbrc_genome_metadata(bvbrc_genome_csv)
    # contig_count + n50 + year are integers
    assert meta["562.5691"]["contig_count"] == 48
    assert meta["562.5691"]["n50"] == 85234
    assert meta["562.5691"]["year"] == 2014
    # MLST + country + assembly_accession are strings
    assert meta["562.5691"]["assembly_accession"] == "GCF_001234567.1"
    assert meta["562.5691"]["mlst"] == "MLST.ecoli.131"
    assert meta["562.5691"]["country"] == "USA"


def test_load_handles_tsv_via_sep_autodetect(tmp_path: Path):
    """sep=None should auto-detect tab as well as comma."""
    tsv = REAL_BVBRC_GENOME_CSV.replace(",", "\t").replace('"', "")
    p = tmp_path / "BVBRC_genome.tsv"
    p.write_text(tsv, encoding="utf-8")
    meta = load_bvbrc_genome_metadata(p)
    assert "562.5691" in meta


# ---- organism filter ----


def test_organism_filter_uses_species_column(bvbrc_genome_csv: Path):
    meta = load_bvbrc_genome_metadata(bvbrc_genome_csv, organism="Escherichia coli")
    assert "1280.3350" not in meta  # Staph excluded by Species column


def test_organism_filter_falls_back_to_genome_name(tmp_path: Path):
    """When Species column is absent, fall back to Genome Name substring match."""
    csv = (
        "Genome ID,Genome Name,MLST,Assembly Accession,Contigs,Contig N50,Collection Year,Isolation Country\n"
        '"562.A","Escherichia coli K-12","ST1","GCF_A","10","200000","2020","USA"\n'
        '"1280.A","Staphylococcus aureus","ST2","GCF_B","10","200000","2020","USA"\n'
    )
    p = tmp_path / "no_species.csv"
    p.write_text(csv, encoding="utf-8")
    meta = load_bvbrc_genome_metadata(p, organism="Escherichia")
    assert "562.A" in meta
    assert "1280.A" not in meta


def test_organism_filter_case_insensitive(bvbrc_genome_csv: Path):
    meta = load_bvbrc_genome_metadata(bvbrc_genome_csv, organism="ESCHERICHIA COLI")
    assert len(meta) == 3


# ---- edge cases ----


def test_load_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_bvbrc_genome_metadata("/tmp/nonexistent_file.csv")


def test_load_missing_genome_id_column_raises(tmp_path: Path):
    csv = "Foo,Bar\n1,2\n"
    p = tmp_path / "no_genome_id.csv"
    p.write_text(csv, encoding="utf-8")
    with pytest.raises(BvBrcGenomeError, match="Genome ID"):
        load_bvbrc_genome_metadata(p)


def test_load_blank_numeric_fields_default_to_zero(tmp_path: Path):
    csv = (
        "Genome ID,Genome Name,Species,Assembly Accession,Contigs,Contig N50,Collection Year\n"
        '"562.X","Escherichia coli X","Escherichia coli","GCF_X.1","","","N/A"\n'
    )
    p = tmp_path / "blank.csv"
    p.write_text(csv, encoding="utf-8")
    meta = load_bvbrc_genome_metadata(p)
    row = meta["562.X"]
    assert row["contig_count"] == 0
    assert row["n50"] == 0
    assert row["year"] == 0


def test_load_handles_float_str_in_int_field(tmp_path: Path):
    """BV-BRC sometimes emits '200000.0' for integer fields."""
    csv = (
        "Genome ID,Genome Name,Species,Assembly Accession,Contigs,Contig N50,Collection Year\n"
        '"562.Y","Escherichia coli Y","Escherichia coli","GCF_Y.1","10","200000.0","2018"\n'
    )
    p = tmp_path / "float_str.csv"
    p.write_text(csv, encoding="utf-8")
    meta = load_bvbrc_genome_metadata(p)
    assert meta["562.Y"]["contig_count"] == 10
    assert meta["562.Y"]["n50"] == 200000


def test_load_extra_unmapped_columns_silently_dropped(bvbrc_genome_csv: Path):
    """`Size`, `CheckM Completeness`, `Plasmids` are not in GENOME_COLUMN_MAP."""
    meta = load_bvbrc_genome_metadata(bvbrc_genome_csv)
    # No error; returned dict has only the GENOME_METADATA_KEYS
    assert set(meta["562.5691"].keys()) == set(GENOME_METADATA_KEYS)


def test_load_empty_after_filter_returns_empty_dict(tmp_path: Path):
    csv = (
        "Genome ID,Genome Name,Species,Contigs,Contig N50,Collection Year\n"
        '"1280.A","Staphylococcus aureus","Staphylococcus aureus","10","200000","2020"\n'
    )
    p = tmp_path / "no_ecoli.csv"
    p.write_text(csv, encoding="utf-8")
    meta = load_bvbrc_genome_metadata(p, organism="Escherichia coli")
    assert meta == {}


def test_load_drops_rows_without_assembly_accession(tmp_path: Path):
    """Rows lacking assembly_accession are dropped (cannot download via NCBI)."""
    csv = (
        "Genome ID,Genome Name,Species,Assembly Accession,Contigs,Contig N50,Collection Year\n"
        '"562.A","Escherichia coli A","Escherichia coli","GCF_AAA.1","10","200000","2020"\n'
        '"562.B","Escherichia coli B","Escherichia coli","","20","150000","2020"\n'
        '"562.C","Escherichia coli C","Escherichia coli","GCF_CCC.1","30","180000","2020"\n'
    )
    p = tmp_path / "mixed_acc.csv"
    p.write_text(csv, encoding="utf-8")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        meta = load_bvbrc_genome_metadata(p)
        assert any("dropped 1" in str(w.message) and "assembly_accession" in str(w.message) for w in caught)
    assert set(meta.keys()) == {"562.A", "562.C"}
    assert "562.B" not in meta


def test_load_duplicate_genome_id_emits_warning(tmp_path: Path):
    csv = (
        "Genome ID,Genome Name,Species,Assembly Accession,Contigs,Contig N50,Collection Year\n"
        '"562.Z","Escherichia coli Z","Escherichia coli","GCF_Z.1","10","200000","2018"\n'
        '"562.Z","Escherichia coli Z again","Escherichia coli","GCF_Z.2","20","100000","2020"\n'
    )
    p = tmp_path / "dupe.csv"
    p.write_text(csv, encoding="utf-8")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        meta = load_bvbrc_genome_metadata(p)
        assert any("duplicate" in str(w.message).lower() for w in caught)
    # Last-write-wins
    assert meta["562.Z"]["contig_count"] == 20
    assert meta["562.Z"]["n50"] == 100000
