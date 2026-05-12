"""Tests for Step 4 — Resistance database loaders."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.data.resistance_db import (
    ResistanceCatalog,
    ResistanceEntry,
    load_amrfinder,
    load_card,
    merge_catalogs,
)


MOCK_CARD_JSON = {
    "ARO_3000027": {
        "model_name": "gyrA",
        "ARO_accession": "3000027",
        "ARO_category": {
            "fam1": {
                "category_aro_class_name": "AMR Gene Family",
                "category_aro_name": "fluoroquinolone resistant gyrA",
            },
            "drug1": {
                "category_aro_class_name": "Drug Class",
                "category_aro_name": "fluoroquinolone",
            },
            "mech1": {
                "category_aro_class_name": "Resistance Mechanism",
                "category_aro_name": "antibiotic target alteration",
            },
        },
    },
    "ARO_3001877": {
        "model_name": "blaCTX-M-15",
        "ARO_accession": "3001877",
        "ARO_category": {
            "fam1": {
                "category_aro_class_name": "AMR Gene Family",
                "category_aro_name": "CTX-M beta-lactamase",
            },
            "drug1": {
                "category_aro_class_name": "Drug Class",
                "category_aro_name": "cephalosporin",
            },
        },
    },
    "non_record_string_value": "this should be skipped (not a dict)",
}

MOCK_AMRFINDER_TSV = """\
gene_symbol\tgene_family\tclass\tresistance_mechanism\taccession
gyrA\tfluoroquinolone-resistance\tFLUOROQUINOLONE\ttarget alteration\tWP_001234
tetA\ttetracycline-efflux\tTETRACYCLINE\tefflux pump\tWP_005678
\textra_field\tNONE\tNONE\tNONE
blaCTX-M-15\tbeta-lactamase\tCEPHALOSPORIN\thydrolysis\tWP_009999
"""


# ---- dataclass + catalog mechanics ----


def test_catalog_starts_empty():
    cat = ResistanceCatalog()
    assert len(cat) == 0
    assert cat.all_gene_symbols() == set()


def test_catalog_add_and_lookup():
    cat = ResistanceCatalog()
    e = ResistanceEntry(
        gene_symbol="gyrA",
        gene_family="fluoroquinolone-resistance",
        drug_class="fluoroquinolone",
        resistance_mechanism="target alteration",
        source_db="CARD",
        source_id="3000027",
    )
    cat.add(e)
    assert len(cat) == 1
    assert cat.map_gene_to_resistance("gyrA") == [e]


def test_catalog_lookup_is_case_insensitive():
    cat = ResistanceCatalog()
    cat.add(
        ResistanceEntry("gyrA", "f", "fluoroquinolone", "alt", "CARD", "X")
    )
    assert cat.map_gene_to_resistance("GYRA")[0].gene_symbol == "gyrA"
    assert cat.map_gene_to_resistance("gyra")[0].gene_symbol == "gyrA"


def test_catalog_filter_by_drug_class():
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "fluoroquinolone", "", "CARD", "1"))
    cat.add(ResistanceEntry("tetA", "t", "tetracycline", "", "CARD", "2"))
    matches = cat.filter_by_drug_class("fluoroquinolone")
    assert len(matches) == 1
    assert matches[0].gene_symbol == "gyrA"


def test_catalog_filter_by_drug_class_case_insensitive():
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "Fluoroquinolone", "", "CARD", "1"))
    matches = cat.filter_by_drug_class("FLUOROQUINOLONE")
    assert len(matches) == 1


# ---- CARD loader ----


def test_load_card_basic(tmp_path: Path):
    p = tmp_path / "card.json"
    p.write_text(json.dumps(MOCK_CARD_JSON), encoding="utf-8")
    cat = load_card(p)
    assert len(cat) == 2  # 2 dict-record entries; non-dict string value skipped
    symbols = cat.all_gene_symbols()
    assert "gyrA" in symbols
    assert "blaCTX-M-15" in symbols


def test_load_card_extracts_drug_class(tmp_path: Path):
    p = tmp_path / "card.json"
    p.write_text(json.dumps(MOCK_CARD_JSON), encoding="utf-8")
    cat = load_card(p)
    gyra = cat.map_gene_to_resistance("gyrA")[0]
    assert gyra.drug_class == "fluoroquinolone"
    assert gyra.source_db == "CARD"
    assert gyra.gene_family == "fluoroquinolone resistant gyrA"


def test_load_card_missing_drug_class_defaults_unknown(tmp_path: Path):
    p = tmp_path / "card.json"
    data = {
        "ARO_1": {
            "model_name": "mystery",
            "ARO_accession": "1",
            "ARO_category": {},  # no drug class
        }
    }
    p.write_text(json.dumps(data), encoding="utf-8")
    cat = load_card(p)
    assert cat.entries[0].drug_class == "unknown"


def test_load_card_missing_file_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_card(tmp_path / "missing.json")


# ---- AMRFinder loader ----


def test_load_amrfinder_basic(tmp_path: Path):
    p = tmp_path / "amrfinder.tsv"
    p.write_text(MOCK_AMRFINDER_TSV, encoding="utf-8")
    cat = load_amrfinder(p)
    # Empty-gene-symbol row skipped; rest accepted = 3
    assert len(cat) == 3
    assert "gyrA" in cat.all_gene_symbols()
    assert "tetA" in cat.all_gene_symbols()
    assert "blaCTX-M-15" in cat.all_gene_symbols()


def test_load_amrfinder_preserves_drug_class(tmp_path: Path):
    p = tmp_path / "amrfinder.tsv"
    p.write_text(MOCK_AMRFINDER_TSV, encoding="utf-8")
    cat = load_amrfinder(p)
    tet = cat.map_gene_to_resistance("tetA")[0]
    assert tet.drug_class == "TETRACYCLINE"
    assert tet.source_db == "AMRFinder"
    assert tet.source_id == "WP_005678"


def test_load_amrfinder_missing_file_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_amrfinder(tmp_path / "missing.tsv")


# ---- merge ----


def test_merge_catalogs_concatenates(tmp_path: Path):
    card_path = tmp_path / "card.json"
    card_path.write_text(json.dumps(MOCK_CARD_JSON), encoding="utf-8")
    amrfinder_path = tmp_path / "amrfinder.tsv"
    amrfinder_path.write_text(MOCK_AMRFINDER_TSV, encoding="utf-8")

    merged = merge_catalogs(load_card(card_path), load_amrfinder(amrfinder_path))
    # CARD: 2 + AMRFinder: 3 = 5
    assert len(merged) == 5
    # Cross-source lookup works
    gyra_hits = merged.map_gene_to_resistance("gyrA")
    sources = {e.source_db for e in gyra_hits}
    assert sources == {"CARD", "AMRFinder"}
