"""Step 2 tests — tool-surface manifest builders (pure, no Docker)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.genome_map_tool_surface import (
    build_amrfinder_inventory,
    build_bakta_inventory,
    build_manifest,
)


def _synth_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"seqid": "c1", "source": "Bakta", "type": "CDS", "start": 1, "end": 9,
             "strand": "+", "gene_id": "g1", "gene_symbol": "gyrA", "locus_tag": "T1",
             "product": "DNA gyrase subunit A"},
            {"seqid": "c1", "source": "Bakta", "type": "CDS", "start": 10, "end": 20,
             "strand": "-", "gene_id": "g2", "gene_symbol": "", "locus_tag": "T2",
             "product": "putative oxidoreductase"},
            {"seqid": "c1", "source": "Bakta", "type": "CDS", "start": 21, "end": 30,
             "strand": "+", "gene_id": "g3", "gene_symbol": "", "locus_tag": "T3",
             "product": "hypothetical protein"},
            {"seqid": "c1", "source": "Bakta", "type": "gene", "start": 1, "end": 9,
             "strand": "+", "gene_id": "g1", "gene_symbol": "gyrA", "locus_tag": "",
             "product": ""},
        ]
    )


# ---- bakta inventory ----


def test_bakta_inventory_counts_and_vocab():
    inv = build_bakta_inventory(_synth_table())
    assert inv["n_features"] == 4
    assert inv["feature_type_counts"]["CDS"] == 3
    assert inv["feature_type_counts"]["gene"] == 1
    assert "DNA gyrase subunit A" in inv["product_vocabulary_sample"]
    # putative -> low-confidence example
    assert any("putative" in p for p in inv["low_confidence_product_examples"])
    # one hypothetical + one empty product
    assert inv["hypothetical_or_empty_count"] == 2
    assert "product" in inv["populated_fields"]


# ---- amrfinder inventory ----


def test_amrfinder_inventory_with_coords(tmp_path: Path):
    tsv = tmp_path / "main.tsv"
    tsv.write_text(
        "Protein identifier\tContig id\tStart\tStop\tElement symbol\tClass\tMethod\n"
        "WP_001\tcontig_1\t100\t900\tgyrA_S83L\tQUINOLONE\tPOINTX\n",
        encoding="utf-8",
    )
    inv = build_amrfinder_inventory(tsv)
    assert inv["present"] is True
    assert inv["has_protein_id"] is True
    assert inv["has_contig"] is True
    assert inv["has_coords"] is True
    assert inv["symbol_column"] == "Element symbol"
    assert inv["n_rows"] == 1


def test_amrfinder_inventory_no_coords(tmp_path: Path):
    tsv = tmp_path / "main.tsv"
    tsv.write_text(
        "Gene symbol\tClass\tMethod\n"
        "blaCTX-M-15\tCEPHALOSPORIN\tEXACTX\n",
        encoding="utf-8",
    )
    inv = build_amrfinder_inventory(tsv)
    assert inv["has_coords"] is False
    assert inv["has_protein_id"] is False
    assert inv["symbol_column"] == "Gene symbol"


def test_amrfinder_inventory_missing_file(tmp_path: Path):
    inv = build_amrfinder_inventory(tmp_path / "absent.tsv")
    assert inv["present"] is False
    assert inv["n_rows"] == 0
    assert inv["has_coords"] is False


# ---- manifest assembly + BLOCKED ----


def test_build_manifest_ok(tmp_path: Path):
    tsv = tmp_path / "main.tsv"
    tsv.write_text("Element symbol\tClass\tMethod\nblaTEM-1\tBETA-LACTAM\tEXACTX\n", encoding="utf-8")
    m = build_manifest("GCA_1", "Escherichia", _synth_table(), tsv, generated="2026-06-18")
    assert m["status"] == "OK"
    assert m["bakta"]["n_features"] == 4
    assert m["amrfinder"]["present"] is True
    assert m["amrfinder_organism"] == "Escherichia"
    assert m["generated"] == "2026-06-18"


def test_build_manifest_bakta_blocked():
    m = build_manifest("GCA_1", "Escherichia", None, None, bakta_status="BAKTA_ANNOTATION_BLOCKED")
    assert m["status"] == "BAKTA_ANNOTATION_BLOCKED"
    assert m["bakta"] is None


def test_build_manifest_amrfinder_blocked():
    m = build_manifest("GCA_1", None, _synth_table(), None, amrfinder_status="AMRFINDER_BLOCKED")
    assert m["status"] == "AMRFINDER_BLOCKED"
    assert m["amrfinder"] is None
    assert m["bakta"] is not None
