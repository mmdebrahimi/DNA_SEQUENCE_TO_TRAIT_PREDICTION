"""Step 4 tests — DeterminantHit parse + the hard join-quality gate."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from dna_decode.genome_map.phenotype_overlay import (
    DeterminantHit,
    JoinedHit,
    all_joins_symbol_fallback,
    build_contig_name_map,
    determinant_phenotype_field,
    join_hits,
    parse_determinant_hits,
)


def _features() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"seqid": "contig_1", "start": 100, "end": 900, "strand": "+",
             "gene_id": "GMAP_001", "gene_symbol": "gyrA", "locus_tag": "TAG_001",
             "product": "DNA gyrase subunit A", "type": "CDS", "source": "Bakta"},
            {"seqid": "contig_1", "start": 2000, "end": 3000, "strand": "-",
             "gene_id": "GMAP_002", "gene_symbol": "", "locus_tag": "TAG_002",
             "product": "blaCTX-M-15", "type": "CDS", "source": "Bakta"},
            {"seqid": "contig_2", "start": 50, "end": 500, "strand": "+",
             "gene_id": "GMAP_003", "gene_symbol": "", "locus_tag": "TAG_003",
             "product": "hypothetical protein", "type": "CDS", "source": "Bakta"},
        ]
    )


# ---- parse_determinant_hits ----


def test_parse_with_coords(tmp_path: Path):
    tsv = tmp_path / "main.tsv"
    tsv.write_text(
        "Protein identifier\tContig id\tStart\tStop\tElement symbol\tElement name\tClass\tSubclass\tMethod\n"
        "GMAP_001\tcontig_1\t100\t900\tgyrA_S83L\tgyrase\tQUINOLONE\tQUINOLONE\tPOINTX\n",
        encoding="utf-8",
    )
    hits = parse_determinant_hits(tsv)
    assert len(hits) == 1
    h = hits[0]
    assert h.symbol == "gyrA_S83L"
    assert h.protein_id == "GMAP_001"
    assert h.contig == "contig_1"
    assert h.start == 100 and h.stop == 900
    assert h.cls == "QUINOLONE"


def test_parse_na_coords_become_none(tmp_path: Path):
    tsv = tmp_path / "main.tsv"
    tsv.write_text(
        "Protein identifier\tContig id\tStart\tStop\tGene symbol\tClass\tMethod\n"
        "NA\tNA\tNA\tNA\tblaTEM-1\tBETA-LACTAM\tEXACTX\n",
        encoding="utf-8",
    )
    h = parse_determinant_hits(tsv)[0]
    assert h.protein_id is None and h.contig is None
    assert h.start is None and h.stop is None
    assert h.symbol == "blaTEM-1"


def test_parse_missing_file(tmp_path: Path):
    assert parse_determinant_hits(tmp_path / "absent.tsv") == []


# ---- join_hits ----


def test_protein_id_join_is_high_confidence():
    hits = [DeterminantHit("gyrA_S83L", "", "QUINOLONE", "", "POINTX",
                           protein_id="GMAP_001", contig=None, start=None, stop=None)]
    joined, counts = join_hits(_features(), hits)
    assert joined[0].join_confidence == "protein_id"
    assert joined[0].is_high_confidence
    assert counts["n_high_confidence_join"] == 1
    assert counts["n_symbol_fallback"] == 0


def test_coord_join_is_high_confidence():
    # no protein-id match, but contig+coords overlap feature 2 (blaCTX-M-15)
    hits = [DeterminantHit("blaCTX-M-15", "", "CEPHALOSPORIN", "CEPHALOSPORIN", "EXACTX",
                           protein_id=None, contig="contig_1", start=2010, stop=2900)]
    joined, counts = join_hits(_features(), hits)
    assert joined[0].join_confidence == "coord"
    assert joined[0].feature_index == 1
    assert counts["n_high_confidence_join"] == 1


def test_coord_join_with_contig_name_map():
    # AMRFinder contig 'CP012345.1' reconciles to Bakta 'contig_1'
    hits = [DeterminantHit("blaCTX-M-15", "", "CEPHALOSPORIN", "", "EXACTX",
                           protein_id=None, contig="CP012345.1", start=2010, stop=2900)]
    cmap = {"CP012345.1": "contig_1"}
    joined, counts = join_hits(_features(), hits, contig_name_map=cmap)
    assert joined[0].join_confidence == "coord"
    assert joined[0].feature_index == 1


def test_symbol_fallback_when_no_coords():
    # coords absent -> falls to symbol_fallback (gyrA), NOT high-confidence
    hits = [DeterminantHit("gyrA_S83L", "", "QUINOLONE", "", "POINTX",
                           protein_id=None, contig=None, start=None, stop=None)]
    joined, counts = join_hits(_features(), hits)
    assert joined[0].join_confidence == "symbol_fallback"
    assert not joined[0].is_high_confidence
    assert counts["n_symbol_fallback"] == 1
    assert counts["n_high_confidence_join"] == 0


def test_unjoined_hit_counted():
    hits = [DeterminantHit("mystery", "", "OTHER", "", "EXACTX",
                           protein_id=None, contig=None, start=None, stop=None)]
    joined, counts = join_hits(_features(), hits)
    assert joined[0].feature_index is None
    assert counts["n_unjoined"] == 1


# ---- all_joins_symbol_fallback (the NO-GO guard) ----


def test_all_symbol_fallback_true():
    counts = {"n_main_rows": 3, "n_high_confidence_join": 0, "n_symbol_fallback": 2, "n_unjoined": 1}
    assert all_joins_symbol_fallback(counts) is True


def test_all_symbol_fallback_false_when_high_conf_present():
    counts = {"n_main_rows": 3, "n_high_confidence_join": 1, "n_symbol_fallback": 2, "n_unjoined": 0}
    assert all_joins_symbol_fallback(counts) is False


def test_all_symbol_fallback_false_when_no_determinants():
    counts = {"n_main_rows": 0, "n_high_confidence_join": 0, "n_symbol_fallback": 0, "n_unjoined": 0}
    assert all_joins_symbol_fallback(counts) is False


# ---- build_contig_name_map ----


def test_build_contig_name_map_by_unique_length():
    fasta = {"CP012345.1": 5000000, "CP012346.1": 90000}
    bakta = {"contig_1": 5000000, "contig_2": 90000}
    cmap = build_contig_name_map(fasta, bakta)
    assert cmap == {"CP012345.1": "contig_1", "CP012346.1": "contig_2"}


def test_build_contig_name_map_skips_ambiguous_lengths():
    fasta = {"a": 100, "b": 100, "c": 300}
    bakta = {"x": 100, "y": 100, "z": 300}
    cmap = build_contig_name_map(fasta, bakta)
    # length 100 is ambiguous on both sides -> only the unique 300 maps
    assert cmap == {"c": "z"}


# ---- determinant_phenotype_field (the phenotype wall + ABSTAIN) ----


def test_phenotype_field_high_confidence():
    jh = JoinedHit(DeterminantHit("gyrA_S83L", "", "QUINOLONE", "", "POINTX",
                                  "GMAP_001", "contig_1", 100, 900),
                   feature_index=0, join_confidence="coord")
    field = determinant_phenotype_field(jh, "ciprofloxacin", {"prediction": "R", "rule": "rule_x"})
    assert field is not None
    assert field["phenotype"] == "R"
    assert field["drug"] == "ciprofloxacin"
    assert field["abstain"] is False


def test_phenotype_field_symbol_fallback_is_none():
    jh = JoinedHit(DeterminantHit("gyrA", "", "QUINOLONE", "", "POINTX",
                                  None, None, None, None),
                   feature_index=0, join_confidence="symbol_fallback")
    # the phenotype wall: symbol-fallback never carries a phenotype
    assert determinant_phenotype_field(jh, "ciprofloxacin", {"prediction": "R"}) is None


def test_phenotype_field_abstain_propagates():
    jh = JoinedHit(DeterminantHit("blaX", "", "CEPHALOSPORIN", "", "EXACTX",
                                  "GMAP_002", "contig_1", 2000, 3000),
                   feature_index=1, join_confidence="protein_id")
    field = determinant_phenotype_field(jh, "ceftriaxone", {"prediction": "ABSTAIN", "rule": "cal"})
    assert field["phenotype"] == "ABSTAIN"
    assert field["abstain"] is True
