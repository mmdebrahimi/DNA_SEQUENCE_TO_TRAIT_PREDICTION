"""Tests for the MaveDB prospective-holdout pure functions (no network)."""
from scripts.mavedb_prospective_holdout import (
    build_manifest,
    is_held_out,
    parse_hgvs_pro,
)


def test_parse_hgvs_pro_missense():
    assert parse_hgvs_pro("p.Val82Ala") == ("V", 82, "A")
    assert parse_hgvs_pro("p.Trp57Gly") == ("W", 57, "G")


def test_parse_hgvs_pro_rejects_non_single_missense():
    assert parse_hgvs_pro("p.=") is None            # synonymous
    assert parse_hgvs_pro("p.Trp82Ter") is None     # nonsense
    assert parse_hgvs_pro("p.[Val82Ala;Trp57Gly]") is None  # multi
    assert parse_hgvs_pro("p.Val82del") is None     # deletion
    assert parse_hgvs_pro("c.45T>C") is None        # nucleotide, not protein
    assert parse_hgvs_pro("") is None
    assert parse_hgvs_pro("p.Xyz82Ala") is None     # non-standard aa


def test_is_held_out():
    pg = {"BRCA1", "TP53", "PTEN"}
    assert is_held_out("SomeNovelGene", pg) is True      # not in benchmark -> held out
    assert is_held_out("BRCA1", pg) is False             # in benchmark -> leaked, excluded
    assert is_held_out("brca1", pg) is False             # case/normalization
    assert is_held_out("", pg) is False                  # empty gene -> not held out


def test_build_manifest_filters_and_flags():
    pg = {"BRCA1"}
    score_sets = [
        {"urn": "u1", "publishedDate": "2025-03-01", "numVariants": 500,
         "license": {"shortName": "CC0"},
         "targetGenes": [{"name": "NovelHuman", "category": "protein_coding",
                          "targetSequence": {"taxonomy": {"organismName": "Homo sapiens"}}}]},
        {"urn": "u2", "publishedDate": "2022-01-01", "numVariants": 100,
         "license": {"shortName": "CC0"},
         "targetGenes": [{"name": "BRCA1", "category": "protein_coding",  # in ProteinGym -> excluded
                          "targetSequence": {"taxonomy": {"organismName": "Homo sapiens"}}}]},
        {"urn": "u3", "publishedDate": "2025-01-01", "numVariants": 50,
         "license": {"shortName": "CC0"},
         "targetGenes": [{"name": "someRNA", "category": "regulatory",   # not protein_coding -> excluded
                          "targetSequence": {"taxonomy": {"organismName": "Homo sapiens"}}}]},
    ]
    m = build_manifest(score_sets, pg, cutoff="2024-01-01")
    assert [r["urn"] for r in m] == ["u1"]        # only the held-out protein_coding one
    assert m[0]["post_cutoff"] is True            # 2025 >= 2024 cutoff
    assert m[0]["organism"] == "Homo sapiens"
