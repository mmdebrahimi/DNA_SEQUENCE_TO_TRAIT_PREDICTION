"""Guards for the label->annotation join — the ladder's integrity crux.

The FBA arm's S. aureus / P. aeruginosa "walls" are model-gene-join walls that do not apply to this
decoder, but the wall MOVES here: label ids must join the exact annotation table. Two recorded failures
set the bar — a cross-strain id join that gave 0% overlap and a confident AUROC 0.000, and the fact that
a join RATE cannot see a biased joined subset.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import pytest

from dna_decode.essentiality import annotation_join as aj

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path("D:/dna_decode_cache/essentiality")
_HAVE_ECOLI = (CACHE / "ecoli_k12_feature_table.txt.gz").exists()


def _ann(n_keys, prefix="g"):
    return {f"{prefix}{i}": f"product {i}" for i in range(n_keys)}


# --------------------------------------------------------------------------------------------------
# join rate + the gate
# --------------------------------------------------------------------------------------------------
def test_full_join_rate():
    ann = _ann(4)
    r = aj.join_labels(ann, ["g0", "g1"], ["g2", "g3"], key_space="gene_symbol")
    assert r.join_rate == 1.0 and r.n_unjoined == 0 and r.wall is None
    assert len(r.rows) == 4


def test_the_floor_is_INCLUSIVE_as_documented():
    """0.70 exactly must PASS; the bar reads 'below it the rung walls'."""
    ann = _ann(7)
    ess = [f"g{i}" for i in range(7)]
    non = [f"absent{i}" for i in range(3)]          # 7/10 = 0.70 exactly
    r = aj.join_labels(ann, ess, non, key_space="gene_symbol")
    assert r.join_rate == 0.7
    assert r.wall is None, "0.70 is the floor and must be inclusive"


def test_just_below_the_floor_walls():
    ann = _ann(6)
    ess = [f"g{i}" for i in range(6)]
    non = [f"absent{i}" for i in range(4)]          # 6/10 = 0.60
    r = aj.join_labels(ann, ess, non, key_space="gene_symbol")
    assert r.join_rate == 0.6
    assert r.wall == aj.WALL_JOIN_RATE_BELOW_FLOOR


def test_a_walled_rung_carries_NO_coverage_key_and_NO_scorable_rows():
    """Absence, not a None beside a number that looks computed. A walled rung must also be unable to
    hand rows downstream, or a caller could score it anyway."""
    ann = _ann(1)
    r = aj.join_labels(ann, ["g0"], [f"absent{i}" for i in range(19)], key_space="gene_symbol")
    assert r.wall == aj.WALL_JOIN_RATE_BELOW_FLOOR
    assert r.rows == []
    d = r.as_dict()
    assert "coverage" not in d and "coverage_lift" not in d
    assert d["wall"] == aj.WALL_JOIN_RATE_BELOW_FLOOR


def test_a_catastrophic_join_reproduces_the_recorded_trap_shape():
    """The gene_id-vs-gene_symbol disaster: ~0% overlap. It must wall, not report a number."""
    ann = {"gene-b0001": "x", "gene-b0002": "y"}         # strain-unique ids
    r = aj.join_labels(ann, ["gyrA", "parC"], ["lacZ"], key_space="gene_symbol")
    assert r.join_rate == 0.0
    assert r.wall == aj.WALL_JOIN_RATE_BELOW_FLOOR


# --------------------------------------------------------------------------------------------------
# composition — the bias a rate cannot see
# --------------------------------------------------------------------------------------------------
def test_composition_detects_a_conserved_enriched_joined_subset_at_an_ACCEPTABLE_rate():
    """THE review's point. Both halves of this test have the same join rate; only one is biased."""
    # 8 of 10 join either way -> rate 0.8 in both cases
    unbiased = aj.join_labels(_ann(8), ["g0", "g1", "g2", "g3", "x1"],
                              ["g4", "g5", "g6", "g7", "x2"], key_space="gene_symbol")
    biased = aj.join_labels(_ann(8), ["g0", "g1", "g2", "g3", "g4", "g5", "g6", "g7"],
                            ["x1", "x2"], key_space="gene_symbol")
    assert unbiased.join_rate == biased.join_rate == 0.8
    assert unbiased.composition["essential_enrichment_in_joined"] == 0.0
    # every joined label is essential, every unjoined one is not -> maximal enrichment
    assert biased.composition["essential_enrichment_in_joined"] == 1.0
    assert biased.composition["joined_essential_fraction"] == 1.0
    assert biased.composition["unjoined_essential_fraction"] == 0.0


def test_composition_is_reported_even_when_nothing_is_unjoined():
    r = aj.join_labels(_ann(2), ["g0"], ["g1"], key_space="gene_symbol")
    c = r.composition
    assert c["unjoined_n"] == 0
    assert c["unjoined_essential_fraction"] is None
    assert c["essential_enrichment_in_joined"] is None, "no comparison is possible, so no number"


# --------------------------------------------------------------------------------------------------
# normalisation — explicit, never fuzzy
# --------------------------------------------------------------------------------------------------
def test_locus_tag_zero_padding():
    assert aj.normalise_key("SAUSA300_1", "sausa300_locus_tag") == "SAUSA300_0001"
    assert aj.normalise_key("sausa300_0001", "sausa300_locus_tag") == "SAUSA300_0001"
    assert aj.normalise_key("PA14_1", "pa14_locus_tag") == "PA14_0001"


def test_orf_names_are_case_folded_but_gene_symbols_are_NOT():
    assert aj.normalise_key("yal001c", "systematic_orf") == "YAL001C"
    # gene symbol case is meaningful -- recA and RecA are not interchangeable
    assert aj.normalise_key("gyrA", "gene_symbol") == "gyrA"
    assert aj.normalise_key("GYRA", "gene_symbol") == "GYRA"


def test_an_unknown_key_space_RAISES_rather_than_guessing():
    """A fuzzy fallback is exactly how a strain-unique id silently joins the wrong gene."""
    with pytest.raises(aj.AnnotationError, match="no normalisation rule"):
        aj.normalise_key("x", "some_new_key_space")


# --------------------------------------------------------------------------------------------------
# annotation loading
# --------------------------------------------------------------------------------------------------
def test_unknown_annotation_kind_raises():
    with pytest.raises(aj.AnnotationError, match="unknown annotation kind"):
        aj.load_annotation("whatever", "not_a_kind")


def test_an_empty_annotation_RAISES_rather_than_joining_nothing(tmp_path):
    """An empty annotation joins nothing and would yield a clean-looking 0.0 coverage on no data."""
    p = tmp_path / "ft.txt"
    p.write_text("# feature\tsymbol\tname\n", encoding="utf-8")
    with pytest.raises(aj.AnnotationError, match="ZERO entries"):
        aj.load_annotation(str(p), "ncbi_feature_table")


def test_feature_table_reader_keeps_only_CDS_rows_with_a_symbol(tmp_path):
    p = tmp_path / "ft.txt"
    p.write_text("# feature\tsymbol\tname\n"
                 "CDS\tgyrA\tDNA gyrase subunit A\n"
                 "gene\tgyrA\tshould be skipped (not CDS)\n"
                 "CDS\t\tno symbol so skipped\n"
                 "CDS\tparC\ttopoisomerase IV subunit A\n", encoding="utf-8")
    ann = aj.load_annotation(str(p), "ncbi_feature_table")
    assert ann == {"gyrA": "DNA gyrase subunit A", "parC": "topoisomerase IV subunit A"}


def test_sgd_reader_keeps_only_ORF_rows(tmp_path):
    p = tmp_path / "sgd.tab"
    cols = [""] * 16
    orf = list(cols); orf[3] = "ORF"; orf[5] = "YAL001C"; orf[4] = "TFC3"; orf[15] = "subunit of TFIIIC"
    notorf = list(cols); notorf[3] = "ncRNA"; notorf[5] = "YNC001"; notorf[15] = "noise"
    p.write_text("\t".join(orf) + "\n" + "\t".join(notorf) + "\n", encoding="utf-8")
    ann = aj.load_annotation(str(p), "sgd_features")
    assert list(ann) == ["YAL001C"]
    assert "TFIIIC" in ann["YAL001C"]


def test_gene_info_reader_requires_its_columns(tmp_path):
    p = tmp_path / "gi.tsv"
    p.write_text("GeneID\tSymbol\n1\tA\n", encoding="utf-8")
    with pytest.raises(aj.AnnotationError, match="required column"):
        aj.load_annotation(str(p), "ncbi_gene_info")


# --------------------------------------------------------------------------------------------------
# the promotion is real
# --------------------------------------------------------------------------------------------------
def test_cli_imports_the_shared_reader_rather_than_defining_its_own():
    """Guards re-divergence: two readers of one file format is how they drift."""
    src = (ROOT / "dna_decode" / "essentiality" / "cli.py").read_text(encoding="utf-8")
    assert "from dna_decode.essentiality.annotation_join import iter_feature_table" in src
    assert "header.index(\"symbol\")" not in src, "the CLI must no longer parse the table itself"


@pytest.mark.skipif(not _HAVE_ECOLI, reason="D: essentiality cache absent")
def test_promoted_reader_matches_the_old_cli_helper_on_the_REAL_table():
    """Byte-identical output on real data, not just on a fixture."""
    from dna_decode.essentiality.cli import _iter_feature_table

    path = str(CACHE / "ecoli_k12_feature_table.txt.gz")
    via_cli = list(_iter_feature_table(path))
    via_shared = list(aj.iter_feature_table(path))
    assert via_cli == via_shared
    assert len(via_shared) > 4000, len(via_shared)
