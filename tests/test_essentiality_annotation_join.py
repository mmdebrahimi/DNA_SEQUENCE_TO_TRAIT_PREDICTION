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


# --------------------------------------------------------------------------------------------------
# the reader paths the real rungs actually use — exercised OFFLINE, since every real-data test above is
# skipped on a host without D: and these branches would then be wholly uncovered
# --------------------------------------------------------------------------------------------------
def test_a_GZIPPED_feature_table_reads_the_same_as_a_plain_one(tmp_path):
    """The committed E. coli annotation is `ecoli_k12_feature_table.txt.gz`, so the gzip branch is the
    one the published 0.3125 was produced through -- and without this it is only covered by a D:-gated
    test."""
    body = ("# feature\tsymbol\tname\n"
            "CDS\tgyrA\tDNA gyrase subunit A\n"
            "CDS\tparC\ttopoisomerase IV subunit A\n")
    plain = tmp_path / "ft.txt"
    plain.write_text(body, encoding="utf-8")
    gz = tmp_path / "ft.txt.gz"
    gz.write_bytes(gzip.compress(body.encode("utf-8")))
    assert aj.load_annotation(str(gz), "ncbi_feature_table") == \
        aj.load_annotation(str(plain), "ncbi_feature_table")
    assert aj.load_annotation(str(gz), "ncbi_feature_table")["gyrA"] == "DNA gyrase subunit A"


def test_gene_info_keeps_protein_coding_only_and_keys_by_GeneID(tmp_path):
    """The human rung's annotation shape. Only the missing-column REFUSAL was pinned; this pins what the
    reader actually returns -- keyed by GeneID (not Symbol, which is what the ENTREZ_ID join depends on)
    and with the description appended to the symbol so the decoder reads both."""
    p = tmp_path / "gi.tsv"
    p.write_text("#tax_id\tGeneID\tSymbol\tdescription\ttype_of_gene\n"
                 "9606\t6122\tRPL3\tribosomal protein L3\tprotein-coding\n"
                 "9606\t999\tMIR1\tmicroRNA 1\tncRNA\n", encoding="utf-8")
    ann = aj.load_annotation(str(p), "ncbi_gene_info")
    assert list(ann) == ["6122"], "keyed by GeneID; the ncRNA row must be dropped"
    assert ann["6122"] == "RPL3 ribosomal protein L3"


def test_a_feature_table_whose_header_LACKS_a_column_raises_rather_than_yielding_partial_rows(tmp_path):
    """`iter_feature_table` returns silently when `symbol`/`name` are absent -- a reader that parsed
    NOTHING is indistinguishable from one that found nothing, which is the recorded trap. The empty-dict
    refusal is what converts it into an error."""
    p = tmp_path / "ft.txt"
    p.write_text("# feature\tsymbol\tproduct\nCDS\tgyrA\tDNA gyrase subunit A\n", encoding="utf-8")
    assert list(aj.iter_feature_table(str(p))) == [], "no `name` column -> the reader yields nothing"
    with pytest.raises(aj.AnnotationError, match="ZERO entries"):
        aj.load_annotation(str(p), "ncbi_feature_table")


# --------------------------------------------------------------------------------------------------
# normalisation — the remaining key spaces and the degenerate inputs
# --------------------------------------------------------------------------------------------------
def test_entrez_ids_lose_leading_zeros_but_an_all_zero_id_is_never_emptied():
    assert aj.normalise_key("00123", "entrez_gene_id") == "123"
    assert aj.normalise_key("123", "entrez_gene_id") == "123"
    # `lstrip("0")` alone would turn these into "", which joins to nothing and reads as a missing label
    assert aj.normalise_key("0", "entrez_gene_id") == "0"
    assert aj.normalise_key("000", "entrez_gene_id") == "000"


def test_locus_tag_padding_only_touches_a_NUMERIC_tail():
    """Zero-padding a non-numeric tail would invent an id. Dashes are folded to underscores because both
    spellings occur in the wild; everything else is uppercased and left alone."""
    assert aj.normalise_key("PA14_abc", "pa14_locus_tag") == "PA14_ABC"
    assert aj.normalise_key("pa14abc", "pa14_locus_tag") == "PA14ABC", "no separator -> no padding"
    assert aj.normalise_key("sausa300-12", "sausa300_locus_tag") == "SAUSA300_0012"


def test_a_blank_identifier_normalises_to_empty_in_EVERY_key_space_and_joins_to_nothing():
    """A blank id must not become a join key. The blank check runs before the key-space dispatch, so it
    is the one input that is handled identically everywhere -- including for a space with no rule."""
    for space in ("gene_symbol", "systematic_orf", "entrez_gene_id", "pa14_locus_tag",
                  "sausa300_locus_tag", "a_key_space_with_no_rule"):
        assert aj.normalise_key("   ", space) == "", space
    r = aj.join_labels(_ann(3), ["g0", "g1"], ["", "  "], key_space="gene_symbol",
                       min_join_rate=0.0)
    assert r.n_joined == 2 and r.n_unjoined == 2, "the two blank ids must count as UNJOINED"


# --------------------------------------------------------------------------------------------------
# join_labels — the remaining refusals, the override, and an asymmetry worth knowing about
# --------------------------------------------------------------------------------------------------
def test_no_labels_at_all_refuses():
    """Sibling of the empty-annotation and zero-id refusals: 0/0 is not a join rate."""
    with pytest.raises(aj.AnnotationError, match="empty label set"):
        aj.join_labels(_ann(2), [], [], key_space="gene_symbol")


def test_an_explicit_min_join_rate_is_HONOURED_in_both_directions():
    """The floor is a parameter, so a caller can tighten or loosen it -- and a parameter that was
    silently ignored would leave the documented default as the only real bar."""
    ann = _ann(8)
    ess = [f"g{i}" for i in range(8)]
    non = ["absent1", "absent2"]                      # 8/10 = 0.80, which clears the 0.70 default
    assert aj.join_labels(ann, ess, non, key_space="gene_symbol").wall is None
    strict = aj.join_labels(ann, ess, non, key_space="gene_symbol", min_join_rate=0.95)
    assert strict.wall == aj.WALL_JOIN_RATE_BELOW_FLOOR and strict.rows == []
    loose = aj.join_labels(_ann(1), ["g0"], ["x1", "x2", "x3"], key_space="gene_symbol",
                           min_join_rate=0.0)
    assert loose.join_rate == 0.25 and loose.wall is None and len(loose.rows) == 1


def test_an_annotation_side_normalisation_COLLISION_is_FIRST_wins():
    """Worth pinning because it is the OPPOSITE of `load_annotation`'s documented last-wins duplicate
    rule: there the later ROW of a duplicated symbol wins, here the first annotation entry that
    normalises to a given key wins (`setdefault`). Two ids that differ only in zero-padding collapse to
    one entry, and which product text survives decides the score."""
    ann = {"SAUSA300_1": "first text", "SAUSA300_0001": "second text"}
    r = aj.join_labels(ann, ["SAUSA300_0001"], ["SAUSA300_9999"],
                       key_space="sausa300_locus_tag", min_join_rate=0.0)
    assert r.rows == [("SAUSA300_0001", "first text", True)]


def test_composition_is_reported_even_when_NOTHING_joined():
    """Mirror of the nothing-unjoined case. A rung that joined zero labels still reports the makeup of
    what it failed to join -- that is the only diagnostic available on a total-miss rung."""
    r = aj.join_labels(_ann(1), ["x1"], ["x2"], key_space="gene_symbol")
    assert r.wall == aj.WALL_JOIN_RATE_BELOW_FLOOR
    c = r.composition
    assert c["joined_n"] == 0 and c["unjoined_n"] == 2
    assert c["joined_essential_fraction"] is None
    assert c["unjoined_essential_fraction"] == 0.5
    assert c["essential_enrichment_in_joined"] is None, "no comparison is possible, so no number"


@pytest.mark.skipif(not _HAVE_ECOLI, reason="D: essentiality cache absent")
def test_promoted_reader_matches_the_old_cli_helper_on_the_REAL_table():
    """Byte-identical output on real data, not just on a fixture."""
    from dna_decode.essentiality.cli import _iter_feature_table

    path = str(CACHE / "ecoli_k12_feature_table.txt.gz")
    via_cli = list(_iter_feature_table(path))
    via_shared = list(aj.iter_feature_table(path))
    assert via_cli == via_shared
    assert len(via_shared) > 4000, len(via_shared)
