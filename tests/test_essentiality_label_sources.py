"""Guards for the label registry — the gate that makes "no parser runs on an unverified schema" real.

Verification Signal 6 of the plan depends entirely on these: if `parse_or_refuse` could be reached
without a W0 record, the probe built in Step 2 would be decoration.
"""
from __future__ import annotations

import json

import pytest

from dna_decode.essentiality import label_sources as ls
from dna_decode.eval import transfer_ladder as tl

BAGEL_TEXT = "GENE\tHGNC_ID\tENTREZ_ID\nAARS\tHGNC:20\t16\nRPL3\tHGNC:10332\t6122\n"


def _rec(source_id="bagel_ceg", verified=True, cands=(0, 2)):
    return {"record": "essentiality-ladder-w0-probe-v1",
            "sources": [{"source_id": source_id, "w0_verified": verified, "reason": "r",
                         "schema": {"key_column_candidates": list(cands)}}]}


# --------------------------------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------------------------------
def test_no_w0_record_at_all_refuses():
    """Tested on `w0_status` directly, because `record=None` MEANS 'look up the committed record' --
    passing None would fall through to the real file and pass for the wrong reason."""
    ok, why = ls.w0_status("bagel_ceg", None)
    assert not ok
    assert "no W0 probe record" in why


def test_an_empty_record_refuses():
    with pytest.raises(ls.SchemaUnverified, match="absent from the W0 record"):
        ls.parse_or_refuse("bagel_ceg", BAGEL_TEXT, record={"sources": []})


def test_a_record_that_says_unverified_refuses():
    with pytest.raises(ls.SchemaUnverified, match="unverified"):
        ls.parse_or_refuse("bagel_ceg", BAGEL_TEXT, record=_rec(verified=False))


def test_a_source_absent_from_the_record_refuses():
    with pytest.raises(ls.SchemaUnverified, match="absent from the W0 record"):
        ls.parse_or_refuse("bagel_neg", BAGEL_TEXT, record=_rec(source_id="bagel_ceg"))


def test_an_unknown_source_id_refuses():
    with pytest.raises(ls.SchemaUnverified, match="unknown label source"):
        ls.parse_or_refuse("not_a_source", BAGEL_TEXT, record=_rec())


def test_a_verified_record_lets_the_parser_run():
    ids = ls.parse_or_refuse("bagel_ceg", BAGEL_TEXT, record=_rec())
    assert ids == {"16", "6122"}


def test_a_key_column_the_schema_does_not_list_refuses():
    """The parser and the pinned schema must agree about where the key lives."""
    with pytest.raises(ls.SchemaUnverified, match="not among the W0 record's candidates"):
        ls.parse_or_refuse("bagel_ceg", BAGEL_TEXT, record=_rec(cands=(0, 1)))


def test_a_zero_id_parse_RAISES_rather_than_returning_an_empty_set():
    """An empty label set yields a clean-looking 0.0 coverage computed on no data."""
    with pytest.raises(ls.LabelParseError, match="ZERO ids"):
        ls.parse_or_refuse("bagel_ceg", "GENE\tHGNC_ID\tENTREZ_ID\n", record=_rec())


# --------------------------------------------------------------------------------------------------
# the measured discrepancy this contract exists for
# --------------------------------------------------------------------------------------------------
def test_bagel_joins_on_ENTREZ_ID_not_on_the_probes_lexical_pick():
    """MEASURED while building Step 5. The W0 probe lexically picks index 0 (`GENE`) because that column
    is named gene-ish; the real join key is index 2 (`ENTREZ_ID`), which is what the report card already
    joins on because gene_info is keyed by GeneID. So the probe CORROBORATES, it does not decide."""
    assert ls.LABEL_SOURCES["bagel_ceg"].key_column_index == 2
    assert ls.parse_bagel(BAGEL_TEXT) == {"16", "6122"}
    # parsing on the probe's lexical pick would give gene SYMBOLS, which join to nothing in gene_info
    assert ls.parse_bagel(BAGEL_TEXT, key_column_index=0) == {"AARS", "RPL3"}


# --------------------------------------------------------------------------------------------------
# reuse by identity — no forked parsers
# --------------------------------------------------------------------------------------------------
def test_the_sgd_parser_IS_the_fba_one_not_a_copy():
    from dna_decode.fba.essentiality_labels import parse_sgd_essential

    assert ls._sgd_parser() is parse_sgd_essential


def test_the_goodall_loader_IS_the_shipped_one_not_a_copy():
    loader = ls._goodall_loader()
    from essentiality_e3_learned import load_labels

    assert loader is load_labels


def test_no_source_reimplements_the_goodall_or_sgd_parse():
    """Checks for the LOGIC, not the words. The first version grepped for 'Unclear' / 'inviable' and
    failed on its own docstrings, which describe the reused parsers — a guard that cannot tell a
    description from an implementation."""
    src = open(ls.__file__, encoding="utf-8").read()
    # the Goodall parse reads boolean call columns by index; this module must not
    assert 'str(r[3])' not in src and 'str(r[4])' not in src and 'str(r[5])' not in src
    assert '== "true"' not in src.lower()
    # the SGD parse tests phenotype membership; this module must not
    assert '"inviable" in' not in src.lower()
    assert "openpyxl" not in src, "loading the Goodall workbook here would fork the reused loader"


# --------------------------------------------------------------------------------------------------
# registry shape
# --------------------------------------------------------------------------------------------------
def test_every_rung_has_a_label_plan_and_every_planned_source_is_registered():
    """A rung cannot exist without declared sources, and no plan may name an unregistered source.

    This caught a real Step-1/Step-5 mismatch: the ladder's human rung named a single source id while
    the human classes come from TWO files, so one id could not express it. RUNG_LABEL_PLAN now carries
    the per-rung shape and this checks both directions.
    """
    rung_keys = {r.key for r in tl.RUNGS}
    assert set(ls.RUNG_LABEL_PLAN) == rung_keys, set(ls.RUNG_LABEL_PLAN) ^ rung_keys
    for rk in rung_keys:
        srcs = ls.sources_for_rung(rk)
        assert srcs, rk
        for s in srcs:
            assert s in ls.LABEL_SOURCES, (rk, s)
    # and each rung's declared primary source is its plan's ESSENTIAL source
    for r in tl.RUNGS:
        assert r.label_source_id == ls.RUNG_LABEL_PLAN[r.key]["essential"], r.key


def test_the_three_class_sourcing_modes_are_all_real_and_disclosed():
    """The modes differ materially: an essential_plus_complement rung has an ANNOTATION-DEFINED negative
    class, whose size depends on the annotation rather than on a screen. That asymmetry must be visible,
    not averaged over."""
    modes = {p["mode"] for p in ls.RUNG_LABEL_PLAN.values()}
    assert modes == {"one_file_two_columns", "two_files", "essential_plus_complement"}
    assert ls.RUNG_LABEL_PLAN["human"]["mode"] == "two_files"
    assert ls.RUNG_LABEL_PLAN["ecoli"]["mode"] == "one_file_two_columns"
    src = open(ls.__file__, encoding="utf-8").read()
    assert "ANNOTATION-DEFINED" in src, "the complement-mode asymmetry must be stated in the source"


def test_every_source_declares_strain_technology_and_provenance():
    for sid, s in ls.LABEL_SOURCES.items():
        assert s.strain, sid
        assert s.technology in ("transposon_insertion", "deletion_collection", "crispr_ko"), sid
        assert s.recorded_in, sid
        assert isinstance(s.key_column_index, int), sid


def test_the_two_unreachable_sources_say_so_in_their_provenance():
    """Neither NTML nor PLOS GOLD has a resolvable URL recorded anywhere in-repo. The registry must say
    HYPOTHESIS rather than imply they are fetchable."""
    for sid in ("ntml_nebraska", "plos_gold_pa14"):
        assert ls.LABEL_SOURCES[sid].url is None, sid
        assert "HYPOTHESIS" in ls.LABEL_SOURCES[sid].recorded_in, sid


def test_the_unpinned_parsers_say_their_shape_is_unpinned():
    """A parser written against a guess must advertise that, or a later reader trusts it."""
    for fn in (ls.parse_ntml, ls.parse_plos_gold):
        assert "UNPINNED" in (fn.__doc__ or "").upper(), fn.__name__


def test_module_is_namespace_separate_from_the_fba_registry():
    src = open(ls.__file__, encoding="utf-8").read()
    assert "LABEL_WALLED" in src and "MODEL_WALLED" in src, \
        "the docstring must name the FBA walls it is distinguishing itself from"
    assert "MODEL-GENE-KEYED" in src


# --------------------------------------------------------------------------------------------------
# the two UNPINNED parsers — their shape is a guess, which is exactly why their BEHAVIOUR should be
# known. Until now only their docstrings were asserted, so neither had been run at all.
# --------------------------------------------------------------------------------------------------
NTML_TEXT = ("locus_tag\tcall\n"
             "SAUSA300_0001\tessential\n"
             "\n"
             "   \n"
             "SAUSA300_0002\tessential\n"
             "SAUSA300_0001\tduplicate row\n")


def test_parse_ntml_skips_the_header_and_blank_lines_and_dedups():
    assert ls.parse_ntml(NTML_TEXT) == {"SAUSA300_0001", "SAUSA300_0002"}
    # the column is a parameter, so a record that pins the key elsewhere is honoured rather than ignored
    assert ls.parse_ntml(NTML_TEXT, 1) == {"essential", "duplicate row"}
    # a key-column index past the end of every row yields nothing -- which `parse_or_refuse` turns into a
    # LabelParseError rather than an empty label set
    assert ls.parse_ntml(NTML_TEXT, 9) == set()


def test_parse_plos_gold_tolerates_real_xlsx_row_shapes():
    """An openpyxl `iter_rows(values_only=True)` sweep yields None for an empty ROW and None for an empty
    CELL, and numeric cells arrive as ints -- so a locus tag that happens to be numeric must be coerced,
    not dropped or left as an int that can never match a string key."""
    rows = [("header", "x"), None, ("PA14_00010", 5), ("", "blank key"), (None, "none cell"),
            (12345, "an int cell")]
    assert ls.parse_plos_gold(rows) == {"PA14_00010", "12345"}
    # a row shorter than the declared key column is skipped, not padded
    assert ls.parse_plos_gold([("h",), ("a",)], 1) == set()


def test_parse_bagel_tolerates_ragged_and_trailing_blank_lines():
    """Every real tsv ends with a newline and some carry short rows; a ragged row must be skipped rather
    than contribute an index error or a stray id."""
    assert ls.parse_bagel("GENE\tHGNC\tENTREZ\nAARS\tHGNC:20\t16\nshort\n\n") == {"16"}


# --------------------------------------------------------------------------------------------------
# which W0 record gates everything
# --------------------------------------------------------------------------------------------------
def test_latest_w0_record_is_None_on_a_directory_with_no_probe_artifact(tmp_path):
    """None is the input that makes `w0_status` refuse, so this is the path that keeps an unprobed host
    fail-closed rather than falling back to some other json in wiki/."""
    (tmp_path / "unrelated_other_artifact.json").write_text("{}", encoding="utf-8")
    assert ls.latest_w0_record(tmp_path) is None


def test_the_NEWEST_w0_record_wins_when_several_exist(tmp_path):
    """Probe artifacts are date-stamped and accrue, so the one that decides whether a parser may run is
    the last by name. A test on this is what makes 'rerun the probe to unblock a rung' true."""
    (tmp_path / "essentiality_ladder_w0_probe_2026-10-05.json").write_text(
        json.dumps({"sources": [{"source_id": "ntml_nebraska", "w0_verified": False,
                                 "reason": "older"}]}), encoding="utf-8")
    (tmp_path / "essentiality_ladder_w0_probe_2026-10-06.json").write_text(
        json.dumps({"sources": [{"source_id": "ntml_nebraska", "w0_verified": True,
                                 "reason": "newer"}]}), encoding="utf-8")
    rec = ls.latest_w0_record(tmp_path)
    ok, why = ls.w0_status("ntml_nebraska", rec)
    assert ok and why == "newer"


def test_a_record_with_NO_schema_block_still_parses_and_that_is_deliberate():
    """The key-column consistency check is skipped when the record lists no candidates -- it corroborates
    a declared index, it cannot manufacture one. Pinned because it is the one fail-OPEN path in an
    otherwise fail-closed gate: `w0_verified` is still required, so the record must have been verified by
    the probe; what is not required is that it enumerated key columns."""
    bare = {"sources": [{"source_id": "bagel_ceg", "w0_verified": True, "reason": "verified"}]}
    assert ls.parse_or_refuse("bagel_ceg", BAGEL_TEXT, record=bare) == {"16", "6122"}
    explicit_none = {"sources": [{"source_id": "bagel_ceg", "w0_verified": True, "reason": "r",
                                  "schema": {"key_column_candidates": None}}]}
    assert ls.parse_or_refuse("bagel_ceg", BAGEL_TEXT, record=explicit_none) == {"16", "6122"}


def test_real_w0_record_lets_bagel_through_if_present(tmp_path):
    """End-to-end against the committed record rather than a fixture, when one exists."""
    rec = ls.latest_w0_record()
    if rec is None:
        pytest.skip("no committed W0 record on this host")
    ok, why = ls.w0_status("bagel_ceg", rec)
    assert ok, why
    for sid in ("ntml_nebraska", "plos_gold_pa14", "sgd_phenotype"):
        ok2, _ = ls.w0_status(sid, rec)
        assert not ok2, f"{sid} should be unverified in the committed record"
