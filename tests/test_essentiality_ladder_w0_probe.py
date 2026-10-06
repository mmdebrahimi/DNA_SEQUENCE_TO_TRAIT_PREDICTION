"""Guards for the W0 schema probe — the gate that exists because two prose-written parsers failed silently.

Both failures produced a NUMBER rather than an error, which is why a probe and not a try/except is the
answer:
  * `uniprot_ecoli_p1.tsv` was assumed a product annotation; it is 501 rows of `Disruption phenotype`.
  * `goodall_TableS1_essential.xlsx` row 0 is a TITLE; reading it as the header mis-joined 3936/372
    against the committed 351/3432.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
import sys  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
import essentiality_ladder_w0_probe as w0  # noqa: E402


GOODALL_LIKE = [
    ["Table S1. Essentiality classification for genes of the TL data", None, None, None, None, None],
    ["Gene", "Insertion Index Score", "Log Likelihood Ratio", "Essential", "Non-essential", "Unclear"],
    ["thrL", "0.39", "31.4", "False", "True", "False"],
    ["thrA", "0.22", "17.2", "False", "True", "False"],
    ["thrB", "0.26", "21.2", "True", "False", "False"],
]
UNIPROT_LIKE = [
    ["Gene Names (primary)", "Gene Names (ordered locus)", "Disruption phenotype"],
    ["ndh", "b1109 JW1095", "DISRUPTION PHENOTYPE: retains activity ..."],
]
PLAIN_TSV = [
    ["GENE", "ENTREZ_ID"],
    ["RPL3", "6122"],
    ["RPL4", "6124"],
]


# --------------------------------------------------------------------------------------------------
# header detection — THE Goodall fix
# --------------------------------------------------------------------------------------------------
def test_a_title_row_is_not_mistaken_for_a_header():
    assert w0.detect_header_row(GOODALL_LIKE) == 1
    assert w0.detect_header_row(PLAIN_TSV) == 0


def test_reading_row_zero_as_header_WOULD_have_produced_a_different_label_count():
    """NON-VACUITY: proves the fix is load-bearing by reproducing the mis-parse it prevents.

    With the true header (row 1) there are 3 data rows. Treating the title as the header makes row 1
    itself a data row -- 4 'labels', and the 'Essential' column is no longer where it is thought to be.
    A guard that could not show the two differ would be testing nothing.
    """
    right = w0.summarize_table(GOODALL_LIKE)
    assert right["header_row_index"] == 1 and right["n_data_rows"] == 3

    # simulate the bug: force row 0 to be the header
    wrong_data = [r for r in GOODALL_LIKE[1:] if any(c not in (None, "") for c in r)]
    assert len(wrong_data) == 4
    assert len(wrong_data) != right["n_data_rows"], "the mis-parse must differ, else nothing is guarded"


def test_the_label_column_detector_is_FOOLED_by_a_descriptive_title():
    """A measured weakness, recorded rather than assumed away.

    I first asserted that a title row "yields NO label column -- the tell". That is FALSE on the real
    file: the title is "Table S1. ESSENTIALITY classification for genes of the TL data", and the lexical
    detector matches `essential` inside `Essentiality`, so it happily nominates column 0 as a label
    column. The consequence matters for how this probe is trusted: `verify()`'s label-column check would
    NOT have rejected a row-0-as-header read. `detect_header_row` is the real protection; the
    label-column check is corroboration only, and must never be relied on as the header guard.
    """
    hdr = w0.normalise_header(GOODALL_LIKE[0])
    assert w0.label_column_candidates(hdr) == [0], \
        "if this ever becomes [] the detector changed and the recorded weakness should be re-read"
    # ...and the header detector still gets it right, which is the part that is load-bearing
    assert w0.detect_header_row(GOODALL_LIKE) == 1


def test_a_title_only_sheet_refuses_rather_than_guessing_row_zero():
    """'Cannot tell' is the honest answer; defaulting to 0 is how the Goodall bug happened."""
    assert w0.detect_header_row([["Just a title", None, None]]) is None
    s = w0.summarize_table([["Just a title", None, None]])
    assert s["header_row_index"] is None and "problem" in s
    assert not w0.verify(s, {"kind": "xlsx"})[0]


# --------------------------------------------------------------------------------------------------
# wrong-shape detection — THE uniprot trap
# --------------------------------------------------------------------------------------------------
def test_a_phenotype_description_table_is_rejected_as_wrong_shape():
    s = w0.summarize_table(UNIPROT_LIKE)
    ok, why = w0.verify(s, {"kind": "tsv"})
    assert not ok
    assert "disruption phenotype" in why.lower()


def test_wrong_shape_fires_on_the_header_not_on_row_count():
    """501 rows was a symptom; the header is the diagnosis. A 50,000-row version must still be caught."""
    big = [UNIPROT_LIKE[0]] + [UNIPROT_LIKE[1]] * 5000
    assert not w0.verify(w0.summarize_table(big), {"kind": "tsv"})[0]


# --------------------------------------------------------------------------------------------------
# summary fields parsers depend on
# --------------------------------------------------------------------------------------------------
def test_duplicate_keys_are_counted_not_silently_deduped():
    rows = [["GENE", "ID"], ["A", "1"], ["B", "2"], ["A", "3"]]
    s = w0.summarize_table(rows)
    assert s["n_keys"] == 3 and s["n_keys_distinct"] == 2 and s["n_duplicate_keys"] == 1


def test_key_columns_are_reported_by_INDEX_so_duplicate_names_cannot_collapse_them():
    rows = [["gene", "gene"], ["a", "b"]]
    assert w0.key_column_candidates(w0.normalise_header(rows[0])) == [0, 1]


def test_zero_data_rows_is_unverified():
    s = w0.summarize_table([["GENE", "ID"]])
    assert s["n_data_rows"] == 0
    assert not w0.verify(s, {"kind": "tsv"})[0]


def test_an_xlsx_label_table_without_a_call_column_is_unverified():
    """A label table must carry a call column; a key-only sheet cannot label anything."""
    s = w0.summarize_table([["Gene", "Score"], ["thrA", "0.2"]])
    ok, why = w0.verify(s, {"kind": "xlsx"})
    assert not ok and "label-column" in why


# --------------------------------------------------------------------------------------------------
# sweep behaviour
# --------------------------------------------------------------------------------------------------
def test_an_absent_source_is_unverified_and_does_not_abort_the_sweep():
    rec = w0.probe_source("nope", {"local": "definitely_absent_xyz.tsv", "kind": "tsv",
                                   "organism": "X", "strain": "s", "key_space": "k",
                                   "technology": "transposon_insertion", "url": None,
                                   "recorded_in": "r"})
    assert rec["w0_verified"] is False
    assert rec["fetch_outcome"] == "absent_locally"
    assert "HYPOTHESIS" in rec["reason"], "an absent file with no URL must say reachability is unproven"


def test_self_check_passes_offline():
    assert w0._self_check() == 0


def test_every_declared_source_names_where_it_is_recorded():
    """A source with no in-repo provenance pointer is a remembered claim."""
    for sid, d in w0.SOURCES.items():
        assert d["recorded_in"], sid
        assert d["strain"], sid
        assert d["technology"] in ("transposon_insertion", "deletion_collection", "crispr_ko"), sid


def test_the_committed_artifact_records_the_three_new_rungs_as_unverified():
    """The honest state of the plan: the 3 sources the ladder NEEDS are exactly the 3 not yet on disk."""
    hits = sorted(ROOT.glob("wiki/essentiality_ladder_w0_probe_*.json"))
    if not hits:
        pytest.skip("probe artifact not emitted on this host")
    art = json.loads(hits[-1].read_text(encoding="utf-8"))
    assert art["record"] == "essentiality-ladder-w0-probe-v1"
    unver = {s["source_id"] for s in art["sources"] if not s["w0_verified"]}
    assert {"sgd_phenotype", "ntml_nebraska", "plos_gold_pa14"} <= unver, unver
    # and the verified ones are the two rungs that already score
    ver = {s["source_id"] for s in art["sources"] if s["w0_verified"]}
    assert {"goodall_tradis", "bagel_ceg", "bagel_neg"} <= ver, ver
    assert art["consumed_by"].startswith("dna_decode/essentiality/label_sources.py")
    assert any("HYPOTHESIS" in x for x in art["honest_limits"])


def test_the_goodall_record_pins_header_row_one_on_the_REAL_file():
    """Guards the live fact, not just the fixture."""
    hits = sorted(ROOT.glob("wiki/essentiality_ladder_w0_probe_*.json"))
    if not hits:
        pytest.skip("probe artifact not emitted on this host")
    art = json.loads(hits[-1].read_text(encoding="utf-8"))
    g = next((s for s in art["sources"] if s["source_id"] == "goodall_tradis"), None)
    if not g or not g.get("schema"):
        pytest.skip("Goodall workbook not present on this host")
    assert g["schema"]["header_row_index"] == 1
    assert g["schema"]["title_rows_skipped"] == 1
