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


# --------------------------------------------------------------------------------------------------
# the remaining refusals in `verify` and `detect_header_row`
# --------------------------------------------------------------------------------------------------
def test_a_table_with_no_KEY_column_candidate_is_unverified():
    """The sibling of the no-label-column refusal, and the one that applies to tsv sources too: a table
    whose header names nothing gene-ish cannot be joined to an annotation, so no parser may be written
    against it. Without this the branch is unreachable in tests, because every other fixture here happens
    to carry a `GENE` column."""
    s = w0.summarize_table([["Score", "Value"], ["1", "2"]])
    assert s["key_column_candidates"] == [] and s["key_column_chosen"] is None
    assert s["n_keys"] == 0, "no key column -> no keys extracted, rather than column 0 by default"
    ok, why = w0.verify(s, {"kind": "tsv"})
    assert not ok and "no key-column candidate" in why


def test_header_detection_skips_None_rows_and_honours_its_window_and_min_cols():
    """`openpyxl` in read_only mode yields None for an entirely empty row, which is the real shape this
    has to survive. And the scan window is a bound, not a search: a header below it is reported as
    'cannot tell' rather than found, which is the honest failure for a sheet with a long preamble."""
    assert w0.detect_header_row([None, ["Gene", "Call"], ["a", "b"]]) == 1
    buried = [[None, None]] * 11 + [["Gene", "Call"], ["a", "b"]]
    assert w0.detect_header_row(buried) is None, "past the 10-row scan window -> cannot tell"
    assert w0.detect_header_row(buried, scan=20) == 11, "...and the window is a parameter"
    # min_cols is what separates a title from a header, so raising it must be able to reject a real one
    assert w0.detect_header_row([["Gene", "Call"], ["a", "b"]], min_cols=3) is None


# --------------------------------------------------------------------------------------------------
# IO + the sweep's promise never to abort
# --------------------------------------------------------------------------------------------------
def test_an_UNREADABLE_file_is_reported_as_read_failed_and_carries_no_schema(tmp_path, monkeypatch):
    """"report, never abort the sweep" -- an on-disk file that cannot be parsed must not take down the
    other five sources, and it must not leave a schema block behind that a parser could be written
    against. A present-but-broken file is a different state from an absent one and is named differently."""
    monkeypatch.setattr(w0, "CACHE", tmp_path)
    (tmp_path / "broken.xlsx").write_text("this is not a workbook", encoding="utf-8")
    rec = w0.probe_source("broken", {"local": "broken.xlsx", "kind": "xlsx", "organism": "O",
                                     "strain": "s", "key_space": "k",
                                     "technology": "transposon_insertion", "url": None,
                                     "recorded_in": "r"})
    assert rec["fetch_outcome"] == "read_failed"
    assert rec["w0_verified"] is False
    assert "BadZipFile" in rec["reason"] or "zip" in rec["reason"].lower()
    assert "schema" not in rec, "a file that could not be read has no pinned schema"


def test_read_rows_handles_a_GZIPPED_tsv(tmp_path):
    """`Homo_sapiens.gene_info.gz` is the real gzipped source in this set."""
    import gzip as _gz

    p = tmp_path / "t.tsv.gz"
    p.write_bytes(_gz.compress(b"GENE\tID\nRPL3\t6122\n"))
    rows, sheets = w0.read_rows(p, "tsv")
    assert rows == [["GENE", "ID"], ["RPL3", "6122"]]
    assert sheets == [], "sheet names are an xlsx-only concept"


# --------------------------------------------------------------------------------------------------
# the CLI exit contract — stated in the module docstring, previously unexercised
# --------------------------------------------------------------------------------------------------
def _one_absent_source(monkeypatch, tmp_path, url=None):
    monkeypatch.setattr(w0, "CACHE", tmp_path / "nonexistent_cache")
    monkeypatch.setattr(w0, "SOURCES", {
        "only": {"local": "absent.tsv", "kind": "tsv", "organism": "O", "strain": "s",
                 "key_space": "k", "technology": "crispr_ko", "url": url, "recorded_in": "r"}})


def test_the_probe_exits_ZERO_on_an_unverified_source_because_it_is_a_REPORT(tmp_path, monkeypatch):
    """Exit 0 always for the local pass. A report that failed the build on an unfetched file would make
    'absent_locally' look like an error, when it means 'not fetched yet'."""
    _one_absent_source(monkeypatch, tmp_path)
    assert w0.main(["--no-emit"]) == 0


def test_require_all_is_the_GATE_and_exits_one(tmp_path, monkeypatch):
    """The same run, opted into gate semantics. Both halves are needed: an exit code that never varies
    cannot gate anything, and one that always fails cannot be a report."""
    _one_absent_source(monkeypatch, tmp_path)
    assert w0.main(["--no-emit", "--require-all"]) == 1


def test_fetch_remote_records_that_it_did_NOT_fetch_and_never_invents_one(tmp_path, monkeypatch):
    """The flag exists so a reader can see the route was considered; it must not imply a fetch happened.
    Asserted on the EMITTED artifact, which also covers the emit path and the artifact's own honesty
    fields on a host where the committed artifact may be absent."""
    _one_absent_source(monkeypatch, tmp_path, url="http://example.invalid/labels.tsv")
    monkeypatch.setattr(w0, "WIKI", tmp_path)
    monkeypatch.setattr(w0, "ROOT", tmp_path)     # main prints the artifact path relative to ROOT
    assert w0.main(["--fetch-remote"]) == 0
    art = json.loads(next(tmp_path.glob("essentiality_ladder_w0_probe_*.json"))
                     .read_text(encoding="utf-8"))
    assert art["record"] == "essentiality-ladder-w0-probe-v1"
    assert art["n_verified"] == 0 and art["n_sources"] == 1
    assert "NOT IMPLEMENTED" in art["sources"][0]["remote_attempt"]
    assert "never invents one" in art["sources"][0]["remote_attempt"]
    assert any("HYPOTHESIS" in x for x in art["honest_limits"])


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
