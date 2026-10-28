"""Guards for the human lentiMPRA loader.

The load-bearing test is `test_WTC11_REFUSES_rather_than_returning_its_3_percent_overlap`. Two of
three cell types join at 1.0000 and one joins at 0.0000, so a loader that returned "what it could"
would hand back a 3%-complete WTC11 cohort indistinguishable in shape from a complete one. That is
the failure this module exists to prevent, and it is also why `min_join_rate` is tested
BEHAVIOURALLY: relaxing the bar must not manufacture rows.

The real tables live on D: and are gitignored, so every data-dependent test skips when they are
absent rather than failing. The pure-logic tests always run.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.lentimpra import (  # noqa: E402
    CELL_TYPES,
    DEFAULT_DIR,
    EXPECTED_LEN,
    MIN_JOIN_RATE,
    S3_FILE,
    S6_FILE,
    LentiMpraDataError,
    LentiStats,
    load_activity,
    load_cell,
    load_design,
    published_incumbents,
    unparsed_rows,
    usable_cells,
)

pytest.importorskip("openpyxl", reason="openpyxl not installed")

HAVE_DATA = (DEFAULT_DIR / S3_FILE).is_file() and (DEFAULT_DIR / S6_FILE).is_file()
needs_data = pytest.mark.skipif(not HAVE_DATA, reason="lentiMPRA tables absent (D: / gitignored)")

#: Measured 2026-10-27 against the MD5-verified deposit.
JOINED = {"HepG2": 61463, "K562": 98332}
DESIGNED = {"HepG2": 164307, "K562": 243778, "WTC11": 75542}
#: The S3 header row MOVES between sheets of one workbook. Measured per sheet.
S3_HEADER_ROW = {"HepG2": 4, "K562": 2, "WTC11": 1}


# --------------------------------------------------------------------- pure logic

def test_the_join_bar_matches_the_sibling_loader_convention():
    """`expression.load_pairs` enforces 0.95 after the sigma70 key gap. Same bar here."""
    assert MIN_JOIN_RATE == 0.95


def test_an_unknown_cell_type_refuses_instead_of_guessing():
    with pytest.raises(LentiMpraDataError, match="unknown cell"):
        load_cell("Jurkat")


def test_the_cell_type_list_is_the_three_the_deposit_covers():
    assert CELL_TYPES == ("HepG2", "K562", "WTC11")


def test_stats_round_trip_through_their_dict():
    st = LentiStats(cell="HepG2", n_measured=10, n_joined=10, join_rate=1.0)
    d = st.as_dict()
    assert d["cell"] == "HepG2" and d["join_rate"] == 1.0


# --------------------------------------------------------------------- the refusal

@needs_data
def test_WTC11_REFUSES_rather_than_returning_its_3_percent_overlap():
    """THE load-bearing guard.

    WTC11's activity rows are keyed by Ensembl gene id and its design rows by element name, so the
    intersection is empty; Table S7 rescues only 699 of 23,064. A loader that returned the overlap
    would produce a cohort that LOOKS like a cohort.
    """
    with pytest.raises(LentiMpraDataError) as exc:
        load_cell("WTC11")
    msg = str(exc.value)
    assert "0.0000" in msg or "joins at 0" in msg
    assert "KEY problem" in msg, "the message must say WHY, or the next session re-derives it"
    assert "gene id" in msg.lower()


@needs_data
def test_relaxing_the_join_bar_does_NOT_manufacture_WTC11_rows():
    """Behavioural: the bar is a REPORT of a key mismatch, not the cause of the empty join."""
    with pytest.raises(LentiMpraDataError):
        load_cell("WTC11", min_join_rate=0.0001)
    # The design and activity sides each load fine on their own -- the failure is the JOIN.
    seqs, _ = load_design("WTC11")
    acts, _ = load_activity("WTC11")
    assert len(seqs) == DESIGNED["WTC11"]
    assert len(acts) == 23064
    assert not (set(seqs) & set(acts)), "the intersection must be empty, not merely small"


# --------------------------------------------------------------------- the usable cells

@needs_data
@pytest.mark.parametrize("cell", sorted(JOINED))
def test_the_two_usable_cells_join_completely(cell):
    recs, st = load_cell(cell)
    assert len(recs) == JOINED[cell]
    assert st.join_rate == pytest.approx(1.0, abs=1e-9), "every measured row must find its sequence"
    assert st.n_designed == DESIGNED[cell]
    assert st.n_measured == JOINED[cell]


@needs_data
@pytest.mark.parametrize("cell", sorted(JOINED))
def test_every_sequence_is_the_documented_uniform_length(cell):
    """Uniform length is why this panel needs no masked pooling, unlike the E. coli fragment rung."""
    recs, st = load_cell(cell)
    assert st.length_min == st.length_max == EXPECTED_LEN
    assert all(len(r.seq) == EXPECTED_LEN for r in recs[:500])


@needs_data
@pytest.mark.parametrize("cell", sorted(JOINED))
def test_the_published_ten_fold_assignment_comes_through(cell):
    """The authors' folds are what make this leakage-safe AND comparable to their benchmark."""
    recs, st = load_cell(cell)
    assert st.n_folds == 10
    assert set(st.fold_sizes) == set(range(1, 11))
    assert sum(st.fold_sizes.values()) == len(recs)
    # Near-equal by construction; no fold may be degenerate.
    assert min(st.fold_sizes.values()) > 0.5 * max(st.fold_sizes.values())


@needs_data
@pytest.mark.parametrize("cell", sorted(S3_HEADER_ROW))
def test_the_header_row_is_DISCOVERED_because_it_moves_between_sheets(cell):
    """A fixed min_row would read a 327-char design note as a column name on the HepG2 sheet."""
    _, st = load_design(cell)
    assert st.s3_header_row == S3_HEADER_ROW[cell]


@needs_data
def test_the_header_rows_are_genuinely_different_so_the_discovery_is_not_decorative():
    rows = {c: load_design(c)[1].s3_header_row for c in CELL_TYPES}
    assert len(set(rows.values())) > 1, (
        "if every sheet shared a header row, discovery would be ceremony: %r" % rows)


@needs_data
def test_the_sequence_column_RECORDS_that_adapters_are_included():
    """The 230 nt is 15 + 200 + 15. A silent trim would change what every number is about."""
    _, st = load_design("HepG2")
    col = st.sequence_column
    assert "230nt" in col and "adaptor" in col
    assert "200nt element" in col, "the biological element length must stay legible in the record"


@needs_data
def test_usable_cells_REPORTS_the_failure_rather_than_omitting_it():
    u = usable_cells()
    assert set(u) == set(CELL_TYPES), "a cell that fails must still appear"
    assert u["HepG2"]["usable"] and u["K562"]["usable"]
    assert not u["WTC11"]["usable"]
    assert "reason" in u["WTC11"]
    assert u["HepG2"]["n"] + u["K562"]["n"] == 159795


# --------------------------------------------------------------------- the incumbents

@needs_data
def test_the_published_metric_is_declared_UNSPECIFIED_rather_than_guessed():
    """A readout-type mismatch is not a comparison.

    The sheet gives method names and numbers and never says whether they are Pearson or Spearman.
    Guessing the likely answer would license a cross-readout comparison, which this repo has
    already been bitten by.
    """
    inc = published_incumbents()
    assert inc["metric"] == "UNSPECIFIED_IN_SOURCE"


@needs_data
def test_the_incumbent_table_carries_all_three_cells_and_ten_folds_each():
    inc = published_incumbents()["by_cell"]
    assert set(inc) == set(CELL_TYPES)
    for cell, methods in inc.items():
        assert len(methods) >= 5, "five benchmarked methods were published"
        for name, folds in methods.items():
            assert set(folds) == set(range(1, 11)), "%s/%s" % (cell, name)


@needs_data
def test_the_best_published_incumbent_is_named_per_cell():
    """Recorded so the arm's bar is written against a real incumbent, not an invented one."""
    import statistics

    inc = published_incumbents()["by_cell"]
    best = {c: max(ms.items(), key=lambda kv: statistics.fmean(kv[1].values()))[0]
            for c, ms in inc.items()}
    assert best["HepG2"] == "MPRALegNet"
    assert best["K562"] == "MPRALegNet"
    assert best["WTC11"] == "Enformer", "MPRALegNet does NOT win on WTC11; that asymmetry is real"


@needs_data
def test_an_absent_directory_refuses_with_the_fetch_command():
    with pytest.raises(LentiMpraDataError, match="fetch_lentimpra"):
        load_design("HepG2", data_dir=REPO / "no_such_dir_for_lentimpra")


# ---------------------------------------------------------------------------------------------
# The consumption rule, and the rows it discards.
#
# These run WITHOUT the D: volume: they build a synthetic Table S6 with openpyxl. That matters
# beyond convenience. The metric refusal rested on "the sheet never states a coefficient", which
# was a claim about the PARSED SUBSET being read as a claim about the sheet. Making the discarded
# rows reachable is what turns that from prose into something checkable, so the check itself must
# not depend on the volume being attached.
# ---------------------------------------------------------------------------------------------

#: A caption of exactly the shape that defeats the parser: it NAMES a coefficient, and it is
#: discarded. Everything these tests assert turns on this row being invisible to one accessor and
#: visible to the other.
_SYNTH_CAPTION = ("Supplementary Table 6. Per-fold benchmarking performance "
                  "(Pearson r) across methods.")


@pytest.fixture
def synth_s6(tmp_path):
    """A minimal `Benchmarking results` workbook carrying every row shape the parser meets."""
    import openpyxl

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Benchmarking results"
    for row in [
        [_SYNTH_CAPTION, None, None, None],
        [None, None, None, None],
        ["HepG2", None, None, None],
        ["test fold", "MPRALegNet", "Enformer", "Legnet"],
        [1, 0.8301, 0.7702, 0.8104],
        [2, 0.8365, 0.7811, 0.8150],
        ["Note: fold 3 withheld pending re-run.", None, None, None],
        ["K562", None, None, None],
        ["test fold", "MPRALegNet", "Enformer", "Legnet"],
        [1, 0.7905, 0.7301, 0.7788],
        [2, 0.7951, 0.7355, 0.7802],
    ]:
        ws.append(row)
    wb.save(tmp_path / S6_FILE)
    return tmp_path


def test_the_consumption_rule_is_defined_in_exactly_one_place():
    """Source-shape, because value equality cannot see a reintroduced duplicate.

    Asserted on the AST rather than on substrings: the corrected `decisive_next_check` prose has to
    QUOTE the header (`test fold | <method names>`) in order to explain why reading it cannot
    settle the metric, so a substring check would fail on the correction itself. Only a literal in
    an executable COMPARISON counts as a predicate; one inside a docstring or a returned dict is
    prose and is exempt structurally.
    """
    import ast

    from dna_decode.glm import lentimpra

    src = Path(lentimpra.__file__).read_text(encoding="utf-8")
    fns = {n.name: n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef)}

    def comparison_operands(fn):
        out = []
        for node in ast.walk(fn):
            if not isinstance(node, ast.Compare):
                continue
            for part in [node.left] + list(node.comparators):
                for sub in ast.walk(part):
                    if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                        out.append(sub.value)
                    elif isinstance(sub, ast.Name):
                        out.append("name:" + sub.id)
        return out

    walk = comparison_operands(fns["_walk_rows"])
    assert "test fold" in walk and "name:CELL_TYPES" in walk, "_walk_rows must own both predicates"

    for consumer in ("published_incumbents", "unparsed_rows"):
        ops = comparison_operands(fns[consumer])
        assert "test fold" not in ops, "%s re-tests the header string" % consumer
        assert "name:CELL_TYPES" not in ops, "%s re-tests CELL_TYPES" % consumer

    assert src.count('first.lower() == "test fold"') == 1
    assert src.count("first in CELL_TYPES") == 1


def test_extracting_the_walk_did_not_change_what_the_incumbents_are(synth_s6):
    """Equivalence proof for a rewrite of the function that decides the published incumbents.

    The expected payload is written out literally rather than snapshotted, so the test carries its
    own reference instead of depending on a file some earlier run produced. Fold keys are INTS --
    a JSON round-trip coerces them to strings, which makes a snapshot comparison fail for a reason
    that has nothing to do with the refactor.
    """
    got = published_incumbents(data_dir=synth_s6)
    assert got["by_cell"] == {
        "HepG2": {"MPRALegNet": {1: 0.8301, 2: 0.8365},
                  "Enformer": {1: 0.7702, 2: 0.7811},
                  "Legnet": {1: 0.8104, 2: 0.8150}},
        "K562": {"MPRALegNet": {1: 0.7905, 2: 0.7951},
                 "Enformer": {1: 0.7301, 2: 0.7355},
                 "Legnet": {1: 0.7788, 2: 0.7802}},
    }
    assert got["metric"] == "UNSPECIFIED_IN_SOURCE", "the refusal must survive the refactor"


def test_unparsed_rows_returns_the_caption_the_parse_discards(synth_s6):
    """The load-bearing pair. A reader that returns everything, or nothing, fails one half.

    The caption names a coefficient. `published_incumbents` cannot see it, `unparsed_rows` must --
    which is precisely why "the sheet never states a metric" could not be concluded from the parse.
    """
    skipped = unparsed_rows(data_dir=synth_s6)
    first_cells = [str(r["cells"][0]) for r in skipped]

    assert _SYNTH_CAPTION in first_cells, "the discarded caption must be recoverable"
    assert "Pearson" not in str(published_incumbents(data_dir=synth_s6)["by_cell"]), (
        "the parse must NOT see it -- if it did, the two accessors would not be complementary")

    # Returns the skipped rows and ONLY the skipped rows: no cell-type marker, header or fold row.
    assert not any(c in CELL_TYPES for c in first_cells), "a consumed cell-type row leaked through"
    assert not any(c.lower() == "test fold" for c in first_cells), "a consumed header leaked through"
    assert not any(c in ("1", "2") for c in first_cells), "a consumed fold row leaked through"

    # Blank rows carry no information and are dropped; 1-based indices match what a human sees.
    assert all(r["cells"] for r in skipped), "a wholly empty row is not a finding"
    assert [r["row"] for r in skipped] == [1, 7], "1-based row indices, as a spreadsheet shows them"


def test_unparsed_rows_refuses_an_absent_workbook_rather_than_reporting_none(tmp_path):
    """An empty list would read as "there are no unparsed rows" when the truth is "could not look".

    That conflation is the whole subject of this correction, so the accessor raises like its sibling.
    """
    with pytest.raises(LentiMpraDataError, match="fetch_lentimpra"):
        unparsed_rows(data_dir=tmp_path / "no_such_dir")
