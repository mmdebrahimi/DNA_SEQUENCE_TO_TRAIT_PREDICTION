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
