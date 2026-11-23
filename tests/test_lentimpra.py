"""Guards for the human lentiMPRA loader.

All THREE cell types now join at 1.0000. WTC11 previously joined at 0.0000 and refused; corrected
2026-11-23 once the mechanism was measured rather than restated -- S6 ids are bare for every cell,
and S3 WTC11 is the only sheet that suffixes the FORWARD strand (`seq1000_F`), so identity alone
could never match it. See the loader docstring for the three measurements that settle the strand.

THE REFUSAL DISCIPLINE IS NOT RELAXED, it is re-pointed. The load-bearing tests are now:

  * `test_the_forward_suffix_is_LOAD_BEARING_and_the_bar_still_bites` -- with the suffix set to a
    value that matches nothing, WTC11 falls back to 0.0000 and REFUSES. One test shows the rule is
    doing real work AND that `min_join_rate` still catches a key gap, so neither claim rests on a
    guard that might have quietly stopped firing.
  * `test_a_bare_id_pairs_with_the_FORWARD_strand_not_the_reverse` -- `_F` and `_R` differ for every
    measured row while the join rate reads 1.0000 either way, so a rate assertion cannot catch a
    wrong-strand cohort. This pins the sequence.
  * `test_the_two_previously_shipped_cells_are_UNTOUCHED_by_the_new_rule` -- the fix may only ADD
    rows; HepG2/K562 must still join by identity with zero suffix hits.

`min_join_rate` stays tested BEHAVIOURALLY: relaxing the bar must not manufacture rows.

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

def _resolve_data_dir():
    """`DEFAULT_DIR`, else the repo-local copy -- mirroring `glm_human_activity_oracle.LOCAL_DIR`.

    This existed as a bug rather than a decision: the gate below read DEFAULT_DIR
    (`D:/dna_decode_cache/mpra_human`), D: is the empty card-reader slot on this host, and the
    tables actually sit at `data/mpra_human`. So EVERY data-dependent test in this module skipped,
    and the loader's load-bearing guards were not guarding anything on the machine the work happens
    on. The scripts already pass an explicit local dir for exactly this reason; the tests now do too.
    """
    for cand in (DEFAULT_DIR, REPO / "data" / "mpra_human"):
        if (cand / S3_FILE).is_file() and (cand / S6_FILE).is_file():
            return cand
    return DEFAULT_DIR


DATA_DIR = _resolve_data_dir()
HAVE_DATA = (DATA_DIR / S3_FILE).is_file() and (DATA_DIR / S6_FILE).is_file()
needs_data = pytest.mark.skipif(not HAVE_DATA, reason="lentiMPRA tables absent (gitignored)")

#: Measured 2026-10-27 against the MD5-verified deposit; WTC11 added 2026-11-23 with the name rule.
JOINED = {"HepG2": 61463, "K562": 98332, "WTC11": 23064}
#: The two that join on a bare id alone. Kept separate from JOINED so the regression guard below
#: cannot silently start covering WTC11 and report "untouched" about the cell the rule changed.
IDENTITY_CELLS = ("HepG2", "K562")
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
def test_identity_alone_still_cannot_join_WTC11():
    """The PREMISE of the name rule, kept as a standing measurement rather than prose.

    If this ever starts passing on bare ids, the deposit changed its naming and the rule below is no
    longer describing the substrate -- which is the thing the old, misstated docstring got wrong by
    asserting a mechanism nobody re-measured.
    """
    seqs, _ = load_design("WTC11", DATA_DIR)
    acts, _ = load_activity("WTC11", DATA_DIR)
    assert len(seqs) == DESIGNED["WTC11"]
    assert len(acts) == JOINED["WTC11"]
    assert not (set(seqs) & set(acts)), "the BARE intersection must be empty; that is why _F exists"
    # ...and it is empty for the whole measured set, not just its 19% ENSG part.
    assert sum(1 for k in acts if k.startswith("ENSG")) == 4397
    assert sum(1 for k in acts if k.startswith("seq")) == 18667


@needs_data
def test_the_forward_suffix_is_LOAD_BEARING_and_the_bar_still_bites(monkeypatch):
    """ONE test for two claims that must not be allowed to rest on each other.

    Point the suffix at something the sheet does not contain: WTC11 must fall back to 0.0000 and
    REFUSE. That shows the rule is doing real work (not coincidentally agreeing with identity) and
    that `min_join_rate` still catches a key gap now that the headline cell no longer trips it.
    """
    import dna_decode.glm.lentimpra as mod

    # No clear_cache: the cache holds PARSED SHEETS keyed by (kind, cell, dir), and the suffix is
    # read at join time in load_cell, which is not cached. Clearing would force a fresh 56 MB parse
    # per test for no behavioural difference.
    monkeypatch.setattr(mod, "ORIENTATION_SUFFIX", "_NOT_A_REAL_SUFFIX")
    with pytest.raises(LentiMpraDataError) as exc:
        mod.load_cell("WTC11", DATA_DIR)
    msg = str(exc.value)
    assert "0.0000" in msg or "joins at 0" in msg
    assert "KEY problem" in msg, "the message must say WHY, or the next session re-derives it"
    assert "do not relax the bar" in msg


@needs_data
def test_relaxing_the_join_bar_does_NOT_manufacture_rows(monkeypatch):
    """Behavioural: the bar REPORTS a key gap, it does not cause it. Relaxing it adds nothing."""
    import dna_decode.glm.lentimpra as mod

    monkeypatch.setattr(mod, "ORIENTATION_SUFFIX", "_NOT_A_REAL_SUFFIX")
    with pytest.raises(LentiMpraDataError):
        mod.load_cell("WTC11", DATA_DIR, min_join_rate=0.0001)


@needs_data
def test_a_bare_id_pairs_with_the_FORWARD_strand_not_the_reverse():
    """THE load-bearing guard now.

    `_F` and `_R` are different sequences for all 23,064 measured rows, so the join rate is 1.0000
    whichever is chosen and CANNOT detect a wrong-strand cohort. Settled on sequence identity: the
    shared promoters the two shipped cells pair by identity equal WTC11's `_F` for 4,397 of 4,397
    and `_R` for 0. This pins that result against the deployed loader.
    """
    recs, st = load_cell("WTC11", DATA_DIR)
    seqs, _ = load_design("WTC11", DATA_DIR)
    assert st.name_rule == "identity_then_forward_suffix"
    assert st.n_joined_forward == JOINED["WTC11"], "every WTC11 row comes via the suffix"
    assert st.n_joined_identity == 0

    checked = 0
    for r in recs[:400]:
        f, rev = seqs.get(r.name + "_F"), seqs.get(r.name + "_R")
        assert f is not None and rev is not None and f != rev
        assert r.seq == f, f"{r.name} paired with the REVERSE strand"
        assert r.seq != rev
        checked += 1
    assert checked == 400


@needs_data
def test_the_forward_and_reverse_inserts_are_reverse_complements():
    """Corroborates that _F/_R are orientations of ONE element, not unrelated oligos -- and pins the
    adaptor coordinates, since comparing the full 230 nt finds zero revcomp pairs (the 15 nt
    adaptors are constant on both strands) and would read as a refutation."""
    seqs, _ = load_design("WTC11", DATA_DIR)
    comp = str.maketrans("ACGTN", "TGCAN")
    pairs = 0
    for name in list(seqs)[:2000]:
        if not name.endswith("_F"):
            continue
        rev = seqs.get(name[:-2] + "_R")
        if rev is None:
            continue
        f_ins, r_ins = seqs[name][15:215], rev[15:215]
        assert f_ins.translate(comp)[::-1] == r_ins, f"{name} inserts are not revcomp partners"
        pairs += 1
    assert pairs > 100, "the sample must actually contain F/R pairs"


@needs_data
@pytest.mark.parametrize("cell", IDENTITY_CELLS)
def test_the_two_previously_shipped_cells_are_UNTOUCHED_by_the_new_rule(cell):
    """The fix may only ADD rows. If a suffix hit appears here, the lookup is reaching rows that
    identity already had -- which is how a 'rescue' silently rewrites a shipped cohort."""
    _, st = load_cell(cell, DATA_DIR)
    assert st.name_rule == "identity"
    assert st.n_joined_forward == 0
    assert st.n_joined_identity == JOINED[cell]


@needs_data
def test_an_ambiguous_id_REFUSES_rather_than_choosing_a_strand(monkeypatch):
    """If an id resolves both bare and `_F` to DIFFERENT sequences, neither is knowably the measured
    construct and BOTH score a complete join. Synthetic: the real deposit has no such id, so this
    guard would never fire on the data and is the one case that must be constructed to be tested."""
    import dna_decode.glm.lentimpra as mod

    seqs = {"dup1": "A" * 230, "dup1_F": "C" * 230}
    monkeypatch.setattr(mod, "load_design", lambda c, d=None: (seqs, mod.LentiStats(cell=c)))
    monkeypatch.setattr(mod, "load_activity",
                        lambda c, d=None: ({"dup1": (1.0, 1)}, mod.LentiStats(cell=c, n_measured=1)))
    with pytest.raises(LentiMpraDataError, match="resolves BOTH"):
        mod.load_cell("WTC11", DATA_DIR)


@needs_data
def test_the_total_usable_cohort_is_the_three_cell_sum():
    """Pins the docstring's 182,859 against the loader, so the two cannot drift apart."""
    assert sum(JOINED.values()) == 182859
    assert sum(v["n"] for v in usable_cells(DATA_DIR).values() if v["usable"]) == 182859


# --------------------------------------------------------------------- the usable cells

@needs_data
@pytest.mark.parametrize("cell", sorted(JOINED))
def test_all_three_cells_join_completely(cell):
    recs, st = load_cell(cell, DATA_DIR)
    assert len(recs) == JOINED[cell]
    assert st.join_rate == pytest.approx(1.0, abs=1e-9), "every measured row must find its sequence"
    assert st.n_designed == DESIGNED[cell]
    assert st.n_measured == JOINED[cell]
    # The two counters must account for every row, or one rule is being credited for the other's.
    assert st.n_joined_identity + st.n_joined_forward == JOINED[cell]


@needs_data
@pytest.mark.parametrize("cell", sorted(JOINED))
def test_every_sequence_is_the_documented_uniform_length(cell):
    """Uniform length is why this panel needs no masked pooling, unlike the E. coli fragment rung."""
    recs, st = load_cell(cell, DATA_DIR)
    assert st.length_min == st.length_max == EXPECTED_LEN
    assert all(len(r.seq) == EXPECTED_LEN for r in recs[:500])


@needs_data
@pytest.mark.parametrize("cell", sorted(JOINED))
def test_the_published_ten_fold_assignment_comes_through(cell):
    """The authors' folds are what make this leakage-safe AND comparable to their benchmark."""
    recs, st = load_cell(cell, DATA_DIR)
    assert st.n_folds == 10
    assert set(st.fold_sizes) == set(range(1, 11))
    assert sum(st.fold_sizes.values()) == len(recs)
    # Near-equal by construction; no fold may be degenerate.
    assert min(st.fold_sizes.values()) > 0.5 * max(st.fold_sizes.values())


@needs_data
@pytest.mark.parametrize("cell", sorted(S3_HEADER_ROW))
def test_the_header_row_is_DISCOVERED_because_it_moves_between_sheets(cell):
    """A fixed min_row would read a 327-char design note as a column name on the HepG2 sheet."""
    _, st = load_design(cell, DATA_DIR)
    assert st.s3_header_row == S3_HEADER_ROW[cell]


@needs_data
def test_the_header_rows_are_genuinely_different_so_the_discovery_is_not_decorative():
    rows = {c: load_design(c, DATA_DIR)[1].s3_header_row for c in CELL_TYPES}
    assert len(set(rows.values())) > 1, (
        "if every sheet shared a header row, discovery would be ceremony: %r" % rows)


@needs_data
def test_the_sequence_column_RECORDS_that_adapters_are_included():
    """The 230 nt is 15 + 200 + 15. A silent trim would change what every number is about."""
    _, st = load_design("HepG2", DATA_DIR)
    col = st.sequence_column
    assert "230nt" in col and "adaptor" in col
    assert "200nt element" in col, "the biological element length must stay legible in the record"


@needs_data
def test_usable_cells_reports_all_three_as_usable():
    u = usable_cells(DATA_DIR)
    assert set(u) == set(CELL_TYPES), "every cell must appear, usable or not"
    assert all(u[c]["usable"] for c in CELL_TYPES)
    assert sum(u[c]["n"] for c in CELL_TYPES) == 182859


@needs_data
def test_usable_cells_REPORTS_a_failure_rather_than_omitting_it(monkeypatch):
    """The discipline this function exists for, kept testable now that no real cell fails.

    Previously this rode on WTC11 genuinely refusing. A guard that depends on a real failure stops
    being a guard the moment the failure is fixed -- and would then pass vacuously -- so the failure
    is induced instead, by pointing the suffix at nothing.
    """
    import dna_decode.glm.lentimpra as mod

    monkeypatch.setattr(mod, "ORIENTATION_SUFFIX", "_NOT_A_REAL_SUFFIX")
    u = mod.usable_cells(DATA_DIR)
    assert set(u) == set(mod.CELL_TYPES), "a cell that fails must still APPEAR"
    assert u["HepG2"]["usable"] and u["K562"]["usable"], "the other two must be unaffected"
    assert not u["WTC11"]["usable"]
    assert "reason" in u["WTC11"] and "KEY problem" in u["WTC11"]["reason"]
    assert "n" not in u["WTC11"], "a failing cell must not report a count"


# --------------------------------------------------------------------- the incumbents

@needs_data
def test_the_published_metric_is_declared_UNSPECIFIED_rather_than_guessed():
    """A readout-type mismatch is not a comparison.

    The sheet gives method names and numbers and never says whether they are Pearson or Spearman.
    Guessing the likely answer would license a cross-readout comparison, which this repo has
    already been bitten by.
    """
    inc = published_incumbents(DATA_DIR)
    assert inc["metric"] == "UNSPECIFIED_IN_SOURCE"


@needs_data
def test_the_incumbent_table_carries_all_three_cells_and_ten_folds_each():
    inc = published_incumbents(DATA_DIR)["by_cell"]
    assert set(inc) == set(CELL_TYPES)
    for cell, methods in inc.items():
        assert len(methods) >= 5, "five benchmarked methods were published"
        for name, folds in methods.items():
            assert set(folds) == set(range(1, 11)), "%s/%s" % (cell, name)


@needs_data
def test_the_best_published_incumbent_is_named_per_cell():
    """Recorded so the arm's bar is written against a real incumbent, not an invented one."""
    import statistics

    inc = published_incumbents(DATA_DIR)["by_cell"]
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
