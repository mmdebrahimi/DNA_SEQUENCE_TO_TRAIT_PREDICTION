"""Tests for the positional frame features.

Three things here carry weight. The **sign-convention** tests are separate per strand, because a
minus-strand sign error is the failure that produces plausible numbers rather than an error -- and
the minus strand's "upstream" runs the opposite way along the genome from the plus strand's. The
**arithmetic reconciliation** test exists because an earlier abstention arm in this repo turned true
positives into an apparent bias through a quietly reduced denominator, and that defect was caught by
the arithmetic rather than by reading code. And the **union-vs-single-window** test pins the
definition the committed 8,239 figure was actually measured with.
"""

from __future__ import annotations

import pathlib

import pytest

from dna_decode.glm.corpus_window import (
    CdsRow,
    corpus_windows,
    window_coverage,
)
from dna_decode.glm.frame_features import (
    FEATURE_NAMES,
    IN_WINDOW_ANY_MIN_OVERLAP_BASES,
    FrameFeatureError,
    FrameFeatureTable,
    frame_features,
)
from dna_decode.glm.tile_locate import LocateResult

REPO = pathlib.Path(__file__).resolve().parents[1]
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
SUBSTRATE = REPO / "data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz"

GENOME_LEN = 10_000
TL = 10

#: One plus gene (window 4800..5000) and one minus gene (window 3000..3200).
CDS = [CdsRow(5000, 6000, "+", "p"), CdsRow(2000, 3000, "-", "m")]


def _mask():
    return window_coverage(corpus_windows(CDS, GENOME_LEN), GENOME_LEN)


def _table(entries: dict[str, LocateResult]):
    return frame_features(entries, _mask(), CDS, tile_length=TL)


def _row(seq: str, pos0: int, strand: str):
    t = _table({seq: LocateResult(1, pos0, strand)})
    return dict(zip(t.names, t.matrix[0]))


def test_the_synthetic_setup_places_the_windows_where_the_tests_assume():
    """Guard the fixture itself, so a window-rule change does not silently void every case."""
    w = corpus_windows(CDS, GENOME_LEN)
    assert sorted((x.start0, x.end0, x.strand) for x in w) == [
        (3000, 3200, "-"), (4800, 5000, "+")
    ]


# ------------------------------------------------------------------------ overlap granularities


def test_a_tile_fully_inside_a_window_is_any_and_full():
    r = _row("A", 4900, "+")
    assert r["in_window_any"] == 1.0
    assert r["in_window_full"] == 1.0
    assert r["overlap_fraction"] == 1.0


def test_a_tile_straddling_a_window_edge_has_a_fraction_strictly_between_zero_and_one():
    r = _row("A", 4995, "+")          # covers 4995..5000 of its 4995..5005 span
    assert 0.0 < r["overlap_fraction"] < 1.0
    assert r["overlap_fraction"] == pytest.approx(0.5)
    assert r["in_window_any"] == 1.0
    assert r["in_window_full"] == 0.0, "a straddling tile is not full"


def test_a_tile_far_from_any_window_is_out_on_every_granularity():
    r = _row("A", 100, "+")
    assert r["in_window_any"] == 0.0
    assert r["in_window_full"] == 0.0
    assert r["overlap_fraction"] == 0.0


def test_in_window_any_fires_at_exactly_one_overlapping_base():
    """The threshold is NAMED, so an epsilon cannot move the count silently."""
    assert IN_WINDOW_ANY_MIN_OVERLAP_BASES == 1
    r = _row("A", 4999, "+")          # exactly one base (4999) inside 4800..5000
    assert r["overlap_fraction"] == pytest.approx(1 / TL)
    assert r["in_window_any"] == 1.0
    assert _row("A", 5000, "+")["in_window_any"] == 0.0, "zero overlap must not fire"


def test_in_window_full_is_UNION_coverage_not_inside_a_single_window():
    """Two overlapping windows jointly cover a tile that neither of them contains.

    This is the definition `scripts/frame_plan_feasibility_probe.py` measured 8,239 with
    (`c == TILE_LEN` against one union mask). The plan's prose says "entirely inside one", which is
    a different predicate; this test is where the two are kept from drifting under one name.
    """
    mask = window_coverage([(100, 200, "+", "a"), (195, 300, "+", "b")], 1000)
    t = frame_features({"S": LocateResult(1, 190, "+")}, mask, CDS, tile_length=20)
    r = dict(zip(t.names, t.matrix[0]))
    assert r["overlap_fraction"] == 1.0
    assert r["in_window_full"] == 1.0
    # neither window alone contains 190..210
    assert not (100 <= 190 and 210 <= 200)
    assert not (195 <= 190 and 210 <= 300)


# ------------------------------------------------------------------------- the sign convention


def test_positive_means_upstream_on_the_PLUS_strand():
    upstream = _row("A", 4900, "+")["signed_distance_to_next_cds_start"]
    downstream = _row("A", 5500, "+")["signed_distance_to_next_cds_start"]
    assert upstream > 0, "a tile in the plus gene's corpus window is upstream of its start"
    assert downstream < 0, "a tile inside the gene body is downstream"
    assert upstream == pytest.approx(5000 - (4900 + TL // 2))


def test_positive_means_upstream_on_the_MINUS_strand_where_upstream_runs_the_OTHER_WAY():
    """The case most likely to be silently wrong: a minus gene's upstream is HIGHER coordinates."""
    upstream = _row("A", 3100, "-")["signed_distance_to_next_cds_start"]
    downstream = _row("A", 2500, "-")["signed_distance_to_next_cds_start"]
    assert upstream > 0, (
        "a tile in the minus gene's corpus window (3000..3200) sits at HIGHER genome coordinates "
        "than the gene start, and must still score POSITIVE"
    )
    assert downstream < 0
    assert upstream == pytest.approx((3100 + TL // 2) - 3000)


def test_a_tile_in_a_corpus_window_scores_positive_on_BOTH_strands():
    """The convention's whole point: in-window implies positive regardless of strand."""
    assert _row("A", 4900, "+")["signed_distance_to_next_cds_start"] > 0
    assert _row("A", 3100, "-")["signed_distance_to_next_cds_start"] > 0


def test_distance_is_refused_rather_than_sentinelled_when_a_strand_has_no_genes():
    with pytest.raises(FrameFeatureError, match="no CDS rows on strand"):
        frame_features(
            {"A": LocateResult(1, 100, "-")},
            window_coverage([(10, 20, "+", "a")], 1000),
            [CdsRow(500, 600, "+", "p")],
            tile_length=TL,
        )


# -------------------------------------------------------------------------------- inside_cds


def test_a_tile_inside_a_cds_is_flagged_and_one_outside_is_not():
    assert _row("A", 5500, "+")["inside_cds"] == 1.0      # inside gene p (5000..6000)
    assert _row("A", 4900, "+")["inside_cds"] == 0.0      # in the promoter window, not the gene


# ------------------------------------------------------------------------------- accounting


def test_ambiguous_and_not_found_sequences_are_excluded_and_counted_separately():
    t = _table({
        "UNIQ": LocateResult(1, 4900, "+"),
        "AMB": LocateResult(3),
        "GONE": LocateResult(0),
    })
    assert t.seqs == ["UNIQ"]
    assert t.excluded == {"AMB": "ambiguous", "GONE": "not_found"}
    assert t.n_excluded_ambiguous == 1 and t.n_excluded_not_found == 1


def test_emitted_plus_excluded_equals_the_input_count_exactly():
    entries = {"U1": LocateResult(1, 4900, "+"), "U2": LocateResult(1, 100, "+"),
               "A": LocateResult(2), "G": LocateResult(0)}
    t = _table(entries)
    acc = t.accounting()
    assert acc["n_emitted"] + acc["n_excluded"] == acc["n_input_distinct"] == len(entries)


def test_a_table_whose_arithmetic_does_not_reconcile_is_REFUSED():
    """Non-vacuity for the guard above: a reduced denominator must not be constructible."""
    with pytest.raises(FrameFeatureError, match="were supplied"):
        FrameFeatureTable(
            names=FEATURE_NAMES, seqs=["a"], matrix=[[0.0] * len(FEATURE_NAMES)],
            excluded={}, n_input=5,
        )


def test_a_row_of_the_wrong_width_is_refused():
    with pytest.raises(FrameFeatureError, match="values for"):
        FrameFeatureTable(
            names=FEATURE_NAMES, seqs=["a"], matrix=[[0.0, 1.0]], excluded={}, n_input=1,
        )


def test_an_ambiguous_tile_never_gets_an_imputed_frame():
    """Imputing here would invent the very quantity under test."""
    t = _table({"AMB": LocateResult(4)})
    assert t.matrix == [] and t.seqs == []
    assert t.excluded == {"AMB": "ambiguous"}


# ------------------------------------------------------------------------------ table surface


def test_the_feature_names_are_fixed_so_the_bar_can_name_a_primary():
    assert FEATURE_NAMES == (
        "in_window_any", "in_window_full", "overlap_fraction",
        "signed_distance_to_next_cds_start", "inside_cds",
    )


def test_column_and_is_constant_report_what_the_bar_needs():
    t = _table({"A": LocateResult(1, 4900, "+"), "B": LocateResult(1, 100, "+")})
    assert sorted(t.column("in_window_any")) == [0.0, 1.0]
    assert not t.is_constant("in_window_any")
    assert t.is_constant("inside_cds"), "neither tile is inside a gene body here"


def test_an_unknown_column_is_refused():
    t = _table({"A": LocateResult(1, 4900, "+")})
    with pytest.raises(FrameFeatureError, match="no feature"):
        t.column("nope")


# ------------------------------------------------------------------------------- real cohort


@pytest.mark.skipif(
    not (FNA.is_file() and GFF.is_file() and SUBSTRATE.is_file()),
    reason="cached genome, annotation or substrate absent (gitignored)",
)
def test_the_real_cohort_reproduces_the_committed_in_window_counts():
    """28,622 tiles overlap a corpus window and 8,239 are fully covered. Neither is constant."""
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.corpus_window import cds_rows_from_annotation_table
    from dna_decode.glm.tile_locate import locate_sequences, read_fasta_single

    _contig, genome = read_fasta_single(FNA)
    cds = cds_rows_from_annotation_table(parse_gff3(GFF))
    mask = window_coverage(corpus_windows(cds, len(genome)), len(genome))

    rows, _ = gw.load_peak_tiles(REPO / "data/mpra", allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == 150]
    located = locate_sequences([t.seq for t in tiles], genome)

    t = frame_features(located, mask, cds, tile_length=150)

    assert len(t.seqs) == 43_061, "one row per uniquely located sequence"
    assert t.n_excluded_ambiguous == 325
    assert t.n_excluded_not_found == 31
    assert sum(t.column("in_window_any")) == 28_622
    assert sum(t.column("in_window_full")) == 8_239

    for name in FEATURE_NAMES:
        assert not t.is_constant(name), "%s is constant and could absorb nothing" % name

    # the review's refuted worry, kept as a live assertion: the flag is NOT dominated by trivial
    # overlaps -- of the in-window tiles, few overlap by 0.2 or less
    fracs = [f for f in t.column("overlap_fraction") if f > 0]
    trivial = sum(1 for f in fracs if f <= 0.2) / len(fracs)
    assert trivial < 0.15, "measured 8.2%%; got %.4f" % trivial
