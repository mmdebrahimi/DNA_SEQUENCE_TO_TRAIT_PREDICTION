"""Tests for the PromoterAtlas corpus-window replica.

**The two quirk tests are the load-bearing ones, and they are written so a tidy-up FAILS.** Each
constructs a case whose emitted window set differs between the upstream rule and the "corrected"
rule, implements the corrected rule inline so the difference is demonstrated rather than asserted,
and then pins the upstream answer. A test that merely asserted a window count could stay green
after someone moved the `last_end` assignment above the `continue`.

The real-data test is the faithfulness check: it asserts the module reproduces every figure
`wiki/frame_plan_window_ordering_2026-11-01.json` committed from
`scripts/frame_window_ordering_probe.py`, under BOTH orderings. Asserting both is what makes the
`order` parameter non-decorative -- a single-ordering assertion would pass whichever one the
default happened to be.
"""

from __future__ import annotations

import pathlib

import pandas as pd
import pytest

from dna_decode.glm.corpus_window import (
    MIN_INTERGENIC_LENGTH,
    REGULATORY_REGION_LENGTH,
    CdsRow,
    CorpusWindowError,
    Window,
    cds_rows_from_annotation_table,
    corpus_windows,
    overlap_fraction,
    window_coverage,
)

REPO = pathlib.Path(__file__).resolve().parents[1]
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"


# --------------------------------------------------------------------- placement and the gates


def test_forward_window_sits_immediately_upstream_of_the_gene_start():
    w = corpus_windows([CdsRow(1000, 2000, "+", "a")], 10_000, order="file")
    assert w == [Window(800, 1000, "+", "a")]


def test_reverse_window_sits_downstream_in_genome_coordinates():
    """Downstream in genome coordinates IS upstream relative to a minus-strand gene."""
    w = corpus_windows([CdsRow(2000, 3000, "-", "r")], 10_000, order="file")
    assert w == [Window(3000, 3200, "-", "r")]


def test_a_gene_too_close_to_the_replicon_start_is_skipped():
    """There is no room for a 200 bp window upstream of a gene starting at 150."""
    assert corpus_windows([CdsRow(150, 900, "+", "a")], 10_000, order="file") == []


def test_a_reverse_gene_too_close_to_the_replicon_end_is_skipped():
    assert corpus_windows([CdsRow(9000, 9900, "-", "r")], 10_000, order="file") == []


def test_the_min_intergenic_gate_rejects_a_tight_neighbour():
    rows = [CdsRow(1000, 2000, "+", "a"), CdsRow(2003, 2500, "+", "b")]
    w = corpus_windows(rows, 10_000, order="file")
    assert w == [Window(800, 1000, "+", "a")], "b sits 3 bp after a, inside the 5 bp gate"


def test_the_gate_boundary_is_inclusive_at_exactly_min_intergenic_length():
    rows = [CdsRow(1000, 2000, "+", "a"), CdsRow(2005, 2500, "+", "b")]
    w = corpus_windows(rows, 10_000, order="file")
    assert [x.locus_tag for x in w] == ["a", "b"], "the comparison is >=, so 5 bp passes"


# ------------------------------------------------------------- the two deliberately-kept quirks


def _forward_windows_CORRECTED(rows, genome_length, rl=REGULATORY_REGION_LENGTH,
                               ig=MIN_INTERGENIC_LENGTH):
    """The forward rule with quirk 1 TIDIED: a skipped gene updates `last_end`.

    Here only so the test can demonstrate the two rules disagree, rather than assert it.
    """
    out, last_end = [], 0
    for r in [x for x in rows if x.strand == "+"]:
        if r.start0 < rl:
            last_end = r.end0          # the tidy-up
            continue
        if r.start0 - last_end >= ig:
            out.append(Window(r.start0 - rl, r.start0, "+", r.locus_tag))
        last_end = r.end0
    return out


def test_a_skipped_gene_does_NOT_update_last_end_and_tidying_it_changes_the_window_set():
    """QUIRK 1. The `continue` precedes the `last_end` assignment in the upstream source.

    Gene `early` starts at 100 so it is skipped, but it ENDS at 5000 -- past the next gene's
    start. Under the upstream rule `last_end` is still 0 when `b` is considered, so `b` emits.
    Under the tidied rule `last_end` is 5000 and `b` is gated out.
    """
    rows = [CdsRow(100, 5000, "+", "early"), CdsRow(4000, 4500, "+", "b")]

    upstream = corpus_windows(rows, 10_000, order="file")
    tidied = _forward_windows_CORRECTED(rows, 10_000)

    assert upstream != tidied, "the case no longer discriminates; it cannot pin the quirk"
    assert tidied == [], "the tidied rule gates b out against last_end=5000"
    assert upstream == [Window(3800, 4000, "+", "b")], (
        "upstream threads last_end=0 past a skipped gene; if this now fails, someone moved the "
        "last_end assignment above the continue and changed which corpus the replica describes"
    )


def _reverse_windows_CORRECTED(rows, genome_length, rl=REGULATORY_REGION_LENGTH,
                               ig=MIN_INTERGENIC_LENGTH):
    """The reverse rule with quirk 2 TIDIED: `last_end_rev` tracks the previous gene's END."""
    out, last_end_rev = [], genome_length
    for r in reversed([x for x in rows if x.strand == "-"]):
        if genome_length - r.end0 < rl:
            continue
        if last_end_rev - r.end0 >= ig:
            out.append(Window(r.end0, r.end0 + rl, "-", r.locus_tag))
        last_end_rev = r.end0          # the tidy-up: END, not START
    return out


def test_the_reverse_branch_compares_against_the_previous_genes_START():
    """QUIRK 2. `last_end_reverse = start` closes the upstream loop body.

    `r1` spans 6000..7000, so its START (6000) and END (7000) straddle `r2`'s end (6500). Under
    the upstream rule the gap is 6000 - 6500 = -500 and `r2` is rejected; under the tidied rule it
    is 7000 - 6500 = 500 and `r2` emits. Because a start is always below its own end, the upstream
    rule is the STRICTER one -- which is why a discriminating case needs a long gene.
    """
    rows = [CdsRow(5800, 6500, "-", "r2"), CdsRow(6000, 7000, "-", "r1")]

    upstream = corpus_windows(rows, 10_000, order="file")
    tidied = _reverse_windows_CORRECTED(rows, 10_000)

    assert upstream != tidied, "the case no longer discriminates; it cannot pin the quirk"
    assert [w.locus_tag for w in tidied] == ["r1", "r2"]
    assert upstream == [Window(7000, 7200, "-", "r1")], (
        "upstream compares r2's end against r1's START and rejects it; if this now fails, the "
        "reverse branch was tidied to track ends and the replica no longer matches the corpus"
    )


def test_the_forward_and_reverse_branches_disagree_about_what_last_end_holds():
    """The asymmetry is itself the finding: one branch tracks ends, the other tracks starts."""
    src = pathlib.Path(
        REPO / "dna_decode/glm/corpus_window.py"
    ).read_text(encoding="utf-8")
    assert "last_end = r.end0" in src
    assert "last_end_rev = r.start0" in src


# ------------------------------------------------------------------- ordering and the adapter


def test_an_unknown_order_is_refused_rather_than_defaulted():
    with pytest.raises(CorpusWindowError, match="order="):
        corpus_windows([CdsRow(1000, 2000, "+", "a")], 10_000, order="whatever")


def test_order_actually_changes_the_result_so_the_parameter_is_not_decorative():
    """File order is not coordinate order, and the threading makes that observable.

    `a` spans 5000..6000 and `b` spans 1000..2000, given in that (non-ascending) file order. Taken
    in file order, `a` emits and advances `last_end` to 6000, which gates `b` out entirely. Sorted,
    `b` emits first and advances `last_end` only to 2000, so BOTH emit. Nothing here depends on
    either preserved quirk -- this is the threading alone, which is why it is the right case for
    this assertion. (An earlier version of this test used a case the quirk happened to make
    order-invariant, and it passed vacuously until it was asserted the other way round.)
    """
    rows = [CdsRow(5000, 6000, "+", "a"), CdsRow(1000, 2000, "+", "b")]

    in_file = corpus_windows(rows, 10_000, order="file")
    in_sorted = corpus_windows(rows, 10_000, order="sorted")

    assert [w.locus_tag for w in in_file] == ["a"]
    assert [w.locus_tag for w in in_sorted] == ["b", "a"]
    assert in_file != in_sorted


def test_the_gff3_one_based_inclusive_conversion_happens_exactly_once():
    """A GFF3 row at 1-based inclusive 1001..2000 is 0-based half-open 1000..2000."""
    tab = pd.DataFrame(
        [
            {"type": "CDS", "start": 1001, "end": 2000, "strand": "+", "locus_tag": "a"},
            {"type": "CDS", "start": 3001, "end": 4000, "strand": "-", "locus_tag": "b"},
        ]
    )
    rows = cds_rows_from_annotation_table(tab)
    assert rows == [CdsRow(1000, 2000, "+", "a"), CdsRow(3000, 4000, "-", "b")]


def test_the_adapter_drops_rows_without_a_locus_tag_like_the_extractor_does():
    tab = pd.DataFrame(
        [
            {"type": "CDS", "start": 1001, "end": 2000, "strand": "+", "locus_tag": ""},
            {"type": "CDS", "start": 3001, "end": 4000, "strand": "+", "locus_tag": "b"},
            {"type": "gene", "start": 1, "end": 500, "strand": "+", "locus_tag": "g"},
        ]
    )
    assert [r.locus_tag for r in cds_rows_from_annotation_table(tab)] == ["b"]


# ----------------------------------------------------------------------- mask and overlap maths


def test_window_coverage_counts_an_overlapped_base_once():
    cov = window_coverage([Window(10, 20, "+", "a"), Window(15, 25, "+", "b")], 100)
    assert sum(cov) == 15
    assert set(cov[10:25]) == {1}


def test_overlap_fraction_is_strictly_between_zero_and_one_on_a_straddling_tile():
    cov = window_coverage([Window(100, 200, "+", "a")], 1000)
    assert overlap_fraction(cov, 150, 100) == pytest.approx(0.5)
    assert overlap_fraction(cov, 100, 100) == 1.0
    assert overlap_fraction(cov, 500, 100) == 0.0


def test_overlap_fraction_refuses_a_tile_running_off_the_mask():
    """A truncated denominator reads as a LOWER overlap rather than as an error."""
    cov = window_coverage([Window(100, 200, "+", "a")], 1000)
    with pytest.raises(CorpusWindowError, match="outside"):
        overlap_fraction(cov, 950, 100)


# ------------------------------------------------------------------------------- real annotation


@pytest.mark.skipif(
    not (GFF.is_file() and FNA.is_file()),
    reason="cached MG1655 annotation absent (gitignored substrate)",
)
def test_the_real_annotation_reproduces_BOTH_committed_orderings():
    """The faithfulness check against `wiki/frame_plan_window_ordering_2026-11-01.json`."""
    from dna_decode.data.annotations import parse_gff3

    genome_length = sum(
        len(ln.strip())
        for ln in FNA.read_text(encoding="utf-8").splitlines()
        if not ln.startswith(">")
    )
    assert genome_length == 4_641_652

    rows = cds_rows_from_annotation_table(parse_gff3(GFF))
    assert len(rows) == 4340, "CDS rows carrying a locus_tag"

    w_sorted = corpus_windows(rows, genome_length, order="sorted")
    w_file = corpus_windows(rows, genome_length, order="file")

    assert len(w_sorted) == 3496
    assert sum(window_coverage(w_sorted, genome_length)) == 625_255
    assert len(w_file) == 3490
    assert sum(window_coverage(w_file, genome_length)) == 624_812

    # The probe recorded the symmetric difference on (start, end, strand); compare on the same
    # projection so the number stays comparable with the committed artifact.
    s3 = {(w.start0, w.end0, w.strand) for w in w_sorted}
    f3 = {(w.start0, w.end0, w.strand) for w in w_file}
    assert len(s3 ^ f3) == 26


@pytest.mark.skipif(
    not (GFF.is_file() and FNA.is_file()),
    reason="cached MG1655 annotation absent (gitignored substrate)",
)
def test_corpus_coverage_is_a_minority_of_the_genome_so_membership_has_variance():
    """13.47% coverage is what makes the absorption test possible; a near-constant could not."""
    from dna_decode.data.annotations import parse_gff3

    genome_length = 4_641_652
    rows = cds_rows_from_annotation_table(parse_gff3(GFF))
    frac = sum(window_coverage(corpus_windows(rows, genome_length), genome_length)) / genome_length
    assert 0.13 < frac < 0.14
