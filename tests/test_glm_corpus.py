"""Guards for the natural-corpus extractor.

The expensive failure here is silent: a strand or off-by-one error yields sequence that still looks like
DNA, still has plausible GC, and quietly poisons every comparison downstream. So both strands are pinned
against hand-computed slices, and the real reference is used where available.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.glm.corpus import (
    CorpusStats,
    as_rows,
    extract_upstream_windows,
    gc_fraction,
    load_fasta,
    revcomp,
)
from dna_decode.glm.corpus import _upstream_slice

ROOT = Path(__file__).resolve().parents[1]
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"


def test_revcomp_is_an_involution_and_correct():
    assert revcomp("ACGT") == "ACGT"          # self-complementary
    assert revcomp("AAAA") == "TTTT"
    assert revcomp("ACCGGT") == "ACCGGT"
    s = "ACGTTGCAATTGC"
    assert revcomp(revcomp(s)) == s


# ---------------------------------------------------------------------------------------------------
# the coordinate conversion -- GFF3 is 1-based inclusive, slices are 0-based half-open
# ---------------------------------------------------------------------------------------------------
def test_plus_strand_upstream_slice_ends_exactly_at_the_cds_start():
    """A CDS starting at 1-based 101 occupies 0-based index 100, so a 10-bp window is [90, 100)."""
    assert _upstream_slice(101, 200, "+", 10) == (90, 100)


def test_minus_strand_upstream_slice_begins_exactly_after_the_cds_end():
    """A minus-strand CDS ending at 1-based 200 occupies 0-based up to 199, so upstream is [200, 210)."""
    assert _upstream_slice(101, 200, "-", 10) == (200, 210)


def test_an_unknown_strand_is_refused_rather_than_guessed():
    with pytest.raises(ValueError, match="strand must be"):
        _upstream_slice(101, 200, ".", 10)


# ---------------------------------------------------------------------------------------------------
# extraction on a synthetic genome where the answer is readable by eye
# ---------------------------------------------------------------------------------------------------
#        idx 0....5....10...15...20...25...30...35
GEN = {"c1": "AAAAAAAAAA" "CCCCCCCCCC" "GGGGGGGGGG" "TTTTTTTTTT"}


def _row(start, end, strand, gid="g"):
    return {"seqid": "c1", "type": "CDS", "start": start, "end": end, "strand": strand,
            "gene_id": gid, "gene_symbol": "", "locus_tag": ""}


def test_plus_strand_window_is_the_bases_immediately_before_the_start():
    # CDS starts 1-based 21 (0-based 20, the G block). 10 bp upstream = the C block.
    wins, _ = extract_upstream_windows(GEN, [_row(21, 30, "+")], length=10)
    assert len(wins) == 1
    assert wins[0].seq == "CCCCCCCCCC"
    assert (wins[0].start, wins[0].end) == (10, 20)


def test_minus_strand_window_is_AFTER_the_end_AND_reverse_complemented():
    """The load-bearing case. Upstream of a minus-strand CDS ending at 1-based 20 is [20,30) = the G
    block, and the extracted promoter must be its REVERSE COMPLEMENT (CCCCCCCCCC), not the G block."""
    wins, _ = extract_upstream_windows(GEN, [_row(11, 20, "-")], length=10)
    assert len(wins) == 1
    assert wins[0].seq == "CCCCCCCCCC", "minus-strand window must be reverse-complemented"
    assert (wins[0].start, wins[0].end) == (20, 30), "forward coords are reported un-flipped"


def test_the_revcomp_step_is_NON_VACUOUS_the_raw_slice_differs():
    """If the raw slice equalled its own reverse complement this test would prove nothing."""
    raw = GEN["c1"][20:30]
    assert raw == "GGGGGGGGGG"
    assert revcomp(raw) != raw, "pick a non-palindromic block or this guard is vacuous"


def test_a_window_running_off_the_contig_start_is_skipped_not_truncated():
    """A short window would change the length distribution the discriminator keys on."""
    wins, stats = extract_upstream_windows(GEN, [_row(5, 30, "+")], length=10)
    assert wins == [] and stats.n_skipped_edge == 1


def test_a_window_running_off_the_contig_end_is_skipped():
    wins, stats = extract_upstream_windows(GEN, [_row(1, 38, "-")], length=10)
    assert wins == [] and stats.n_skipped_edge == 1


def test_ambiguous_bases_are_skipped_and_counted():
    g = {"c1": "ACGTNNNNNN" + "ACGTACGTAC"}
    wins, stats = extract_upstream_windows(g, [_row(11, 20, "+")], length=10)
    assert wins == [] and stats.n_skipped_ambiguous == 1
    # and with the check off, it is admitted
    wins2, _ = extract_upstream_windows(g, [_row(11, 20, "+")], length=10, require_unambiguous=False)
    assert len(wins2) == 1


def test_overlap_with_another_cds_is_measured_and_optionally_excluded():
    rows = [_row(21, 30, "+", "a"), _row(11, 20, "+", "b")]   # b sits inside a's upstream window
    wins, stats = extract_upstream_windows(GEN, rows, length=10)
    assert stats.n_overlapping_cds >= 1
    assert any(w.overlaps_cds for w in wins)
    kept, _ = extract_upstream_windows(GEN, rows, length=10, exclude_overlapping=True)
    assert len(kept) < len(wins)


def test_a_non_positive_length_is_refused():
    with pytest.raises(ValueError, match="length must be positive"):
        extract_upstream_windows(GEN, [_row(21, 30, "+")], length=0)


def test_non_cds_features_are_ignored():
    rows = [dict(_row(21, 30, "+"), type="gene")]
    wins, _ = extract_upstream_windows(GEN, rows, length=10)
    assert wins == []


# ---------------------------------------------------------------------------------------------------
# the DataFrame-vs-list shape trap
# ---------------------------------------------------------------------------------------------------
def test_as_rows_accepts_a_dataframe_because_list_of_a_dataframe_gives_COLUMN_NAMES():
    """`parse_gff3` returns a pandas DataFrame and `list(df)` yields column names, which surfaced only as
    `'str' object has no attribute 'get'` deep inside extraction. Both shapes must work."""
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame([_row(21, 30, "+")])
    assert list(df) != as_rows(df), "list(df) gives columns; as_rows must give records"
    assert as_rows(df)[0]["start"] == 21
    assert as_rows([_row(21, 30, "+")])[0]["start"] == 21
    wins, _ = extract_upstream_windows(GEN, df, length=10)
    assert len(wins) == 1 and wins[0].seq == "CCCCCCCCCC"


def test_gc_fraction():
    assert gc_fraction("GCGC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert gc_fraction("ACGT") == 0.5
    assert gc_fraction("") == 0.0, "empty must not raise"


# ---------------------------------------------------------------------------------------------------
# the real reference
# ---------------------------------------------------------------------------------------------------
def test_extraction_on_the_REAL_mg1655_reference():
    """Real-surface test (R3): thousands of genuine promoters, both strands, exact lengths."""
    if not (REFSEQ / "genome.fna").exists():
        pytest.skip("MG1655 reference not on this host")
    from dna_decode.data.annotations import parse_gff3
    genome = load_fasta(REFSEQ / "genome.fna")
    assert len(genome) == 1 and len(next(iter(genome.values()))) == 4641652
    rows = parse_gff3(REFSEQ / "annotations.gff3")
    wins, stats = extract_upstream_windows(genome, rows, length=150)
    assert stats.n_windows > 3000, f"expected thousands of CDS windows, got {stats.n_windows}"
    assert all(len(w.seq) == 150 for w in wins), "every window must be exactly the requested length"
    assert all(not (set(w.seq) - set("ACGT")) for w in wins)
    # both strands are genuinely represented -- a strand bug often shows up as one strand missing
    assert stats.strand_counts.get("+", 0) > 1000
    assert stats.strand_counts.get("-", 0) > 1000
    # real promoters are AT-richer than the genome average (~0.508 GC for E. coli)
    gc = sum(gc_fraction(w.seq) for w in wins) / len(wins)
    assert 0.35 < gc < 0.50, f"promoter GC {gc:.3f} is implausible for E. coli upstream regions"


# ---------------------------------------------------------------------------------------------------
# the genomic reference set + conditional prompts
# ---------------------------------------------------------------------------------------------------
def test_random_genomic_windows_are_exact_length_both_strands_and_acgt():
    from dna_decode.glm.corpus import random_genomic_windows
    g = {"c1": "ACGT" * 500}
    w, s = random_genomic_windows(g, length=60, n=50, seed=3)
    assert s.n_windows == 50 and all(len(x.seq) == 60 for x in w)
    assert all(not (set(x.seq) - set("ACGT")) for x in w)
    assert set(s.strand_counts) == {"+", "-"}, "both orientations must be sampled"


def test_random_genomic_windows_are_deterministic_under_a_seed():
    from dna_decode.glm.corpus import random_genomic_windows
    g = {"c1": "ACGTTGCA" * 300}
    a, _ = random_genomic_windows(g, length=40, n=20, seed=11)
    b, _ = random_genomic_windows(g, length=40, n=20, seed=11)
    assert [x.seq for x in a] == [x.seq for x in b]


def test_random_genomic_windows_skip_ambiguous_rather_than_emit_N():
    from dna_decode.glm.corpus import random_genomic_windows
    g = {"c1": "N" * 1000}
    w, s = random_genomic_windows(g, length=50, n=10, seed=1)
    assert w == [] and s.n_skipped_ambiguous > 0


def test_the_genomic_set_matches_the_REAL_genome_gc_and_differs_from_promoters():
    """The correctness check that matters, and it is the reason this set exists: the genomic set must
    reproduce E. coli's ~0.508 GC while promoters are AT-richer. If the two sets had the same
    composition, choosing between them would not matter -- and it does."""
    if not (REFSEQ / "genome.fna").exists():
        pytest.skip("MG1655 reference not on this host")
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.glm.corpus import random_genomic_windows
    genome = load_fasta(REFSEQ / "genome.fna")
    gw, _ = random_genomic_windows(genome, length=150, n=400, seed=0)
    gc_gen = sum(gc_fraction(x.seq) for x in gw) / len(gw)
    assert 0.49 < gc_gen < 0.53, f"genomic GC {gc_gen:.4f} should be ~0.508"
    pw, _ = extract_upstream_windows(genome, parse_gff3(REFSEQ / "annotations.gff3"), length=150)
    gc_prom = sum(gc_fraction(x.seq) for x in pw) / len(pw)
    assert gc_prom < gc_gen - 0.03, (
        f"promoters ({gc_prom:.4f}) must be measurably AT-richer than the genome ({gc_gen:.4f}); "
        "if not, the two reference sets are interchangeable and the distinction is pointless"
    )


def test_window_context_is_strand_aware_and_is_the_bases_PRECEDING_the_window():
    """A prompt from the wrong side, or un-complemented, would make the model continue the wrong strand."""
    from dna_decode.glm.corpus import window_context
    # plus strand: window [10,20) -> context [5,10) is the tail of the A block
    w_plus = extract_upstream_windows(GEN, [_row(21, 30, "+")], length=10)[0][0]
    assert window_context(GEN, w_plus, 5) == "AAAAA"
    # minus strand: window forward-coords [20,30); its 5' neighbour is [30,35) revcomp'd
    w_minus = extract_upstream_windows(GEN, [_row(11, 20, "-")], length=10)[0][0]
    assert window_context(GEN, w_minus, 5) == revcomp(GEN["c1"][30:35]) == "AAAAA"


def test_window_context_returns_None_rather_than_a_short_prompt_at_a_contig_edge():
    """A truncated prompt would vary prompt length across the set, which the model conditions on."""
    from dna_decode.glm.corpus import window_context
    w = extract_upstream_windows(GEN, [_row(21, 30, "+")], length=10)[0][0]
    assert window_context(GEN, w, 50) is None      # would run off the start
    assert window_context(GEN, w, 0) is None       # degenerate request


def test_window_context_refuses_ambiguous_context():
    from dna_decode.glm.corpus import window_context
    g = {"c1": "NNNNN" + "AAAAA" + "CCCCCCCCCC" + "GGGGGGGGGG"}
    w = extract_upstream_windows(g, [_row(21, 30, "+")], length=10)[0][0]
    assert window_context(g, w, 10) is None, "an N-containing prompt must be skipped, not passed in"


def test_corpus_stats_round_trip():
    s = CorpusStats(n_cds_total=5, n_windows=3)
    assert s.as_dict()["n_cds_total"] == 5 and s.as_dict()["n_windows"] == 3
