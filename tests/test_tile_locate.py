"""Tests for locating peak tiles by sequence.

Three things here are load-bearing. The **structural** test uses `ast` to assert the module never
reads a coordinate column -- a substring scan would match the module's own docstring, which
discusses those columns at length, so a text search here is guaranteed to be either vacuous or
wrong. The **reverse-strand non-vacuity** test plants a sequence that exists ONLY as a reverse
complement, because every other test would pass with the reverse sweep deleted. And the
**real-data** test asserts the extracted sweep reproduces the committed bucket counts exactly,
which is what makes this an extraction rather than a rewrite.
"""

from __future__ import annotations

import ast
import collections
import gzip
import pathlib

import pytest

from dna_decode.glm.tile_locate import (
    CACHE_SCHEMA,
    LOCATE_FLOOR,
    LocateReport,
    LocateResult,
    TileLocateError,
    load_locate_cache,
    locate_report,
    locate_sequences,
    read_fasta_single,
    reverse_complement,
    save_locate_cache,
    unique_positions,
)

REPO = pathlib.Path(__file__).resolve().parents[1]
MODULE = REPO / "dna_decode/glm/tile_locate.py"
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
SUBSTRATE = REPO / "data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz"

TILE_LEN = 150


def _synthetic_genome(n: int = 2000, seed: int = 7) -> str:
    import random

    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


# ------------------------------------------------------------------------- synthetic placement


def test_a_forward_sequence_is_located_at_its_own_position():
    g = _synthetic_genome()
    want = g[300:300 + 20]
    got = locate_sequences([want], g)[want]
    assert got.is_unique and got.strand == "+" and got.pos0 == 300


def test_a_reverse_complement_sequence_is_normalised_back_to_forward_coordinates():
    """The returned position must index the FORWARD genome, not the reversed string."""
    g = _synthetic_genome()
    planted = g[500:500 + 20]
    query = reverse_complement(planted)
    if query in g:                     # a self-complementary accident would void the test
        pytest.skip("the planted sequence happens to occur forward too")
    got = locate_sequences([query], g)[query]
    assert got.is_unique and got.strand == "-"
    assert got.pos0 == 500, "a reverse hit must map back to the forward coordinate"
    assert g[got.pos0:got.pos0 + 20] == planted


def test_the_reverse_sweep_is_genuinely_exercised():
    """Non-vacuity: this sequence is present ONLY as a reverse complement."""
    g = _synthetic_genome()
    query = reverse_complement(g[800:800 + 25])
    if query in g:
        pytest.skip("the planted sequence happens to occur forward too")
    assert locate_sequences([query], g)[query].is_unique
    # and with the reverse strand not swept, it would be not-found:
    assert query not in g


def test_a_duplicated_sequence_is_excluded_and_counted_never_resolved_to_the_first_hit():
    base = _synthetic_genome(400)
    dup = base[50:50 + 20]
    g = base + dup                     # now `dup` occurs at 50 and at 400
    got = locate_sequences([dup], g)[dup]
    assert not got.is_unique
    assert got.n_occurrences >= 2
    assert got.pos0 is None and got.strand is None, (
        "an ambiguous sequence must carry no position, so using one raises instead of quietly "
        "computing a frame from an arbitrary hit"
    )


def test_an_absent_sequence_is_counted_as_not_found():
    g = _synthetic_genome()
    missing = "A" * 30
    assert missing not in g
    got = locate_sequences([missing], g)[missing]
    assert got.n_occurrences == 0 and got.pos0 is None


def test_a_ragged_cohort_is_refused_rather_than_swept_at_one_width():
    g = _synthetic_genome()
    with pytest.raises(TileLocateError, match="ragged"):
        locate_sequences([g[10:30], g[10:35]], g)


def test_unique_positions_returns_only_the_unambiguous_rows():
    g = _synthetic_genome(400)
    dup = g[50:50 + 20]
    g2 = g + dup
    uniq = g[100:100 + 20]
    located = locate_sequences([dup, uniq], g2)
    assert set(unique_positions(located)) == {uniq}


# -------------------------------------------------------------------------------- accounting


def test_the_report_buckets_sum_to_the_distinct_count():
    g = _synthetic_genome(400)
    dup = g[50:50 + 20]
    g2 = g + dup
    seqs = [g[100:100 + 20], g[200:200 + 20], dup, "T" * 20]
    located = locate_sequences(seqs, g2)
    rep = locate_report(located, tile_length=20, genome_length=len(g2))
    assert rep.n_distinct == 4
    assert rep.n_unique + rep.n_ambiguous + rep.n_not_found == 4
    assert rep.n_located == rep.n_unique + rep.n_ambiguous


def test_a_report_whose_buckets_do_not_reconcile_is_refused():
    """A reduced denominator reads as a BETTER locate rate rather than as an error."""
    with pytest.raises(TileLocateError, match="buckets sum"):
        LocateReport(
            n_distinct=10, n_unique=5, n_ambiguous=1, n_not_found=1,
            tile_length=150, genome_length=1000,
        )


def test_a_rate_below_the_floor_does_not_raise_so_the_graders_branch_stays_reachable():
    """A module that crashed first would make INDETERMINATE_LOCATE_RATE dead code."""
    rep = LocateReport(
        n_distinct=10, n_unique=5, n_ambiguous=3, n_not_found=2,
        tile_length=150, genome_length=1000,
    )
    assert rep.unique_rate == 0.5
    assert not rep.meets_floor and not rep.frame_claim_permitted
    assert rep.verdict == "LOCATE_BELOW_FLOOR"
    with pytest.raises(TileLocateError, match=r"< floor"):
        rep.require_floor()            # available, but only when a caller asks


def test_the_floor_has_headroom_over_the_measured_rate():
    assert LOCATE_FLOOR == 0.95, "0.95 sits below the measured 0.9918, not fitted down to it"


# ------------------------------------------------------- the structural no-coordinates guarantee


def test_the_module_never_reads_a_coordinate_COLUMN_asserted_with_ast_not_text():
    """A substring scan would match this module's own docstring, which discusses those columns."""
    tree = ast.parse(MODULE.read_text(encoding="utf-8"))

    banned_strings = {"peak_start", "peak_end", "tile_start", "tile_end"}
    banned_attrs = {"start", "end", "peak_start", "peak_end"}

    bad_consts, bad_attrs = [], []
    for node in ast.walk(tree):
        # a docstring is an ast.Expr wrapping a Constant; skip those, keep real literals
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if node.value in banned_strings:
                bad_consts.append(node.value)
        if isinstance(node, ast.Attribute) and node.attr in banned_attrs:
            bad_attrs.append(node.attr)

    assert not bad_consts, "module names coordinate columns in executable positions: %s" % bad_consts
    assert not bad_attrs, "module reads coordinate attributes: %s" % bad_attrs


def test_the_structural_guard_is_non_vacuous():
    """Proof the check above can fail: the same walk fires on a snippet that does read them."""
    tree = ast.parse("row['peak_start']\nx = tile.start\n")
    consts = [
        n.value for n in ast.walk(tree)
        if isinstance(n, ast.Constant) and n.value == "peak_start"
    ]
    attrs = [n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)]
    assert consts == ["peak_start"] and "start" in attrs


# -------------------------------------------------------------------------------------- cache


def test_the_cache_round_trips_including_the_ambiguous_and_not_found_rows(tmp_path):
    """Their COUNTS are part of the accounting, so a cache of only usable rows cannot reconcile."""
    located = {
        "AAAA": LocateResult(1, 10, "+"),
        "CCCC": LocateResult(3),
        "GGGG": LocateResult(0),
    }
    rep = locate_report(located, tile_length=4, genome_length=100)
    p = save_locate_cache(tmp_path / "loc.json", located, report=rep)
    back = load_locate_cache(p)
    assert back == located
    assert locate_report(back, tile_length=4, genome_length=100).as_dict() == rep.as_dict()


def test_an_unexpected_cache_schema_is_refused(tmp_path):
    p = tmp_path / "loc.json"
    p.write_text('{"schema": "something-else", "located": {}}', encoding="utf-8")
    with pytest.raises(TileLocateError, match="schema"):
        load_locate_cache(p)
    assert CACHE_SCHEMA == "tile-locate-cache-v1"


# --------------------------------------------------------------------------------- real genome


def test_read_fasta_single_refuses_a_multi_record_file(tmp_path):
    p = tmp_path / "two.fna"
    p.write_text(">a\nACGT\n>b\nTTTT\n", encoding="utf-8")
    with pytest.raises(TileLocateError, match="FASTA records"):
        read_fasta_single(p)


@pytest.mark.skipif(
    not (FNA.is_file() and SUBSTRATE.is_file()),
    reason="cached genome or substrate absent (gitignored)",
)
def test_the_real_cohort_reproduces_the_committed_bucket_counts():
    """The faithfulness check: 43,061 unique of 43,417 distinct, 325 ambiguous, 31 not found."""
    from dna_decode.glm import genomewide as gw

    contig, genome = read_fasta_single(FNA)
    assert contig == "NC_000913.3"
    assert len(genome) == 4_641_652

    rows, _ = gw.load_peak_tiles(REPO / "data/mpra", allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == TILE_LEN]
    assert len(tiles) == 44_106

    located = locate_sequences([t.seq for t in tiles], genome)
    rep = locate_report(located, tile_length=TILE_LEN, genome_length=len(genome))

    assert rep.n_distinct == 43_417
    assert rep.n_unique == 43_061
    assert rep.n_ambiguous == 325
    assert rep.n_not_found == 31
    assert rep.n_located == 43_386
    assert rep.unique_rate == pytest.approx(0.9918, abs=5e-5)
    assert rep.meets_floor


@pytest.mark.skipif(
    not (FNA.is_file() and SUBSTRATE.is_file()),
    reason="cached genome or substrate absent (gitignored)",
)
def test_the_stated_coordinates_disagree_with_the_genome_in_OFFSET_FAMILIES():
    """Why the module locates by sequence. Asserted so the 0-vs-1-based convention cannot matter.

    An off-by-one in the stated column shifts every offset uniformly, so it can change WHICH
    offsets are the modes but not that there are several, nor that near-zero agreement is rare.
    Both asserted facts survive either convention; the exact committed offsets (+1978/+1977/+1976/
    +777/+776/0) live in the probe's artifact rather than being re-pinned here under a convention
    this test deliberately does not assume.
    """
    contig, genome = read_fasta_single(FNA)
    uniq = unique_positions(locate_sequences(
        _substrate_tile_sequences(), genome
    ))

    offsets = []
    for seq, stated in _substrate_stated_starts().items():
        if seq in uniq:
            offsets.append(uniq[seq][0] - stated)
    assert len(offsets) > 40_000, "too few comparable rows to say anything"

    counts = collections.Counter(offsets)
    top_share = counts.most_common(1)[0][1] / len(offsets)
    assert top_share < 0.5, (
        "a single offset explaining most rows would be a correctable constant shift; measured "
        "top share %.4f" % top_share
    )
    families = [o for o, c in counts.items() if c >= 50]
    assert len(families) >= 3, "offsets should form FAMILIES (indel drift), got %r" % families

    agree = sum(counts.get(d, 0) for d in (-1, 0, 1)) / len(offsets)
    assert agree < 0.05, (
        "the stated column agrees at %.4f under either convention -- a coincidence rate, which is "
        "why reading it would produce a self-consistent and entirely wrong frame" % agree
    )


def _substrate_rows():
    with gzip.open(SUBSTRATE, "rt", errors="replace") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) == len(hdr):
                yield q, ix


def _substrate_tile_sequences() -> list[str]:
    out = []
    for q, ix in _substrate_rows():
        seq = q[ix["variant"]].upper()
        if len(seq) == TILE_LEN and not (set(seq) - set("ACGT")):
            if q[ix["category"]] == "tile":
                out.append(seq)
    return out


def _substrate_stated_starts() -> dict[str, int]:
    out: dict[str, int] = {}
    for q, ix in _substrate_rows():
        seq = q[ix["variant"]].upper()
        if len(seq) != TILE_LEN or q[ix["category"]] != "tile":
            continue
        try:
            out[seq] = int(q[ix["start"]])
        except ValueError:
            continue
    return out
