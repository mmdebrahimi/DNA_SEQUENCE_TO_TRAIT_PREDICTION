"""Guards for the genome-wide MPRA substrate.

Four traps are encoded in the module and pinned here. Two of them CONTRADICT what the sibling GSE108535
loader taught, which is the reason this file exists rather than reusing those tests: the delimiter is mixed
within one GEO series, and the expression column is already normalised.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from dna_decode.glm.genomewide import (
    MIN_COLUMNS,
    Fragment,
    GenomeWideDataError,
    condition_effect,
    gc_features,
    kmer_features,
    load_fragments,
    noise_ceiling,
    position_blocked_split,
    sniff_delimiter,
)

REAL = Path("D:/dna_decode_cache/mpra/GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz")
HAVE_REAL = REAL.exists()

COLS = ("fragment RNA_exp_1 RNA_exp_2 RNA_exp_ave DNA_sum_1 DNA_sum_2 DNA_ave "
        "num_mapped_barcodes num_integrated_barcodes start end strand variation").split()


def _write(path: Path, rows, delim=" "):
    with gzip.open(path, "wt") as fh:
        fh.write(delim.join(COLS) + "\n")
        for r in rows:
            fh.write(delim.join(str(x) for x in r) + "\n")


def _row(seq="ACGT" * 30, e1=1.0, e2=1.1, ave=1.05, start=1000, end=1120, strand="+"):
    return [seq, e1, e2, ave, 1, 1, 1, 1, 1, start, end, strand, 0.1]


def _frags(n=400, seed=0):
    import random
    r = random.Random(seed)
    out = []
    for i in range(n):
        s = i * 300
        true = r.random()
        out.append(Fragment(seq="".join(r.choice("ACGT") for _ in range(244)),
                            expression=true, rep1=true + r.gauss(0, .3), rep2=true + r.gauss(0, .3),
                            start=s, end=s + 244, strand="+"))
    return out


# ---------------------------------------------------------------------------------------------------
# TRAP 1 -- the delimiter is MIXED within one GEO series
# ---------------------------------------------------------------------------------------------------
def test_sniff_detects_tab_and_whitespace_separately():
    """Measured: the two frag files are SPACE-delimited and the peak-tile file is TAB-delimited. A global
    assumption parses one of them to zero rows."""
    assert sniff_delimiter("a\tb\tc") == "\t"
    assert sniff_delimiter("a b c") is None          # None => split on any whitespace


def test_load_refuses_when_the_delimiter_guess_collapses_the_header(tmp_path):
    """A wrong delimiter yields one giant column. That must RAISE, naming the file, rather than parse
    zero rows -- a zero-row parse that is then averaged reads as 'no signal', not 'no data'."""
    p = tmp_path / "GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("fragment;RNA_exp_1;RNA_exp_2\n")       # semicolons: no tab, no whitespace
        fh.write("ACGT;1;2\n")
    with pytest.raises(GenomeWideDataError, match="delimiter sniff is wrong"):
        load_fragments(tmp_path, which="LB", allow_download=False)


def test_min_columns_is_a_real_floor():
    assert MIN_COLUMNS >= 8


def test_load_refuses_a_zero_fragment_parse(tmp_path):
    """Every data row rejected => refuse, never return an empty list for a caller to average."""
    p = tmp_path / "GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz"
    _write(p, [_row(seq="ACGTN" * 10)])                  # non-ACGT: all rows dropped
    with pytest.raises(GenomeWideDataError, match="ZERO fragments"):
        load_fragments(tmp_path, which="LB", allow_download=False)


def test_load_reports_what_it_dropped_rather_than_dropping_silently(tmp_path):
    p = tmp_path / "GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz"
    _write(p, [_row(), _row(seq="ACGTN" * 10), _row(e1="notanumber", ave="nope"), _row()])
    frags, st = load_fragments(tmp_path, which="LB", allow_download=False)
    assert len(frags) == 2
    assert st.n_non_acgt == 1
    assert st.n_unparseable_value == 1
    assert st.delimiter == "whitespace"
    assert st.n_header_columns == len(COLS)


def test_load_refuses_a_file_missing_a_required_column(tmp_path):
    p = tmp_path / "GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("fragment RNA_exp_ave start end strand extra1 extra2 extra3\n")
        fh.write("ACGT 1.0 1 4 + x y z\n")
    with pytest.raises(GenomeWideDataError, match="lacks required columns"):
        load_fragments(tmp_path, which="LB", allow_download=False)


def test_unknown_file_key_refuses():
    with pytest.raises(GenomeWideDataError, match="unknown file key"):
        load_fragments(Path("."), which="nonsense", allow_download=False)


# ---------------------------------------------------------------------------------------------------
# the noise ceiling
# ---------------------------------------------------------------------------------------------------
def test_noise_ceiling_is_spearman_brown_then_sqrt():
    """Perfectly reproducible replicates => reliability 1 => ceiling 1. The anchor."""
    f = [Fragment(seq="A", expression=v, rep1=v, rep2=v, start=0, end=1, strand="+")
         for v in (1., 2., 3., 4., 5.)]
    c = noise_ceiling(f)
    assert c["single_replicate_agreement"] == pytest.approx(1.0)
    assert c["ceiling_for_a_sequence_model"] == pytest.approx(1.0)


def test_noise_ceiling_steps_a_noisy_single_replicate_UP_for_the_mean():
    """Spearman-Brown: averaging 2 replicates is more reliable than either. So the ceiling on predicting
    the MEAN must exceed the single-replicate agreement -- which is exactly why comparing a model's rho
    against single-replicate agreement would understate the ceiling."""
    c = noise_ceiling(_frags(500))
    r11 = c["single_replicate_agreement"]
    assert 0.0 < r11 < 1.0, "the fixture must be genuinely noisy for this test to mean anything"
    assert c["reliability_of_2rep_mean_spearman_brown"] > r11
    assert c["ceiling_for_a_sequence_model"] > r11


# ---------------------------------------------------------------------------------------------------
# TRAP: the condition comparison must be NOISE-MATCHED
# ---------------------------------------------------------------------------------------------------
def test_condition_effect_verdict_uses_the_MATCHED_comparison_not_the_unmatched_one():
    """THE load-bearing test. On the real data the unmatched average-vs-average form reads 0.5025 against
    single-replicate agreement 0.458/0.477 -- i.e. it looks like condition does NOT matter, which is
    backwards. The verdict must key off the single-vs-single figures."""
    a = _frags(600, seed=1)
    # b is the SAME underlying truth, independently re-measured: condition is genuinely irrelevant here
    import random
    r = random.Random(9)
    b = [Fragment(seq=f.seq, expression=f.expression,
                  rep1=f.expression + r.gauss(0, .3), rep2=f.expression + r.gauss(0, .3),
                  start=f.start, end=f.end, strand=f.strand) for f in a]
    res = condition_effect(a, b)
    assert res["n_shared"] == 600
    # with condition irrelevant, cross and within should be comparable, so no signal is claimed
    assert res["verdict"] == "CONDITION_NOT_DISTINGUISHABLE_FROM_NOISE", res
    assert res["disattenuated_true_cross_condition"] > 0.8, res


def test_condition_effect_detects_a_REAL_condition_difference():
    """Non-vacuity: a genuinely different second condition must come out as carrying signal."""
    a = _frags(600, seed=2)
    import random
    r = random.Random(3)
    b = [Fragment(seq=f.seq, expression=r.random(),      # unrelated truth in condition b
                  rep1=0.0, rep2=0.0, start=f.start, end=f.end, strand=f.strand) for f in a]
    for f in b:
        f.rep1 = f.expression + r.gauss(0, .3)
        f.rep2 = f.expression + r.gauss(0, .3)
    res = condition_effect(a, b)
    assert res["verdict"] == "CONDITION_CARRIES_SIGNAL", res
    assert res["noise_matched_gap"] < 0


def test_condition_effect_refuses_an_underpowered_overlap():
    a, b = _frags(50), _frags(50, seed=5)
    assert condition_effect(a, b)["status"] == "insufficient_shared_fragments"


# ---------------------------------------------------------------------------------------------------
# TRAP: 91.5% of fragments overlap, so there is NO random split on offer
# ---------------------------------------------------------------------------------------------------
def test_position_blocked_split_is_disjoint_BY_CONSTRUCTION():
    f = _frags(500)
    tr, te, held = position_blocked_split(f, n_blocks=10, held_blocks=2)
    assert len(tr) + len(te) == len(f), "must partition, not sample"
    lo, hi = min(x.start for x in f), max(x.end for x in f)
    w = (hi - lo) / 10
    blk = lambda x: min(9, int(((x.start + x.end) / 2 - lo) / w))
    assert all(blk(x) not in held for x in tr)
    assert all(blk(x) in held for x in te)


def test_position_blocked_split_holds_out_CONTIGUOUS_territory():
    """The point is geography, not a sample: every held fragment's midpoint lies in a held block, so no
    near-duplicate of a training fragment can appear in test except across a block seam."""
    f = _frags(500)
    _tr, te, held = position_blocked_split(f, n_blocks=10, held_blocks=1, seed=3)
    assert len(held) == 1
    mids = sorted((x.start + x.end) / 2 for x in te)
    assert mids[-1] - mids[0] < (max(y.end for y in f) - min(y.start for y in f)) / 5


def test_position_blocked_split_is_seed_stable_and_validates_its_args():
    f = _frags(300)
    assert position_blocked_split(f, seed=7)[2] == position_blocked_split(f, seed=7)[2]
    with pytest.raises(ValueError, match="held_blocks must be"):
        position_blocked_split(f, n_blocks=5, held_blocks=5)
    with pytest.raises(ValueError, match="no fragments"):
        position_blocked_split([])


def test_there_is_deliberately_NO_random_split_in_this_module():
    """A random split leaks here (91.5% overlap measured), so the module must not offer one to reach for."""
    import dna_decode.glm.genomewide as gw
    assert not hasattr(gw, "random_split")


# ---------------------------------------------------------------------------------------------------
# features must be LENGTH-INVARIANT here
# ---------------------------------------------------------------------------------------------------
def test_kmer_features_are_length_invariant():
    """Fragments are 48-475 bp with 280 distinct lengths, so a fixed-width feature cannot apply."""
    v = kmer_features(["ACGT" * 12, "ACGT" * 100], k=4)
    assert len(v[0]) == len(v[1]) == 256
    assert sum(v[0]) == pytest.approx(1.0)
    assert sum(v[1]) == pytest.approx(1.0)


def test_positional_onehot_REFUSES_ragged_input_rather_than_padding():
    """REWRITTEN, not renumbered. This test originally asserted `positional_onehot` was ABSENT from the
    module, which was right while the module held only variable-length FRAGMENTS. Once peak TILES arrived
    the premise changed -- tiles are a fixed 150 bp, so the grid's winning feature genuinely applies to
    them. Leaving the old assertion in place (or deleting it) would have left the real invariant untested.

    The real invariant: it must work on fixed-length input and REFUSE ragged input, because silently
    padding sheared fragments to a common width would invent an anchor and produce a plausible,
    meaningless number.
    """
    from dna_decode.glm.genomewide import positional_onehot
    v = positional_onehot(["ACGT", "TGCA"])
    assert len(v) == 2 and len(v[0]) == 16 and sum(v[0]) == 4.0
    with pytest.raises(GenomeWideDataError, match="requires a FIXED length"):
        positional_onehot(["ACGT", "ACG"])


def test_positional_onehot_sees_position_where_a_kmer_histogram_cannot():
    """Why it is worth having for tiles at all: same composition, reversed layout."""
    from dna_decode.glm.genomewide import positional_onehot
    a, b = "AACC", "CCAA"
    assert kmer_features([a], k=1) == kmer_features([b], k=1)
    assert positional_onehot([a]) != positional_onehot([b])


def test_gc_features_width_and_range():
    v = gc_features(["AAAA", "GCGC", "ACGT"])
    assert [len(x) for x in v] == [1, 1, 1]
    assert v[0][0] == 0.0 and v[1][0] == 1.0 and v[2][0] == 0.5


# ---------------------------------------------------------------------------------------------------
# the real substrate
# ---------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not HAVE_REAL, reason="GSE144621 not cached on this host")
def test_the_REAL_frag_file_loads_with_the_measured_shape():
    frags, st = load_fragments(REAL.parent, which="LB", allow_download=False)
    assert st.n_parsed == 321123, st.as_dict()
    assert st.delimiter == "whitespace", "the frag files are SPACE-delimited (the peak-tile file is tab)"
    assert st.n_header_columns == 13
    assert st.length_min < 100 and st.length_max > 400, "variable length is the point"
    assert st.n_distinct_lengths > 200


@pytest.mark.skipif(not HAVE_REAL, reason="GSE144621 not cached on this host")
def test_the_REAL_expression_column_is_ALREADY_dna_normalised():
    """Trap 2, pinned on real data: if RNA_exp_ave were raw counts it would track DNA abundance strongly.
    Measured spearman is -0.039, so dividing by DNA_ave again would double-normalise."""
    import gzip as _g
    from dna_decode.glm.genomewide import spearman
    rna, dna = [], []
    with _g.open(REAL, "rt") as fh:
        hdr = next(fh).split()
        ix = {c: i for i, c in enumerate(hdr)}
        for i, line in enumerate(fh):
            q = line.split()
            if len(q) == len(hdr):
                rna.append(float(q[ix["RNA_exp_ave"]])); dna.append(float(q[ix["DNA_ave"]]))
            if i >= 20000:
                break
    assert abs(spearman(rna, dna)) < 0.15, "RNA_exp_ave must be independent of DNA_ave"
