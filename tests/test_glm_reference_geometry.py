"""Guards for the reference-geometry triage and the coordinate-key fix.

Both were produced by an adversarial review of the GLM decomposition. The review's claim about the
sequence-only key was CONFIRMED by measurement but its stated magnitude was wrong, and that distinction is
pinned here so neither half gets lost.
"""
from __future__ import annotations

import random
from pathlib import Path

import pytest

from dna_decode.glm.genomewide import Fragment, condition_effect

CACHE = Path("D:/dna_decode_cache/mpra")
HAVE_DATA = (CACHE / "GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz").exists()


def _frags(n=400, seed=0, shift=0):
    r = random.Random(seed)
    out = []
    for i in range(n):
        t = r.random()
        out.append(Fragment(seq="".join(r.choice("ACGT") for _ in range(60)),
                            expression=t, rep1=t + r.gauss(0, .2), rep2=t + r.gauss(0, .2),
                            start=i * 100 + shift, end=i * 100 + 60 + shift, strand="+"))
    return out


# ---------------------------------------------------------------------------------------------------
# the coordinate-key fix
# ---------------------------------------------------------------------------------------------------
def test_coordinate_key_is_the_DEFAULT():
    """A sequence-only key collapses a repeated genomic element to one arbitrary locus. The default must be
    the safe one, so a caller who does not know about the trap is not exposed to it."""
    a = _frags(300, seed=1)
    b = [Fragment(seq=f.seq, expression=f.expression, rep1=f.rep1, rep2=f.rep2,
                  start=f.start, end=f.end, strand=f.strand) for f in a]
    assert condition_effect(a, b)["key"] == "coordinate"
    assert condition_effect(a, b, coordinate_key=False)["key"] == "sequence_only"


def test_a_repeated_sequence_at_TWO_loci_collapses_under_the_sequence_key_only():
    """THE mechanism, on a synthetic repeat. Under the coordinate key both loci survive as distinct
    fragments; under the sequence key one is silently discarded."""
    dup = "ACGT" * 15
    base = _frags(200, seed=3)
    # same sequence, two different loci -- a transposon-like repeat
    rep = [Fragment(seq=dup, expression=0.5, rep1=0.5, rep2=0.5, start=s, end=s + 60, strand="+")
           for s in (10_000, 2_000_000)]
    a = base + rep
    b = [Fragment(seq=f.seq, expression=f.expression, rep1=f.rep1, rep2=f.rep2,
                  start=f.start, end=f.end, strand=f.strand) for f in a]
    coord = condition_effect(a, b, coordinate_key=True)
    seqonly = condition_effect(a, b, coordinate_key=False)
    assert coord["n_key_collisions_dropped"]["a"] == 0, "distinct loci must not collide on a coordinate key"
    assert seqonly["n_key_collisions_dropped"]["a"] == 1, "the repeat must register as a collision"
    assert coord["n_shared"] == seqonly["n_shared"] + 1


def test_collision_count_is_REPORTED_not_just_avoided():
    """The original defect was silence, not just the wrong key. The count must surface either way so the
    condition can never go unobserved again."""
    r = condition_effect(_frags(200), _frags(200), coordinate_key=False)
    assert "n_key_collisions_dropped" in r and "key" in r


@pytest.mark.skipif(not HAVE_DATA, reason="GSE144621 not cached on this host")
def test_on_REAL_data_the_key_changes_n_but_NOT_the_verdict():
    """The review claimed the sequence-only key could 'corrupt n_shared' and 'mispair condition labels'.
    Measured: 344 LB / 306 M9 rows collide, n moves 297,599 -> 297,868, and the noise-matched gap is
    -0.0614 under BOTH keys. So the defect is real and the stated magnitude was not -- both halves pinned.
    """
    from dna_decode.glm.genomewide import load_fragments
    lb, _ = load_fragments(CACHE, which="LB", allow_download=False)
    m9, _ = load_fragments(CACHE, which="M9", allow_download=False)
    coord = condition_effect(lb, m9, coordinate_key=True)
    seqonly = condition_effect(lb, m9, coordinate_key=False)
    assert coord["n_shared"] == 297868 and seqonly["n_shared"] == 297599
    assert seqonly["n_key_collisions_dropped"] == {"a": 344, "b": 306}
    assert coord["n_key_collisions_dropped"] == {"a": 0, "b": 0}
    # the headline is unchanged -- the fix is prospective (split integrity), not retrospective
    assert abs(coord["noise_matched_gap"] - seqonly["noise_matched_gap"]) < 1e-3
    assert coord["verdict"] == seqonly["verdict"] == "CONDITION_CARRIES_SIGNAL"


# ---------------------------------------------------------------------------------------------------
# the triage's pure logic
# ---------------------------------------------------------------------------------------------------
def _load_script():
    import importlib.util
    spec = importlib.util.spec_from_file_location("rg", Path("scripts/glm_reference_geometry.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_bands_are_monotone_and_cover_the_range():
    rg = _load_script()
    assert rg._band(0.95) == "DISTINCT"
    assert rg._band(0.80) == "DISTINCT"
    assert rg._band(0.70) == "PARTIAL"
    assert rg._band(0.65) == "PARTIAL"
    assert rg._band(0.50) == "NOT_MEANINGFULLY_DIFFERENT"


def test_the_validity_gate_is_evaluated_BEFORE_any_band():
    """Decision D1 in code: a broken instrument must not emit a graded number. A set separable from itself
    means the discriminator fabricates signal, and no pairwise AUROC is interpretable."""
    rg = _load_script()
    self_ctrl = {"grid": {"kmer": {"max": 0.95, "median": 0.93, "values": []}}}
    pairwise = {"grid_vs_tile": {"kmer": {"median": 1.0}}}
    v, d = rg._verdict(False, self_ctrl, pairwise, ["kmer"])
    assert v == "INDETERMINATE_INSTRUMENT_SEPARATES_A_SET_FROM_ITSELF"
    assert "Nothing graded" in d
    assert "DISTINCT" not in v, "a failed gate must not produce a band"


def test_a_clean_gate_DOES_produce_a_band():
    """Non-vacuity for the test above: the gate must not block everything."""
    rg = _load_script()
    sc = {"grid": {"kmer": {"max": 0.51, "median": 0.50, "values": []}}}
    pw = {"grid_vs_tile": {"kmer": {"median": 1.0}}, "tile_vs_natural": {"kmer": {"median": 0.62}}}
    v, d = rg._verdict(True, sc, pw, ["kmer"])
    assert v == "DISTINCT" and "1.0000" in d


def test_verdict_refuses_when_the_primary_pair_is_missing():
    rg = _load_script()
    sc = {"grid": {"kmer": {"max": 0.5, "median": 0.5, "values": []}}}
    v, _ = rg._verdict(True, sc, {"tile_vs_natural": {"kmer": {"median": 0.6}}}, ["kmer"])
    assert v == "INDETERMINATE_PRIMARY_PAIR_MISSING"


# ---------------------------------------------------------------------------------------------------
# the mechanism measurement
# ---------------------------------------------------------------------------------------------------
def test_diversity_separates_a_COMBINATORIAL_set_from_a_random_one():
    """The measured mechanism, on synthetic data. A set recombined from a small element vocabulary shares
    long identical blocks between members while having NO constant scaffold -- which is exactly the
    signature found on the real designed grid (median 13bp shared block, 0bp common prefix/suffix)."""
    rg = _load_script()
    r = random.Random(5)
    vocab = ["".join(r.choice("ACGT") for _ in range(30)) for _ in range(6)]
    combo = ["".join(r.choice(vocab) for _ in range(5)) for _ in range(500)]      # 150bp from 6 elements
    rand = ["".join(r.choice("ACGT") for _ in range(150)) for _ in range(500)]
    dc, dr = rg.sequence_diversity(combo), rg.sequence_diversity(rand)
    assert dc["median_longest_shared_block_bp"] > dr["median_longest_shared_block_bp"] + 10, (
        f"combinatorial {dc['median_longest_shared_block_bp']} vs random "
        f"{dr['median_longest_shared_block_bp']}")
    assert dc["longest_common_prefix_bp"] == 0, "recombination gives NO fixed scaffold -- that is the point"
    assert dr["mean_per_position_entropy_bits"] > 1.9


def test_diversity_reports_a_real_scaffold_when_there_IS_one():
    """Non-vacuity the other way: if a set DOES have a constant prefix, it must be reported. This is the
    hypothesis that was tested against the real grid and refuted."""
    rg = _load_script()
    r = random.Random(6)
    scaffolded = ["GGGGGGGGGG" + "".join(r.choice("ACGT") for _ in range(140)) for _ in range(300)]
    d = rg.sequence_diversity(scaffolded)
    assert d["longest_common_prefix_bp"] == 10
    assert d["n_invariant_positions"] >= 10
