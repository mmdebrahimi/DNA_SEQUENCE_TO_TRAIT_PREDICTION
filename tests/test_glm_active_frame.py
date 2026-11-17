"""Guards for `dna_decode/glm/active_frame.py`.

The load-bearing test here is `test_the_curve_is_flat_on_three_unrelated_distributions`. The
statistic this module replaces was retracted because its value was a closed form in the split
proportion and four unrelated distributions all returned it. So the replacement has to clear
the opposite bar: hold the design fixed, vary the data, and the curve must MOVE -- flat where
reliability really is constant, sloped where it really is not. A statistic that cannot come
out flat is measuring the design, and this file fails before any real data is read.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path

import numpy as np
import pytest

from dna_decode.glm.active_frame import (
    ActiveFrameError,
    SubReplicatePair,
    attenuation_null,
    frame_decomposition,
    join_sub_replicates,
    read_sub_replicates_with_stats,
    reproducibility_curve,
    restrict_predictions,
    shape_summary,
    spearman_brown_ceiling,
    sub_replicate_key,
    threshold_audit,
    tile_key,
)
from dna_decode.glm.genomewide import noise_ceiling

SUBSTRATE = Path("data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz")

# Design parameters held IDENTICAL across every synthetic case below. Varying the data while
# holding these fixed is the whole test.
N = 1200
N_BINS = 8
N_DRAWS = 30
SEED = 7


class _Tile:
    """Minimal stand-in carrying only the attributes this module touches."""

    def __init__(self, seq, peak, expression, rep1, rep2, active="active", strand="+"):
        self.seq = seq
        self.peak = peak
        self.expression = expression
        self.rep1 = rep1
        self.rep2 = rep2
        self.active = active
        self.strand = strand
        self.category = "tile"


def _three_column_sample(truth, noise_sd, rng):
    """One binning column and one scored PAIR, all noisy reads of the same latent truth.

    This mirrors the real design: `bin_by` stands in for group 1's summary and the pair for
    group 2's two sub-replicates, so the binning variable shares no noise with the pair.
    `noise_sd` may be a scalar or per-row.

    CAUTION, and `test_constant_additive_noise_on_a_heavy_tail_really_is_a_gradient` records
    why: a SCALAR `noise_sd` does NOT give constant reliability across the rank range unless
    the marginal is light-tailed. On a heavy-tailed truth the standard deviation is dominated
    by the tail, so one fixed absolute error swamps the bulk and barely touches the tail --
    a real reliability gradient, planted by accident. Use `_constant_reliability_sample` when
    constant reliability is what the test needs.
    """
    t = np.asarray(truth, dtype=float)
    sd = np.broadcast_to(np.asarray(noise_sd, dtype=float), t.shape)
    bin_by = t + rng.normal(0.0, sd)
    pair_a = t + rng.normal(0.0, sd)
    pair_b = t + rng.normal(0.0, sd)
    return bin_by, pair_a, pair_b


def _constant_reliability_sample(marginal, rng, *, a=0.80, b=0.85):
    """Three columns whose RANK reliability is constant across the range, on any marginal.

    Built on ranks rather than by adding an absolute error, which is the only way to hold
    reliability fixed while the marginal varies. A single shared latent factor gives one
    constant correlation everywhere; rank-mapping onto `marginal` then changes the
    distribution WITHOUT touching the dependence structure.

    Deliberately does NOT reuse the module's correlation-matrix or Spearman-to-Gaussian code:
    the construction here is a two-line shared-factor model, so a bug in the module's null
    cannot be cancelled out by the same bug in the data it is judged against.
    """
    m = np.asarray(marginal, dtype=float)
    n = m.size
    common = rng.standard_normal(n)
    z_v = a * common + math.sqrt(1 - a * a) * rng.standard_normal(n)
    z_a = b * common + math.sqrt(1 - b * b) * rng.standard_normal(n)
    z_b = b * common + math.sqrt(1 - b * b) * rng.standard_normal(n)
    srt = np.sort(m)
    out = []
    for z in (z_v, z_a, z_b):
        ranks = np.argsort(np.argsort(z, kind="stable"), kind="stable")
        out.append(srt[ranks])
    return out[0], out[1], out[2]


# ----------------------------------------------------------------------------------------
# THE decisive pair of tests
# ----------------------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["uniform", "gaussian", "heavy_tailed"])
def test_the_curve_is_flat_on_three_unrelated_distributions(name):
    """Design fixed, MARGINAL varied, reliability genuinely constant -> trend inside the band.

    The retracted statistic failed because four unrelated distributions all returned the same
    value. This is the mirror bar: three unrelated marginals, identical design parameters, and
    a dependence structure with no real gradient must all land inside the null band.

    HONEST SCOPE: the sample and the null are both rank-based, so this checks that the curve
    does not manufacture a trend from the design -- it does NOT prove the curve is unbiased,
    because a bias shared by the curve and its own null would cancel. The claim that the curve
    can MOVE is carried by `test_a_planted_reliability_gradient_beats_the_null_maximum` and by
    `test_constant_additive_noise_on_a_heavy_tail_really_is_a_gradient`, neither of which
    could pass if `_trend` were pinned near zero.
    """
    rng = np.random.default_rng(SEED)
    marginal = {
        "uniform": lambda: rng.uniform(0.0, 1.0, N),
        "gaussian": lambda: rng.normal(0.0, 1.0, N),
        "heavy_tailed": lambda: rng.pareto(2.0, N),
    }[name]()
    bin_by, pair_a, pair_b = _constant_reliability_sample(marginal, rng)

    curve = reproducibility_curve(
        bin_by, pair_a, pair_b, n_bins=N_BINS, binning_is_independent=True,
    )
    null = attenuation_null(
        bin_by, pair_a, pair_b, n_bins=N_BINS, n_draws=N_DRAWS, seed=SEED,
    )
    lo, hi = null["trend_null"]["min"], null["trend_null"]["max"]
    assert lo <= curve["trend"] <= hi, (
        f"{name}: observed trend {curve['trend']} fell OUTSIDE the null band [{lo}, {hi}] "
        f"even though reliability is constant by construction. Either the curve is picking "
        f"up the design or the null does not hold the design fixed"
    )


def test_constant_additive_noise_on_a_heavy_tail_really_is_a_gradient():
    """A trap this file hit while being written, kept as a test because it is instructive.

    Adding ONE fixed-variance error to a Pareto truth looks homoscedastic and is not: the
    standard deviation is set by the tail, so the same absolute error swamps the bulk and
    barely touches the tail. Reliability therefore climbs with the value, and the curve is
    CORRECT to clear the null -- the first draft of the flatness test above called this a
    failure of the curve when it was a failure of the generator.
    """
    rng = np.random.default_rng(SEED)
    truth = rng.pareto(2.0, N)
    bin_by, pair_a, pair_b = _three_column_sample(truth, 0.6 * float(np.std(truth)), rng)

    curve = reproducibility_curve(
        bin_by, pair_a, pair_b, n_bins=N_BINS, binning_is_independent=True,
    )
    null = attenuation_null(
        bin_by, pair_a, pair_b, n_bins=N_BINS, n_draws=N_DRAWS, seed=SEED,
    )
    assert curve["trend"] > null["trend_null"]["max"], (
        f"trend {curve['trend']} did not clear the null max {null['trend_null']['max']}; the "
        f"documented heavy-tail gradient should be visible"
    )


def test_a_planted_reliability_gradient_beats_the_null_maximum():
    """The other half: when reliability really DOES move, the curve must see it.

    Without this, a curve that always returns ~0 would pass the flatness test and be useless.
    """
    rng = np.random.default_rng(SEED)
    truth = rng.normal(0.0, 1.0, N)
    # noise SHRINKS as the truth rises, so agreement must climb with the binning variable
    sd = np.interp(truth, (truth.min(), truth.max()), (1.6, 0.12))
    bin_by, pair_a, pair_b = _three_column_sample(truth, sd, rng)

    curve = reproducibility_curve(
        bin_by, pair_a, pair_b, n_bins=N_BINS, binning_is_independent=True,
    )
    null = attenuation_null(
        bin_by, pair_a, pair_b, n_bins=N_BINS, n_draws=N_DRAWS, seed=SEED,
    )
    assert curve["trend"] > null["trend_null"]["max"], (
        f"planted gradient: trend {curve['trend']} did not clear the null MAX "
        f"{null['trend_null']['max']}, so the curve cannot detect a real gradient"
    )
    assert curve["trend"] > 0.5, f"expected a strong positive trend, got {curve['trend']}"


def test_the_null_reports_that_it_holds_the_design_fixed():
    rng = np.random.default_rng(SEED)
    truth = rng.normal(0.0, 1.0, N)
    bin_by, pair_a, pair_b = _three_column_sample(truth, 0.6, rng)
    null = attenuation_null(bin_by, pair_a, pair_b, n_bins=N_BINS, n_draws=5, seed=SEED)
    held = null["design_held_fixed"]
    assert "CONSTANT by construction" in held["reliability_across_the_range"]
    assert held["marginals"].startswith("observed")
    # the null must reproduce the observed global correlations it was built from
    from dna_decode.glm.genomewide import spearman
    assert held["pooled_spearman_pair"] == pytest.approx(
        round(spearman(list(pair_a), list(pair_b)), 4), abs=1e-9)


# ----------------------------------------------------------------------------------------
# the refusals that keep the curve honest
# ----------------------------------------------------------------------------------------

def test_curve_refuses_the_same_object_as_binning_variable():
    rng = np.random.default_rng(SEED)
    a = rng.normal(0, 1, N)
    b = rng.normal(0, 1, N)
    with pytest.raises(ActiveFrameError, match="SAME OBJECT"):
        reproducibility_curve(a, a, b, n_bins=N_BINS, binning_is_independent=True)


def test_curve_refuses_an_element_wise_identical_binning_variable():
    rng = np.random.default_rng(SEED)
    a = rng.normal(0, 1, N)
    b = rng.normal(0, 1, N)
    with pytest.raises(ActiveFrameError, match="ELEMENT-WISE IDENTICAL"):
        reproducibility_curve(a.copy(), a, b, n_bins=N_BINS, binning_is_independent=True)


def test_binning_is_independent_has_no_default_so_the_caller_must_state_it():
    rng = np.random.default_rng(SEED)
    v, a, b = _three_column_sample(rng.normal(0, 1, N), 0.6, rng)
    with pytest.raises(TypeError):
        reproducibility_curve(v, a, b, n_bins=N_BINS)  # type: ignore[call-arg]


def test_declaring_a_confounded_binning_withholds_the_trend():
    """A confounded curve may be inspected; its slope may not be quoted."""
    rng = np.random.default_rng(SEED)
    v, a, b = _three_column_sample(rng.normal(0, 1, N), 0.6, rng)
    out = reproducibility_curve(v, a, b, n_bins=N_BINS, binning_is_independent=False)
    assert out["trend"] is None
    assert "WITHHELD" in out["reads_as"]
    assert out["conditions_on_the_outcome"] is True
    assert len(out["bins"]) == N_BINS


def test_curve_refuses_too_few_rows_for_the_requested_bins():
    rng = np.random.default_rng(SEED)
    v, a, b = _three_column_sample(rng.normal(0, 1, 10), 0.6, rng)
    with pytest.raises(ActiveFrameError, match="cannot support"):
        reproducibility_curve(v, a, b, n_bins=10, binning_is_independent=True)


def test_curve_refuses_mismatched_lengths():
    rng = np.random.default_rng(SEED)
    v, a, b = _three_column_sample(rng.normal(0, 1, N), 0.6, rng)
    with pytest.raises(ActiveFrameError, match="length mismatch"):
        reproducibility_curve(v[:-1], a, b, n_bins=N_BINS, binning_is_independent=True)


# ----------------------------------------------------------------------------------------
# one definition of the ceiling arithmetic, not two
# ----------------------------------------------------------------------------------------

@pytest.mark.parametrize("r", [0.0, 0.1, 0.3927, 0.7230, 0.7283, 0.7945, 0.8849, 0.99])
def test_spearman_brown_ceiling_agrees_with_noise_ceiling_to_the_decimal(r):
    """The two doors into the same arithmetic must not drift apart."""
    rng = np.random.default_rng(3)
    n = 4000
    # build a pair with approximately the target correlation, then compare BOTH routes on it
    t = rng.normal(0, 1, n)
    e1, e2 = rng.normal(0, 1, n), rng.normal(0, 1, n)
    w = math.sqrt(max(r, 1e-9) / max(1 - r, 1e-9)) if r < 1 else 50.0
    rep1, rep2 = t * w + e1, t * w + e2
    tiles = [_Tile(f"S{i}", f"p{i}", 0.0, float(rep1[i]), float(rep2[i])) for i in range(n)]

    from dna_decode.glm.genomewide import spearman as _sp

    via_tiles = noise_ceiling(tiles)
    # feed the RAW correlation, not noise_ceiling's rounded `single_replicate_agreement`:
    # the step-up is computed from the unrounded value and rounded once, so passing the
    # rounded value through it rounds twice and can shift the last digit by 1e-4. That is a
    # property of the rounding, not a drift between the two routes, and comparing them on
    # the same raw input is what actually pins the arithmetic.
    raw = _sp([t.rep1 for t in tiles], [t.rep2 for t in tiles])
    via_bare = spearman_brown_ceiling(raw)
    assert via_bare["single_replicate_agreement"] == via_tiles["single_replicate_agreement"]
    assert via_bare["reliability_of_2rep_mean_spearman_brown"] == pytest.approx(
        via_tiles["reliability_of_2rep_mean_spearman_brown"], abs=1e-9)
    assert via_bare["ceiling_for_a_sequence_model"] == pytest.approx(
        via_tiles["ceiling_for_a_sequence_model"], abs=1e-9)


def test_spearman_brown_reproduces_the_published_anchor_by_hand():
    """0.7945 -> 0.941 is the published ceiling; the arithmetic must land on it."""
    out = spearman_brown_ceiling(0.7945)
    assert out["ceiling_for_a_sequence_model"] == pytest.approx(0.941, abs=0.001)


# ----------------------------------------------------------------------------------------
# the join: the sequence is NOT a key
# ----------------------------------------------------------------------------------------

def _tiles_and_subs(n=6):
    tiles = [_Tile(f"SEQ{i}", f"{i}_{i + 150}_+", float(i), float(i), float(i)) for i in range(n)]
    subs = {tile_key(t): (1.0, 2.0, 3.0, 4.0) for t in tiles}
    return tiles, subs


def test_join_succeeds_and_preserves_tile_order():
    tiles, subs = _tiles_and_subs()
    pairs = join_sub_replicates(tiles, subs)
    assert [p.key for p in pairs] == [tile_key(t) for t in tiles]
    assert pairs[0].g1_mean == pytest.approx(1.5)
    assert pairs[0].g2_mean == pytest.approx(3.5)


def test_join_raises_on_a_dropped_key():
    tiles, subs = _tiles_and_subs()
    subs.pop(tile_key(tiles[2]))
    with pytest.raises(ActiveFrameError, match="key sets DIFFER"):
        join_sub_replicates(tiles, subs)


def test_join_raises_on_an_extra_key():
    tiles, subs = _tiles_and_subs()
    subs[("EXTRA", "9_9_+")] = (1.0, 1.0, 1.0, 1.0)
    with pytest.raises(ActiveFrameError, match="key sets DIFFER"):
        join_sub_replicates(tiles, subs)


def test_join_raises_on_a_shifted_key():
    """A shift keeps the COUNT identical and changes the pairing -- the dangerous case."""
    tiles, subs = _tiles_and_subs()
    keys = [tile_key(t) for t in tiles]
    shifted = {keys[i]: subs[keys[(i + 1) % len(keys)]] for i in range(len(keys))}
    # a pure value shift on identical keys cannot be caught by key comparison, so the guard
    # that matters is that a shifted KEY SET is caught; rename one key to simulate it
    shifted.pop(keys[0])
    shifted[("SHIFTED", "0_150_+")] = (1.0, 1.0, 1.0, 1.0)
    with pytest.raises(ActiveFrameError, match="key sets DIFFER"):
        join_sub_replicates(tiles, shifted)


def test_a_duplicated_SEQUENCE_with_distinct_peaks_joins_correctly():
    """688 sequences recur on the real substrate; the key is the PAIR, not the sequence."""
    tiles = [
        _Tile("SAME", "1_151_+", 1.0, 1.0, 1.0),
        _Tile("SAME", "2_152_+", 2.0, 2.0, 2.0),
        _Tile("SAME", "3_153_-", 3.0, 3.0, 3.0),
    ]
    subs = {tile_key(t): (float(i), float(i), float(i), float(i)) for i, t in enumerate(tiles)}
    assert len({t.seq for t in tiles}) == 1, "the test is vacuous unless the sequence repeats"
    pairs = join_sub_replicates(tiles, subs)
    assert [p.g1_1 for p in pairs] == [0.0, 1.0, 2.0]


def test_join_refuses_when_the_TILE_side_has_a_duplicate_key():
    tiles = [_Tile("SAME", "1_151_+", 1.0, 1.0, 1.0)] * 2
    subs = {tile_key(tiles[0]): (1.0, 1.0, 1.0, 1.0)}
    with pytest.raises(ActiveFrameError, match="share a .seq, peak. key"):
        join_sub_replicates(tiles, subs)


def test_sub_replicate_key_matches_the_loaders_own_construction():
    """`Tile.peak` is f"{start}_{end}_{strand}" and `Tile.seq` is upper-cased by the loader."""
    t = _Tile("acgt", "10_160_-", 0.0, 0.0, 0.0)
    t.seq = "ACGT"
    assert sub_replicate_key("acgt", "10", "160", "-") == tile_key(t)


# ----------------------------------------------------------------------------------------
# shape_summary must NOT regrow the retracted field
# ----------------------------------------------------------------------------------------

def test_shape_summary_has_no_threshold_spearman_field():
    """Asserted by KEY SET, so a reintroduction under any name fails this test."""
    rng = np.random.default_rng(SEED)
    out = shape_summary(rng.normal(0, 1, 500))
    assert set(out.keys()) == {
        "n", "n_bins", "min", "max", "median", "modal_bin_share",
        "n_distinct_values", "distinct_value_share", "conditions_on_the_outcome",
        "deliberately_absent",
    }
    blob = " ".join(k.lower() for k in out.keys())
    assert "threshold" not in blob and "ceiling" not in blob


def test_shape_summary_separates_a_degenerate_from_a_spread_sample():
    rng = np.random.default_rng(SEED)
    spread = shape_summary(rng.normal(0, 1, 2000))
    piled = shape_summary(np.r_[np.zeros(1800), rng.normal(0, 1, 200)])
    assert piled["modal_bin_share"] > spread["modal_bin_share"]
    assert piled["n_distinct_values"] < spread["n_distinct_values"]


def test_shape_summary_handles_an_all_identical_sample_without_claiming_structure():
    out = shape_summary(np.zeros(100))
    assert out["n_distinct_values"] == 1
    assert out["modal_bin_share"] == 1.0


def test_shape_summary_refuses_an_empty_sample():
    with pytest.raises(ActiveFrameError, match="no finite values"):
        shape_summary([])


# ----------------------------------------------------------------------------------------
# the active flag, measured
# ----------------------------------------------------------------------------------------

def test_threshold_audit_recovers_a_planted_clean_cut():
    tiles = (
        [_Tile(f"a{i}", f"{i}_x_+", 0.1 * i, 0, 0, active="inactive") for i in range(1, 6)]
        + [_Tile(f"b{i}", f"{i}_y_+", 1.0 + 0.1 * i, 0, 0, active="active") for i in range(1, 6)]
    )
    out = threshold_audit(tiles)
    assert out["overlap_either_side"] == 0
    assert out["flag_is_a_clean_cut_on_the_label"] is True
    assert out["implied_cut"] == pytest.approx(1.1)
    assert "RANGE-RESTRICTED" in out["reads_as"]


def test_threshold_audit_reports_overlap_and_refuses_to_claim_a_cut():
    tiles = (
        [_Tile(f"a{i}", f"{i}_x_+", 0.1 * i, 0, 0, active="inactive") for i in range(1, 6)]
        + [_Tile(f"b{i}", f"{i}_y_+", 1.0 + 0.1 * i, 0, 0, active="active") for i in range(1, 6)]
    )
    tiles.append(_Tile("rogue", "99_z_+", 0.05, 0, 0, active="active"))
    out = threshold_audit(tiles)
    assert out["overlap_either_side"] > 0
    assert out["flag_is_a_clean_cut_on_the_label"] is False
    assert "LABEL rather than a threshold" in out["reads_as"]


def test_threshold_audit_refuses_a_single_class():
    tiles = [_Tile(f"a{i}", f"{i}_x_+", float(i), 0, 0, active="active") for i in range(5)]
    with pytest.raises(ActiveFrameError, match="needs both classes"):
        threshold_audit(tiles)


# ----------------------------------------------------------------------------------------
# frame_decomposition stamps the conditioned blocks
# ----------------------------------------------------------------------------------------

def test_frame_decomposition_stamps_selection_attenuated_on_both_conditioned_blocks():
    rng = np.random.default_rng(SEED)
    n = 400
    t = rng.normal(0, 1, n)
    tiles = []
    for i in range(n):
        e = float(t[i])
        tiles.append(_Tile(
            f"s{i}", f"{i}_p_+", e,
            e + float(rng.normal(0, 0.5)), e + float(rng.normal(0, 0.5)),
            active="active" if e > 0 else "inactive",
        ))
    out = frame_decomposition(tiles)
    assert out["pooled"]["conditions_on_the_outcome"] is False
    assert "selection_attenuated" not in out["pooled"]
    for k in ("active", "inactive"):
        assert out[k]["selection_attenuated"] is True
        assert out[k]["conditions_on_the_outcome"] is True
        assert "reproducibility_curve" in out[k]["see_instead"]


def test_frame_decomposition_pooled_block_reproduces_noise_ceiling_exactly():
    rng = np.random.default_rng(SEED)
    n = 300
    t = rng.normal(0, 1, n)
    tiles = [
        _Tile(f"s{i}", f"{i}_p_+", float(t[i]),
              float(t[i] + rng.normal(0, 0.4)), float(t[i] + rng.normal(0, 0.4)),
              active="active" if i % 2 else "inactive")
        for i in range(n)
    ]
    out = frame_decomposition(tiles)
    direct = noise_ceiling(tiles)
    assert out["pooled"]["ceiling_for_a_sequence_model"] == direct["ceiling_for_a_sequence_model"]
    assert out["pooled"]["single_replicate_agreement"] == direct["single_replicate_agreement"]


# ----------------------------------------------------------------------------------------
# restrict_predictions fails loudly rather than scoring a shifted vector
# ----------------------------------------------------------------------------------------

def test_restrict_predictions_keeps_only_the_subset_and_pairs_correctly():
    tiles = [
        _Tile(f"s{i}", f"{i}_p_+", float(i), 0, 0, active="active" if i % 2 else "inactive")
        for i in range(6)
    ]
    y, p = restrict_predictions([0, 1, 2, 3, 4, 5], [10, 11, 12, 13, 14, 15], tiles,
                                lambda t: t.active == "active")
    assert y == [1.0, 3.0, 5.0]
    assert p == [11.0, 13.0, 15.0]


def test_restrict_predictions_raises_on_an_out_of_range_index():
    tiles = [_Tile(f"s{i}", f"{i}_p_+", float(i), 0, 0) for i in range(3)]
    with pytest.raises(ActiveFrameError, match="outside the 3 loaded tiles"):
        restrict_predictions([0, 1, 7], [1, 2, 3], tiles, lambda t: True)


def test_restrict_predictions_raises_on_a_length_mismatch():
    tiles = [_Tile(f"s{i}", f"{i}_p_+", float(i), 0, 0) for i in range(3)]
    with pytest.raises(ActiveFrameError, match="disagree"):
        restrict_predictions([0, 1, 2], [1, 2], tiles, lambda t: True)


# ----------------------------------------------------------------------------------------
# real data, gated
# ----------------------------------------------------------------------------------------

@pytest.mark.skipif(not SUBSTRATE.exists(), reason=f"absent substrate: {SUBSTRATE}")
def test_real_substrate_reproduces_every_published_anchor():
    from dna_decode.glm.genomewide import load_peak_tiles, spearman

    loaded, _st = load_peak_tiles("data/mpra", allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    assert len(tiles) == 44106

    audit = threshold_audit(tiles)
    assert audit["n_active"] == 14913
    assert audit["n_inactive"] == 29193
    assert audit["overlap_either_side"] == 0
    assert audit["implied_cut"] == pytest.approx(0.5778, abs=1e-4)

    dec = frame_decomposition(tiles)
    assert dec["pooled"]["single_replicate_agreement"] == pytest.approx(0.7945, abs=1e-3)
    assert dec["pooled"]["ceiling_for_a_sequence_model"] == pytest.approx(0.9410, abs=1e-3)

    subs, st = read_sub_replicates_with_stats(SUBSTRATE)
    assert st.n_kept == 44106, f"expected all 44,106 rows numeric on all four columns, got {st.n_kept}"
    assert st.n_unparseable == 0

    pairs = join_sub_replicates(tiles, subs)
    assert len(pairs) == 44106

    # the group summary IS the mean of its own two sub-columns -- verified, not assumed
    g1 = [p.g1_mean for p in pairs]
    g2 = [p.g2_mean for p in pairs]
    assert all(abs(p.g1_mean - (p.g1_1 + p.g1_2) / 2) < 1e-9 for p in pairs[:2000])
    assert [round(t.rep1, 6) for t in tiles[:2000]] == [round(v, 6) for v in g1[:2000]]
    assert [round(t.rep2, 6) for t in tiles[:2000]] == [round(v, 6) for v in g2[:2000]]

    within_1 = spearman([p.g1_1 for p in pairs], [p.g1_2 for p in pairs])
    within_2 = spearman([p.g2_1 for p in pairs], [p.g2_2 for p in pairs])
    assert within_1 == pytest.approx(0.7283, abs=1e-3)
    assert within_2 == pytest.approx(0.7230, abs=1e-3)
    assert spearman(g1, g2) == pytest.approx(0.7945, abs=1e-3)


@pytest.mark.skipif(not SUBSTRATE.exists(), reason=f"absent substrate: {SUBSTRATE}")
def test_the_real_curve_bins_on_one_group_and_scores_the_other():
    """End-to-end on real data: the binning variable shares no measurement with the pair."""
    from dna_decode.glm.genomewide import load_peak_tiles

    loaded, _st = load_peak_tiles("data/mpra", allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    pairs = join_sub_replicates(tiles, read_sub_replicates_with_stats(SUBSTRATE)[0])

    curve = reproducibility_curve(
        [p.g1_mean for p in pairs], [p.g2_1 for p in pairs], [p.g2_2 for p in pairs],
        n_bins=10, binning_is_independent=True,
    )
    assert curve["n"] == 44106
    assert len(curve["bins"]) == 10
    assert sum(b["n"] for b in curve["bins"]) == 44106
    assert curve["trend"] is not None
    assert not any(b["thin"] for b in curve["bins"])


# --------------------------------------------------------------------------------------
# ownership: this module is a CONSUMER of genomewide, never a second reader
# --------------------------------------------------------------------------------------
# Added 2026-11-16. `tests/test_substrate_loader_ownership.py` flagged active_frame as a second
# reader of GSE144621 and it was right: the reader hardcoded `split("\t")` where the owner has
# `sniff_delimiter` for this deposit's space-vs-tab trap, and its row-admission predicate differed
# from `load_peak_tiles` while `join_sub_replicates` requires their key sets to match EXACTLY.
# The reader moved into the owner. That module is now on this substrate's ALLOWED_MENTIONS list,
# which is a PER-MODULE exemption -- so if file-reading code is ever added back here, the ownership
# guard will NOT catch it. These tests are what covers that hole, and they work at the SYNTAX level
# because a substring search would match the docstrings that explain the arrangement.

AF_SOURCE = Path(__file__).resolve().parent.parent / "dna_decode" / "glm" / "active_frame.py"


def _af_tree():
    return ast.parse(AF_SOURCE.read_text(encoding="utf-8"))


def test_active_frame_imports_no_file_reading_machinery():
    """Not a substring check: the docstrings legitimately say the word `gzip`-adjacent things."""
    offenders = []
    for node in ast.walk(_af_tree()):
        if isinstance(node, ast.Import):
            offenders += [a.name for a in node.names if a.name.split(".")[0] in {"gzip", "csv"}]
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in {"gzip", "csv"}:
                offenders.append(node.module)
    assert not offenders, (
        f"active_frame imports file-reading machinery {offenders}. It is a CONSUMER of "
        f"genomewide, which owns GSE144621; a reader here re-opens the space-vs-tab trap and "
        f"forks the row-admission predicate that join_sub_replicates depends on."
    )


def test_active_frame_calls_no_open_and_no_gzip_open():
    """`open(...)` / `gzip.open(...)` anywhere in this module is the regression."""
    bad = []
    for node in ast.walk(_af_tree()):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        if isinstance(f, ast.Name) and f.id == "open":
            bad.append("open(")
        if isinstance(f, ast.Attribute) and f.attr == "open":
            base = f.value
            if isinstance(base, ast.Name):
                bad.append(f"{base.id}.open(")
    assert not bad, f"active_frame opens files directly: {sorted(set(bad))}"


def test_active_frame_actually_reaches_the_owner_for_the_read():
    """Non-vacuity for the two tests above: 'opens no file' must not mean 'reads nothing'.

    A module that had simply dropped the capability would pass both guards above, so pin that
    the call still resolves INTO genomewide and that the four column names are the owner's.
    """
    from dna_decode.glm import genomewide

    assert hasattr(genomewide, "read_sub_replicates_with_stats")
    assert genomewide.SUB_REPLICATE_COLUMNS == (
        "RNA_exp_sum_1_1", "RNA_exp_sum_1_2", "RNA_exp_sum_2_1", "RNA_exp_sum_2_2",
    )
    # the names active_frame re-exports must BE the owner's objects, not copies
    from dna_decode.glm import active_frame

    assert active_frame.SUB_REPLICATE_COLUMNS is genomewide.SUB_REPLICATE_COLUMNS
    assert active_frame.sub_replicate_key is genomewide.sub_replicate_key
    assert active_frame.SubReplicateReadStats is genomewide.SubReplicateReadStats
