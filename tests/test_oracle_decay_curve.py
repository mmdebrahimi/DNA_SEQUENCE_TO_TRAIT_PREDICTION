"""Tests for the oracle decay-curve harness (`scripts/oracle_decay_curve.py`).

All offline and torch-free: every test here exercises the PURE helpers. The torch path
(`load_finetuned_teacher`) is deliberately untested from here -- it needs a 22 MB checkpoint and the
atlas architecture, and it already has its own coverage via the degeneracy probe.

The properties worth pinning are the ones that could break SILENTLY and still produce a
publishable-looking curve:

  * chunking is a memory optimisation, so it must not change a single distance;
  * the per-bin floor must be exact at the boundary, because an off-by-one there turns a refused bin
    into a reported correlation;
  * a refused bin must not be able to contribute to the monotonicity verdict, which is the specific
    way a decay "trend" could be assembled out of bins that each declined to be measured.
"""
from __future__ import annotations

import numpy as np
import pytest

from scripts.oracle_decay_curve import (
    BINS,
    MIN_BIN_N,
    MIN_N_FOR_CI,
    SEQ_LENGTH,
    bin_index,
    encode,
    fisher_z_ci,
    hamming_min_to_set,
    monotonicity_report,
    summarise_bins,
)


# --------------------------------------------------------------------------- encode


def test_encode_shape_is_always_fixed_length():
    out = encode(["ACGT", "A" * (SEQ_LENGTH + 50)])
    assert out.shape == (2, SEQ_LENGTH)
    assert out.dtype == np.uint8


def test_non_acgt_becomes_a_fifth_symbol_that_always_mismatches():
    """An N must never silently read as a match against a real base."""
    n_row = encode(["N" * SEQ_LENGTH])
    for base in "ACGT":
        d = hamming_min_to_set(n_row, encode([base * SEQ_LENGTH]))
        assert d[0] == SEQ_LENGTH, "N matched %s somewhere" % base


def test_short_sequence_pads_with_the_mismatching_symbol_not_with_adenine():
    """A 4-base sequence padded with A would read as distance 0 from a poly-A train row."""
    d = hamming_min_to_set(encode(["AAAA"]), encode(["A" * SEQ_LENGTH]))
    assert d[0] == SEQ_LENGTH - 4


# ------------------------------------------------------------- hamming_min_to_set


def test_exact_min_distance_picks_the_nearest_row_not_the_first():
    far = "T" * SEQ_LENGTH
    near = "A" * (SEQ_LENGTH - 2) + "TT"
    d = hamming_min_to_set(encode(["A" * SEQ_LENGTH]), encode([far, near]))
    assert d[0] == 2


def test_identical_sequence_present_gives_zero():
    d = hamming_min_to_set(encode(["ACGT" * 10]), encode(["ACGT" * 10, "T" * SEQ_LENGTH]))
    assert d[0] == 0


@pytest.mark.parametrize("chunk", [1, 2, 3, 7, 64, 1000])
def test_chunking_never_changes_a_distance(chunk):
    """The whole point of the chunk parameter is bounded memory at identical output."""
    rng = np.random.default_rng(0)
    test = rng.integers(0, 4, size=(23, SEQ_LENGTH), dtype=np.uint8)
    train = rng.integers(0, 4, size=(17, SEQ_LENGTH), dtype=np.uint8)
    reference = hamming_min_to_set(test, train, chunk=10_000)
    assert np.array_equal(hamming_min_to_set(test, train, chunk=chunk), reference)


def test_chunk_smaller_than_input_still_covers_every_row():
    """A chunk loop that drops its tail would leave later rows at their uninitialised value."""
    rng = np.random.default_rng(1)
    test = rng.integers(0, 4, size=(10, SEQ_LENGTH), dtype=np.uint8)
    train = test.copy()  # every row has an exact match, so every answer must be 0
    assert hamming_min_to_set(test, train, chunk=3).tolist() == [0] * 10


def test_empty_train_set_raises_rather_than_returning_infinity():
    with pytest.raises(ValueError):
        hamming_min_to_set(encode(["A" * SEQ_LENGTH]), np.empty((0, SEQ_LENGTH), dtype=np.uint8))


# ------------------------------------------------------------------------ binning


def test_bins_are_exhaustive_over_every_reachable_distance():
    """bin_index's own error message claims the bins are meant to be exhaustive. Hold it to that."""
    for d in range(0, SEQ_LENGTH + 1):
        bin_index(d)


def test_bins_do_not_overlap():
    seen: dict[int, int] = {}
    for d in range(0, SEQ_LENGTH + 1):
        k = bin_index(d)
        seen.setdefault(k, 0)
        seen[k] += 1
    assert sum(seen.values()) == SEQ_LENGTH + 1
    assert len(seen) == len(BINS), "a declared bin captured no distance at all"


def test_bin_boundaries_land_where_declared():
    assert (bin_index(0), bin_index(1)) == (0, 0)
    assert (bin_index(2), bin_index(3)) == (1, 1)
    assert bin_index(7) == 2
    assert bin_index(8) == 3


def test_out_of_range_distance_raises_rather_than_being_dropped():
    with pytest.raises(ValueError):
        bin_index(SEQ_LENGTH + 1)
    with pytest.raises(ValueError):
        bin_index(-1)


# ------------------------------------------------------------------- fisher_z_ci


def test_interval_refused_below_the_fixed_n_floor():
    for n in range(0, MIN_N_FOR_CI):
        assert fisher_z_ci(0.5, n) is None


def test_interval_exists_at_the_floor_itself():
    assert fisher_z_ci(0.5, MIN_N_FOR_CI) is not None


def test_interval_brackets_the_point_estimate():
    lo, hi = fisher_z_ci(0.5, 100)
    assert lo < 0.5 < hi


def test_interval_widens_as_n_shrinks():
    wide = fisher_z_ci(0.5, 10)
    narrow = fisher_z_ci(0.5, 500)
    assert (wide[1] - wide[0]) > (narrow[1] - narrow[0])


def test_perfect_correlation_has_no_finite_interval():
    """Fisher z is undefined at +/-1; returning a number there would be fabricated precision."""
    assert fisher_z_ci(1.0, 100) is None
    assert fisher_z_ci(-1.0, 100) is None


def test_interval_stays_inside_the_correlation_range():
    lo, hi = fisher_z_ci(0.95, 5)
    assert -1.0 <= lo <= hi <= 1.0


# ---------------------------------------------------------------- summarise_bins


def _distances_for(counts: dict[int, int]) -> list[int]:
    """Build a distance list placing `n` points in each requested bin index."""
    out: list[int] = []
    for k, n in counts.items():
        lo, _hi = BINS[k]
        out.extend([lo] * n)
    return out


def test_bin_below_the_floor_refuses_and_reports_no_correlation():
    d = _distances_for({0: MIN_BIN_N - 1})
    rows = summarise_bins(d, list(range(len(d))), list(range(len(d))), spearman_fn=lambda a, b: 1.0)
    row = rows[0]
    assert row["status"] == "refused:below_floor"
    assert row["spearman"] is None and row["ci95"] is None
    assert str(MIN_BIN_N) in row["why"]


def test_the_floor_is_exact_at_the_boundary():
    """n == MIN_BIN_N scores; n == MIN_BIN_N - 1 refuses. An off-by-one here reports a refused bin."""
    at = _distances_for({0: MIN_BIN_N})
    below = _distances_for({0: MIN_BIN_N - 1})
    f = lambda a, b: 0.4  # noqa: E731
    assert summarise_bins(at, list(range(len(at))), list(range(len(at))), spearman_fn=f)[0]["status"] == "scored"
    assert summarise_bins(below, list(range(len(below))), list(range(len(below))), spearman_fn=f)[0]["status"] == "refused:below_floor"


def test_every_declared_bin_appears_in_the_output_even_when_empty():
    """A silently omitted empty bin would make a curve look denser than it is."""
    rows = summarise_bins([0] * MIN_BIN_N, list(range(MIN_BIN_N)), list(range(MIN_BIN_N)),
                          spearman_fn=lambda a, b: 0.1)
    assert len(rows) == len(BINS)
    assert [r["bin_index"] for r in rows] == list(range(len(BINS)))
    assert rows[1]["n"] == 0 and rows[1]["status"] == "refused:below_floor"


def test_only_the_selected_points_reach_the_metric():
    """Bin selection must subset both vectors identically, or predictions pair with wrong measurements."""
    seen: dict[str, list] = {}

    def spy(p, m):
        seen["p"], seen["m"] = list(p), list(m)
        return 0.0

    # MIN_BIN_N points in bin 0 (distance 0) then MIN_BIN_N in bin 2 (distance 4).
    d = [0] * MIN_BIN_N + [4] * MIN_BIN_N
    preds = [10.0] * MIN_BIN_N + [99.0] * MIN_BIN_N
    meas = [1.0] * MIN_BIN_N + [2.0] * MIN_BIN_N
    summarise_bins(d, preds, meas, spearman_fn=spy)
    # Bin 2 is summarised last, so the spy holds its slice.
    assert seen["p"] == [99.0] * MIN_BIN_N
    assert seen["m"] == [2.0] * MIN_BIN_N


# ----------------------------------------------------------- monotonicity_report


def test_fewer_than_two_scored_bins_is_indeterminate_not_a_trend():
    rows = [{"status": "scored", "bin_index": 0, "spearman": 0.5},
            {"status": "refused:below_floor", "bin_index": 1, "spearman": None}]
    rep = monotonicity_report(rows)
    assert rep["verdict"] == "indeterminate:fewer_than_two_scored_bins"
    assert rep["n_scored_bins"] == 1


def test_refused_bins_cannot_manufacture_a_monotone_verdict():
    """The specific artefact this harness exists to avoid: a 'decay' built from unmeasured bins."""
    rows = [{"status": "refused:below_floor", "bin_index": k, "spearman": None} for k in range(len(BINS))]
    rep = monotonicity_report(rows)
    assert rep["verdict"] == "indeterminate:fewer_than_two_scored_bins"
    assert rep["n_scored_bins"] == 0


def test_declining_scored_bins_are_reported_as_monotone():
    rows = [{"status": "scored", "bin_index": 0, "spearman": 0.6},
            {"status": "scored", "bin_index": 1, "spearman": 0.4},
            {"status": "scored", "bin_index": 2, "spearman": 0.1}]
    rep = monotonicity_report(rows)
    assert rep["verdict"] == "declines_monotonically"
    assert rep["successive_deltas"] == [-0.2, -0.3]
    assert "caveat" in rep, "a monotone verdict must ship with its own weakness stated"


def test_a_single_rise_breaks_monotonicity():
    rows = [{"status": "scored", "bin_index": 0, "spearman": 0.6},
            {"status": "scored", "bin_index": 1, "spearman": 0.7},
            {"status": "scored", "bin_index": 2, "spearman": 0.1}]
    assert monotonicity_report(rows)["verdict"] == "not_monotonic"


def test_a_flat_run_is_not_a_decline():
    """Strict decline is required; equal successive values must not count as decay."""
    rows = [{"status": "scored", "bin_index": 0, "spearman": 0.5},
            {"status": "scored", "bin_index": 1, "spearman": 0.5}]
    assert monotonicity_report(rows)["verdict"] == "not_monotonic"


def test_monotonicity_uses_only_scored_bins_and_keeps_their_indices():
    rows = [{"status": "scored", "bin_index": 0, "spearman": 0.8},
            {"status": "refused:below_floor", "bin_index": 1, "spearman": None},
            {"status": "scored", "bin_index": 2, "spearman": 0.3}]
    rep = monotonicity_report(rows)
    assert rep["scored_bin_indices"] == [0, 2]
    assert rep["spearman_by_scored_bin"] == [0.8, 0.3]


# ------------------------------------------------------- pre-registered constants


def test_the_floor_and_bins_are_fixed_constants_not_runtime_derived():
    """These are the pre-registration. If they become data-dependent the refusals mean nothing."""
    assert isinstance(MIN_BIN_N, int) and MIN_BIN_N > 1
    assert BINS[0][0] == 0, "the first bin must start at an exact match"
    assert BINS[-1][1] == SEQ_LENGTH, "the last bin must reach the maximum possible distance"
    assert all(lo <= hi for lo, hi in BINS)
    assert all(BINS[k][1] + 1 == BINS[k + 1][0] for k in range(len(BINS) - 1)), "bins must be contiguous"
