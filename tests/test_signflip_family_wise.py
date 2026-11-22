"""Pin RULE_2's family-wise sign-flip permutation.

The property that justifies the whole construction is that flipping WHOLE ROWS carries the
cross-candidate correlation into the null instead of assuming it away. Two consequences are
testable and both are pinned here:

  * on a PERFECTLY CORRELATED family the family-wise p equals the per-candidate p exactly -- the
    correction costs nothing, because there was only ever one distinct candidate. Bonferroni, which
    assumes independence, multiplies by k and throws that power away.
  * on an INDEPENDENT family the family-wise p rises well above the smallest per-candidate p -- the
    correction bites, which is the under-correction the rule exists to prevent.

Also pinned: p_family >= p of the candidate attaining the family statistic. That is a theorem, not
an empirical observation -- the event {max_j |t_j| >= M} contains {|t_j*| >= M} for the j* attaining
M -- so a violation means the implementation is wrong rather than the data unusual.

No torch, no substrate, no network.
"""
from __future__ import annotations

import math
import random

import pytest

from dna_decode.eval.signflip import (
    EXACT_MAX_ROWS,
    SignFlipRefusal,
    SignFlipResult,
    bonferroni_p,
    family_wise_signflip,
)


# --------------------------------------------------------------------- refusals

def test_one_replicate_row_is_refused_as_a_unit_of_replication_problem():
    with pytest.raises(SignFlipRefusal) as e:
        family_wise_signflip([[0.1, 0.2]])
    msg = str(e.value)
    assert "RULE_1" in msg or "UNIT-OF-REPLICATION" in msg, (
        "the refusal must name the unit-of-replication question; a caller who hits this is "
        "usually about to pass re-scorings of one generated set as replicates")


def test_a_one_candidate_family_is_refused():
    """A family-wise correction over a single candidate computes a per-candidate p and invites it
    to be reported as though a selection had been priced."""
    with pytest.raises(SignFlipRefusal):
        family_wise_signflip([[0.1], [0.2], [0.3]])


def test_ragged_rows_are_refused():
    with pytest.raises(SignFlipRefusal):
        family_wise_signflip([[0.1, 0.2], [0.3]])


def test_a_non_finite_cell_is_refused_rather_than_propagated():
    with pytest.raises(SignFlipRefusal):
        family_wise_signflip([[0.1, float("nan")], [0.2, 0.3], [0.1, 0.1]])


def test_a_zero_variance_candidate_is_refused_under_t_and_allowed_under_mean():
    rows = [[0.05, 0.1], [0.05, 0.3], [0.05, -0.2]]
    with pytest.raises(SignFlipRefusal) as e:
        family_wise_signflip(rows, statistic="max_abs_t")
    assert "statistic='mean'" in str(e.value), "the refusal must name the usable alternative"
    r = family_wise_signflip(rows, statistic="mean")
    assert r.statistic == "mean"


def test_raw_scores_are_FLAGGED_not_refused_and_the_reason_is_undecidability():
    """Sign-flipping asks whether the values are symmetric about ZERO, so raw scores near 0.5 would
    return a p that answers nothing -- and the function says so in a note.

    It must NOT refuse. A first version did, and it rejected the single most valuable legitimate
    input: a strong consistent effect with small variance is numerically indistinguishable from
    raw scores. Warning beats overruling when the two cannot be separated.
    """
    raw = [[0.50, 0.52], [0.49, 0.53], [0.51, 0.50], [0.50, 0.51]]
    r = family_wise_signflip(raw, statistic="mean")
    assert any("UNCENTRED" in n for n in r.notes)
    assert any("cannot tell the two apart" in n for n in r.notes), (
        "the note must concede it is a heuristic rather than implying a verdict")
    strong = [[0.30, 0.31], [0.29, 0.30], [0.31, 0.29], [0.30, 0.30]]
    assert family_wise_signflip(strong, statistic="mean").family_wise_p > 0
    quiet = family_wise_signflip(raw, statistic="mean", flag_uncentred=False)
    assert not any("UNCENTRED" in n for n in quiet.notes)


def test_an_unknown_statistic_is_refused():
    with pytest.raises(SignFlipRefusal):
        family_wise_signflip([[0.1, 0.2], [0.3, 0.1]], statistic="median")


# --------------------------------------------------------------------- the correlation property

def _identical_columns(n=10, k=4, effect=0.04, noise=0.02, seed=1):
    g = random.Random(seed)
    rows = []
    for _ in range(n):
        v = effect + g.gauss(0, noise)
        rows.append([v] * k)          # perfect correlation: one candidate, repeated
    return rows


def test_a_perfectly_correlated_family_costs_nothing_to_correct():
    """THE load-bearing property. Four identical candidates are one candidate; flipping whole rows
    keeps them identical in the null, so the family-wise p equals the per-candidate p EXACTLY."""
    r = family_wise_signflip(_identical_columns())
    assert r.exact, "10 rows should enumerate exactly, leaving no sampling error to hide behind"
    assert r.family_wise_p == pytest.approx(r.per_candidate_p[0], abs=1e-12)
    assert all(p == pytest.approx(r.per_candidate_p[0], abs=1e-12) for p in r.per_candidate_p)


def test_bonferroni_overcorrects_the_same_correlated_family():
    """Bonferroni assumes independence, so on this family it multiplies a correct p by 4 for
    nothing. This is what the sign-flip test replaces."""
    r = family_wise_signflip(_identical_columns())
    bon = bonferroni_p(r.per_candidate_p)
    assert bon > r.family_wise_p, "with k=4 identical candidates Bonferroni must be strictly worse"
    assert bon == pytest.approx(min(1.0, 4 * r.family_wise_p))


def test_an_independent_family_does_get_corrected():
    """The other direction: with genuinely independent candidates the family-wise p rises well
    above the smallest per-candidate p. A correction that never bites would be useless."""
    g = random.Random(7)
    rows = [[g.gauss(0, 1) for _ in range(6)] for _ in range(12)]
    r = family_wise_signflip(rows, statistic="mean")
    assert r.family_wise_p > min(r.per_candidate_p) * 1.5, (
        "selecting the best of 6 independent candidates must cost something")


@pytest.mark.parametrize("builder,seed", [("identical", 1), ("identical", 2), ("independent", 3)])
def test_family_p_is_never_below_the_p_of_the_candidate_attaining_it(builder, seed):
    """A theorem, so a violation is an implementation bug and not unusual data."""
    if builder == "identical":
        rows = _identical_columns(seed=seed)
        r = family_wise_signflip(rows)
    else:
        g = random.Random(seed)
        rows = [[g.gauss(0, 1) for _ in range(5)] for _ in range(11)]
        r = family_wise_signflip(rows, statistic="mean")
    jstar = max(range(r.n_candidates), key=lambda j: abs(r.observed[j]))
    assert r.family_wise_p >= r.per_candidate_p[jstar] - 1e-12


def test_flipping_whole_rows_is_not_the_same_as_flipping_cells():
    """Pin that the row structure is actually used. Flipping cells independently destroys the
    within-row correlation; on a correlated family the two nulls differ, so a future refactor to
    per-cell flipping cannot pass silently."""
    rows = _identical_columns(n=9, k=4)
    r = family_wise_signflip(rows)

    def cellwise_family_p(mat):
        import itertools
        n, k = len(mat), len(mat[0])

        def stats(m):
            out = []
            for j in range(k):
                col = [row[j] for row in m]
                mean = sum(col) / n
                var = sum((v - mean) ** 2 for v in col) / (n - 1)
                out.append(mean / math.sqrt(var / n) if var > 0 else 0.0)
            return out

        obs = max(abs(v) for v in stats(mat))
        g = random.Random(0)
        ge = 0
        draws = 4000
        for _ in range(draws):
            flipped = [[(1 if g.random() < 0.5 else -1) * v for v in row] for row in mat]
            if max(abs(v) for v in stats(flipped)) >= obs - 1e-12:
                ge += 1
        return ge / draws

    # Direction, not mere difference. With identical columns, cell-wise flipping makes the four
    # columns near-independent in the null, so the null MAXIMUM is stochastically larger and the p
    # strictly more conservative. A tolerance check at abs=0.02 passed by accident because the gap
    # is only ~0.012, which is why this asserts the inequality instead.
    cw = cellwise_family_p(rows)
    assert cw > r.family_wise_p, (
        f"cell-wise p {cw} should exceed row-wise {r.family_wise_p}: destroying the within-row "
        f"correlation inflates the null maximum. Coinciding means the row structure is unused.")


# --------------------------------------------------------------------- behaviour + honesty

def test_a_strong_consistent_effect_reaches_the_attainable_floor():
    rows = [[0.10, 0.11, 0.09, 0.12] for _ in range(8)]
    r = family_wise_signflip(rows, statistic="mean")
    # TWO-SIDED: |mean| is unchanged by flipping EVERY row, so all-plus AND all-minus both attain
    # the observed statistic. The floor is 2/2**n. A first version asserted 1/2**n -- a p this
    # test can never produce.
    assert r.family_wise_p == pytest.approx(2 / 2 ** 8, abs=1e-12), (
        "exactly the all-plus and all-minus vectors attain the statistic")
    assert r.family_wise_p > 0


def test_a_null_family_is_not_significant():
    g = random.Random(11)
    rows = [[g.gauss(0, 1) for _ in range(4)] for _ in range(12)]
    r = family_wise_signflip(rows, statistic="mean")
    assert r.family_wise_p > 0.05


def test_p_is_never_zero_in_either_enumeration_mode():
    """The observed configuration is included deliberately, so a p of 0 -- which would read as
    impossible rather than as unresolved -- cannot occur."""
    small = family_wise_signflip([[0.3, 0.31], [0.29, 0.30], [0.31, 0.29]], statistic="mean")
    assert small.exact and small.family_wise_p > 0
    g = random.Random(5)
    big = [[0.2 + g.gauss(0, 0.01) for _ in range(3)] for _ in range(EXACT_MAX_ROWS + 4)]
    mc = family_wise_signflip(big, statistic="mean", n_draws=2000)
    assert not mc.exact and mc.family_wise_p > 0
    assert mc.family_wise_p >= 1 / mc.n_draws


def test_small_row_counts_are_flagged_as_underpowered_not_as_absence():
    r = family_wise_signflip([[0.1, 0.2], [0.15, 0.1], [0.2, 0.25]], statistic="mean")
    joined = " ".join(r.notes)
    assert "UNDERPOWERED" in joined and "not evidence of absence" in joined


def test_exact_mode_is_deterministic_and_ignores_the_seed():
    a = family_wise_signflip(_identical_columns(), seed=0)
    b = family_wise_signflip(_identical_columns(), seed=999)
    assert a.exact and a.family_wise_p == b.family_wise_p
    assert a.seed is None, "an exact enumeration has no seed to report; claiming one is misleading"


def test_monte_carlo_is_reproducible_for_a_given_seed():
    g = random.Random(3)
    rows = [[g.gauss(0.05, 1) for _ in range(3)] for _ in range(EXACT_MAX_ROWS + 6)]
    a = family_wise_signflip(rows, statistic="mean", n_draws=1500, seed=42)
    b = family_wise_signflip(rows, statistic="mean", n_draws=1500, seed=42)
    assert a.family_wise_p == b.family_wise_p and a.seed == 42


def test_the_record_states_what_it_does_not_cover():
    r = family_wise_signflip(_identical_columns())
    d = r.as_dict()
    assert "per_candidate_p_UNCORRECTED" in d, (
        "the uncorrected per-candidate p must be labelled as such wherever it is published")
    assert "unit of replication" in d["what_it_does_not_cover"]
    assert d["enumeration"] == "exact"
    assert isinstance(r, SignFlipResult)
