"""Guards for the G-C conditioning arm.

The thing most worth pinning here is not a number -- it is the REASON the primary arm is `partial_pooling`
and not `interaction`. That was a measured correction mid-build: the first design made the primary arm a
medium INDICATOR, which can only learn a global per-medium offset and is strictly LESS expressive than the
two-model comparator; correcting it to a FULL interaction over-shot to exactly EQUIVALENT (a full
interaction with a binary medium spans the same hypothesis space as two independent per-medium models), so
neither arm could ever show a gain. `partial_pooling` sits strictly between the two endpoints and is the
only one of the three that can beat both.

A future edit that "simplifies" the scale back to 1.0, or re-points the verdict at `interaction`, would
silently restore a rigged design that can only return NO_GAIN. These tests fail loudly on both.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path("scripts/glm_condition_conditioning.py")


def _mod():
    spec = importlib.util.spec_from_file_location("cc", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _arms(two_model, partial_pooling, shuffled, data_matched, *, interaction=None, indicator=0.10):
    """Synthetic arm aggregate in the shape `verdict` consumes."""
    k = "per_medium_mean"
    return {
        "two_model": {k: two_model},
        "partial_pooling": {k: partial_pooling},
        "interaction": {k: two_model if interaction is None else interaction},
        "partial_pooling_shuffled": {k: shuffled},
        "data_matched": {k: data_matched},
        "indicator": {k: indicator},
    }


# ---------------------------------------------------------------------------------------------------
# THE MECHANISM: why partial_pooling exists at all
# ---------------------------------------------------------------------------------------------------
def test_full_interaction_IS_the_two_model_comparator_re_parameterised():
    """The measured degeneracy, pinned on synthetic data whose truth genuinely HAS medium-specific
    coefficients -- the regime where an interaction arm is supposed to win if it ever can.

    If this ever stops holding, the recorded justification for `partial_pooling` is wrong and the arm
    choice has to be re-derived rather than inherited.
    """
    import numpy as np
    cc = _mod()
    rng = np.random.default_rng(0)
    n, w = 400, 8
    S = rng.random((n, w))
    w_lb, w_m9 = rng.normal(size=w), rng.normal(size=w)  # genuinely different per medium

    y_lb = S @ w_lb
    y_m9 = S @ w_m9
    # two independent per-medium models
    p_lb = cc.fit_ridge(cc.build_design(S, [0.0] * n, arm="sequence_only"), y_lb,
                        cc.build_design(S, [0.0] * n, arm="sequence_only"))
    p_m9 = cc.fit_ridge(cc.build_design(S, [0.0] * n, arm="sequence_only"), y_m9,
                        cc.build_design(S, [0.0] * n, arm="sequence_only"))
    two_model = np.asarray(p_lb + p_m9)

    def conditioned(scale):
        X = np.vstack([cc.build_design(S, [0.0] * n, arm="interaction"),
                       cc.build_design(S, [1.0] * n, arm="interaction")])
        y = np.concatenate([y_lb, y_m9])
        cs = None if scale == 1.0 else cc.interaction_col_scale(X.shape[1], scale)
        return np.asarray(cc.fit_ridge(X, y, X, col_scale=cs))

    full = conditioned(1.0)
    partial = conditioned(0.1)  # a mid-grid pooling strength

    r_full = float(np.corrcoef(two_model, full)[0, 1])
    r_partial = float(np.corrcoef(two_model, partial)[0, 1])
    assert r_full > 0.999, f"full interaction should be ~identical to two_model, got r={r_full}"
    assert r_partial < r_full, (
        "partial pooling must be a GENUINELY different model from the comparator, otherwise the primary "
        f"arm is degenerate too (r_partial={r_partial}, r_full={r_full})")


def test_the_pooling_grid_spans_both_endpoints_without_degenerating_onto_zero():
    """1.0 must be IN the grid -- if no pooling genuinely is best, the selection has to be able to say so
    rather than being forced to claim pooling helps. 0.0 must NOT be, because it would kill the interaction
    block outright and silently re-enter the (already-failing) indicator arm under another name."""
    cc = _mod()
    assert 1.0 in cc.POOLING_GRID
    assert 0.0 not in cc.POOLING_GRID
    assert all(0.0 < c <= 1.0 for c in cc.POOLING_GRID)
    assert len(set(cc.POOLING_GRID)) == len(cc.POOLING_GRID)


def test_the_asserted_scale_is_recorded_as_SUPERSEDED_and_not_wired_to_anything():
    """The old constant stays visible in the record, but a grep for it must not find it driving a fit --
    otherwise the superseded value is still silently in the pipeline."""
    cc = _mod()
    assert cc.SUPERSEDED_ASSERTED_SCALE == 0.3
    src = SCRIPT.read_text(encoding="utf-8")
    assert "PARTIAL_POOLING_SCALE" not in src, "the superseded constant is still referenced somewhere"
    uses = [ln for ln in src.splitlines()
            if "SUPERSEDED_ASSERTED_SCALE" in ln and "interaction_scale" in ln]
    assert uses == [], f"superseded value is still feeding a fit: {uses}"


def test_the_pooling_scale_hits_ONLY_the_interaction_block():
    """The shared sequence block and the medium column must stay at 1.0, or partial pooling shrinks the
    shared model too and stops being a statement about condition-specificity."""
    cc = _mod()
    cs = cc.interaction_col_scale(129, 0.3)
    assert list(cs[:65]) == [1.0] * 65
    assert list(cs[65:]) == [0.3] * 64


def test_scaling_BEFORE_standardisation_would_be_a_NO_OP_which_is_why_it_is_applied_after():
    """The exact defect this design had, pinned so it cannot come back: `StandardScaler` divides out a
    constant column factor, so a pre-scaler multiply changes nothing. The post-scaler `col_scale` must
    change predictions; a pre-scaler one must not."""
    import numpy as np
    cc = _mod()
    rng = np.random.default_rng(3)
    X = rng.random((200, 9)).astype(np.float32)
    y = rng.random(200)

    pre = X.copy()
    pre[:, 5:] *= 0.3  # the broken placement
    assert cc.fit_ridge(pre, y, pre) == pytest.approx(cc.fit_ridge(X, y, X), abs=1e-6), \
        "pre-standardisation scaling must be inert -- that is the whole bug"

    post = cc.fit_ridge(X, y, X, col_scale=cc.interaction_col_scale(9, 0.3))
    assert post != pytest.approx(cc.fit_ridge(X, y, X), abs=1e-6), \
        "post-standardisation scaling must actually change the fit"


def test_col_scale_length_mismatch_raises_rather_than_broadcasting():
    """A wrong-length scale that silently broadcast would penalise the wrong columns and still return a
    plausible number."""
    import numpy as np
    cc = _mod()
    X = np.random.default_rng(1).random((40, 9)).astype(np.float32)
    with pytest.raises(ValueError):
        cc.fit_ridge(X, np.random.default_rng(2).random(40), X, col_scale=np.ones(5))


def test_a_scale_of_one_reproduces_the_unscaled_fit_exactly():
    """Non-vacuity anchor for the endpoint claim: c=1 must BE the no-pooling comparator, not merely near it."""
    import numpy as np
    cc = _mod()
    X = np.random.default_rng(4).random((120, 9)).astype(np.float32)
    y = np.random.default_rng(5).random(120)
    assert cc.fit_ridge(X, y, X, col_scale=np.ones(9)) == pytest.approx(cc.fit_ridge(X, y, X), abs=1e-9)


def test_the_interaction_block_is_inert_for_LB_by_construction():
    """medium=0 zeroes the interaction block, which is WHY the arm can give a feature a different
    coefficient per medium (LB uses w, M9 uses w+v) rather than merely an offset."""
    import numpy as np
    cc = _mod()
    S = np.asarray([[0.25, 0.5, 0.25]])
    lb = cc.build_design(S, [0.0], arm="interaction")
    assert np.all(np.asarray(lb[0, 4:]) == 0.0)
    m9 = cc.build_design(S, [1.0], arm="interaction")
    assert np.any(np.asarray(m9[0, 4:]) != 0.0)


def test_indicator_can_only_add_an_offset_column():
    """The arm retained to DEMONSTRATE failure. One extra column, not a per-feature interaction."""
    import numpy as np
    cc = _mod()
    S = np.asarray([[0.1] * 64])
    assert cc.build_design(S, [1.0], arm="indicator").shape[1] == 65
    assert cc.build_design(S, [1.0], arm="interaction").shape[1] == 129


def test_an_unknown_arm_raises_rather_than_silently_returning_sequence_only():
    """A typo'd arm name must not quietly score as the sequence-only baseline and look like a null."""
    import numpy as np
    cc = _mod()
    with pytest.raises(ValueError):
        cc.build_design(np.asarray([[0.5, 0.5]]), [1.0], arm="interation")


# ---------------------------------------------------------------------------------------------------
# the selection must never see the test split
# ---------------------------------------------------------------------------------------------------
def _frags(n=600):
    """Synthetic fragments laid out along a genome, so position blocking has something to block on."""
    import random
    from dna_decode.glm.genomewide import Fragment
    r = random.Random(0)
    out = []
    for i in range(n):
        start = i * 100
        out.append(Fragment(seq="".join(r.choice("ACGT") for _ in range(150)),
                            expression=r.random(), rep1=r.random(), rep2=r.random(),
                            start=start, end=start + 150, strand="+"))
    return out


def test_the_inner_split_is_a_PARTITION_of_the_training_indices_only():
    """The whole point of selecting on an inner split: if a single test index leaked into the selection,
    the reported gain would be partly selected ON the data it is scored against."""
    cc = _mod()
    frag = {"LB": _frags(), "M9": _frags()}
    idx_tr = list(range(0, 400))
    idx_te = list(range(400, 600))
    inner_tr, inner_va, _held = cc._inner_split(frag, idx_tr, seed=0, n_blocks=10, held_blocks=2)

    assert set(inner_tr) | set(inner_va) == set(idx_tr), "inner split must cover the training set exactly"
    assert set(inner_tr) & set(inner_va) == set(), "inner train and inner val must be disjoint"
    assert set(inner_tr).isdisjoint(idx_te), "TEST INDEX LEAKED into inner training"
    assert set(inner_va).isdisjoint(idx_te), "TEST INDEX LEAKED into inner validation"
    assert inner_tr and inner_va, "a degenerate inner split would make selection vacuous"


def test_the_inner_split_holds_out_CONTIGUOUS_blocks_not_scattered_rows():
    """Blocked for the same measured reason the outer split is: 91.5% of consecutive fragments overlap, so
    a scattered inner holdout would be validating against near-duplicates of its own training rows."""
    cc = _mod()
    frag = {"LB": _frags(), "M9": _frags()}
    _itr, iva, _h = cc._inner_split(frag, list(range(400)), seed=0, n_blocks=10, held_blocks=2)
    runs = 1 + sum(1 for a, b in zip(sorted(iva), sorted(iva)[1:]) if b != a + 1)
    assert runs <= 2 + 1, f"inner validation fell into {runs} fragments; expected contiguous blocks"


# ---------------------------------------------------------------------------------------------------
# the verdict rule: VALIDITY FIRST (decision D1)
# ---------------------------------------------------------------------------------------------------
def test_a_dirty_null_blocks_grading_even_with_a_spectacular_gain():
    """D1's whole point: a run that fabricates signal must not be graded, however good the headline looks.
    Give it a huge real gain AND a shuffled arm that also gains -- it must still refuse."""
    cc = _mod()
    v, d = cc.verdict(_arms(two_model=0.20, partial_pooling=0.80, shuffled=0.50, data_matched=0.80))
    assert v == "INDETERMINATE_NULL_NOT_CLEAN"
    assert "SHUFFLED" in d


def test_a_clean_null_then_grades():
    """Non-vacuity for the test above -- the validity gate must not block legitimate runs."""
    cc = _mod()
    v, _d = cc.verdict(_arms(two_model=0.20, partial_pooling=0.30, shuffled=0.20, data_matched=0.30))
    assert v == "CONDITIONING_HELPS"


def test_no_gain_is_the_verdict_when_the_margin_is_missed():
    cc = _mod()
    v, d = cc.verdict(_arms(two_model=0.2013, partial_pooling=0.2050, shuffled=0.2010,
                            data_matched=0.2050))
    assert v == "NO_GAIN", d


def test_a_gain_that_vanishes_under_row_matching_is_called_DATA_VOLUME_not_conditioning():
    """The shared model sees ~2x the rows. A gain that does not survive matching the row count is
    attributable to data volume, and saying so is the difference between a finding and an artifact."""
    cc = _mod()
    v, _d = cc.verdict(_arms(two_model=0.20, partial_pooling=0.25, shuffled=0.20, data_matched=0.205))
    assert v == "GAIN_IS_DATA_VOLUME"


def test_the_verdict_is_driven_by_partial_pooling_and_NOT_by_the_degenerate_interaction_arm():
    """Pins the redesign. `interaction` is reported for its demonstrative value only; if the verdict ever
    reads it as primary, the grade comes from an arm that cannot beat the comparator by construction."""
    cc = _mod()
    assert cc.PREREGISTERED["primary_arm"] == "partial_pooling"
    # a spectacular full-interaction number must change NOTHING
    base = _arms(two_model=0.20, partial_pooling=0.205, shuffled=0.20, data_matched=0.205)
    rigged = _arms(two_model=0.20, partial_pooling=0.205, shuffled=0.20, data_matched=0.205,
                   interaction=0.90)
    assert cc.verdict(base)[0] == cc.verdict(rigged)[0] == "NO_GAIN"


def test_the_shuffled_control_is_the_PARTIAL_POOLING_one_not_the_full_interaction_one():
    """The null has to be the null OF THE GRADED ARM. A shuffled control built on a different arm would
    license the primary number while testing something else."""
    cc = _mod()
    arms = _arms(two_model=0.20, partial_pooling=0.30, shuffled=0.20, data_matched=0.30)
    del arms["partial_pooling_shuffled"]
    with pytest.raises(KeyError):
        cc.verdict(arms)


# ---------------------------------------------------------------------------------------------------
# substrate integrity
# ---------------------------------------------------------------------------------------------------
def test_the_shared_fragment_floor_is_a_REFUSAL_not_a_warning():
    """A degraded intersection must not be scored. Pinned by reading the threshold rather than by
    downloading 90 MB: the floor exists and is high enough to matter."""
    cc = _mod()
    assert cc.PREREGISTERED["min_shared_fragments"] >= 100_000


def test_per_medium_mean_is_primary_and_pooled_is_secondary():
    """Measured: pooled 0.310 vs per-medium mean 0.197 for the same baseline. The gap is the between-media
    LEVEL difference, which a global offset captures without any condition-specific sequence response -- so
    grading on pooled would reward exactly the thing this family is not about."""
    cc = _mod()
    assert cc.PREREGISTERED["primary_metric"] == "per_medium_mean"
    assert cc.PREREGISTERED["secondary_metric"] == "pooled"


def test_both_metrics_are_reported_together_so_an_offset_effect_is_visible():
    cc = _mod()
    r = cc.per_medium_and_pooled({"LB": [1.0, 2.0, 3.0], "M9": [1.0, 2.0, 3.0]},
                                 {"LB": [1.0, 2.0, 3.0], "M9": [3.0, 2.0, 1.0]})
    assert set(r) == {"per_medium", "per_medium_mean", "pooled"}
    assert r["per_medium"]["LB"] == pytest.approx(1.0)
    assert r["per_medium"]["M9"] == pytest.approx(-1.0)
    assert r["per_medium_mean"] == pytest.approx(0.0)
