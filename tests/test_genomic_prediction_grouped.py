"""Grouped-CV seam in `eval/genomic_prediction.py` (F3 Step 2).

The load-bearing property is that the permutation null consumes the SAME `folds` object as the
measured run, so a grouped run has a grouped null by construction rather than via a second code
path. That is asserted DIRECTLY here (by capturing the argument) rather than inferred from a number.
"""
from __future__ import annotations

import numpy as np
import pytest

from dna_decode.eval import genomic_prediction as gp
from dna_decode.eval.genomic_prediction import (
    GroupedFoldError,
    _grouped_fold_indices,
    _inner_split,
    _resolve_folds,
    cv_model_gp,
    cv_ridge_gp,
)


def _fixture(n=90, m=40, seed=7, with_nan=True):
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 2, size=(n, m)).astype(float)
    beta = np.zeros(m)
    beta[:5] = [1.4, -1.1, 0.9, 0.7, -0.6]
    y = X @ beta + rng.normal(0, 1.0, n)
    if with_nan:
        y[3] = np.nan
    return X, y


# --------------------------------------------------------------------------- no-op proof

# Produced by the PRE-change implementation (git HEAD before the F3 Step 2 edit), via
# scripts-free direct execution of the original module on `_fixture()`. Hardcoded on purpose:
# recomputing it from the implementation under test would make the assertion vacuous.
GOLDEN_RIDGE_SEED0_R = 0.4723343055050292


def test_groups_none_reproduces_the_pre_change_value_exactly():
    X, y = _fixture()
    got = cv_ridge_gp(X, y, trait="t", seed=0, n_perm=25)
    assert got.predictive_r == GOLDEN_RIDGE_SEED0_R


def test_cv_model_gp_groups_none_reproduces_the_pre_change_PATH_exactly():
    """`cv_ridge_gp`'s no-op is pinned by a hardcoded golden value above, but `cv_model_gp` had
    its NaN filter and fold construction replaced too and has no such pin.

    Rather than a second golden constant, the pre-change path is RECONSTRUCTED here from the two
    primitives the edit did not touch (`_kfold_indices`, `_cv_r_with`). That is the old code,
    not a re-run of the new code, so the comparison is not circular.
    """
    from dna_decode.eval.genomic_prediction import _cv_r_with, _kfold_indices, _ridge_fit_predict

    X, y = _fixture()
    Xa, ya = np.asarray(X, float), np.asarray(y, float)
    keep = ~np.isnan(ya)
    Xa, ya = Xa[keep], ya[keep]
    folds = _kfold_indices(len(ya), 5, 0)
    want_r, want_r2 = _cv_r_with(Xa, ya, folds,
                                 lambda a, b, c: _ridge_fit_predict(a, b, c, 100.0))

    got = cv_model_gp(X, y, model="ridge", n_splits=5, seed=0, n_perm=0)
    assert got.predictive_r == want_r
    assert got.r2_score == want_r2
    assert got.n == int(keep.sum())
    assert got.fold_scheme == "random" and got.n_groups is None


def test_cv_model_gp_grouped_path_works_for_the_tree_model_too():
    """The grouped seam is model-agnostic by design (`_resolve_folds` never sees the model), so
    the non-ridge arm must report the scheme rather than silently staying random."""
    X, y = _fixture(n=60, m=12)
    res = cv_model_gp(X, y, model="gbm", n_splits=5, seed=0, n_perm=0,
                      groups=[f"g{i % 10}" for i in range(len(y))])
    assert res.fold_scheme == "grouped"
    assert res.n_groups == 10


def test_groups_none_reports_random_scheme_and_no_group_count():
    X, y = _fixture()
    res = cv_ridge_gp(X, y, seed=0, n_perm=5)
    assert res.fold_scheme == "random"
    assert res.n_groups is None
    assert res.as_dict()["fold_scheme"] == "random"


def test_groups_is_keyword_only():
    X, y = _fixture()
    with pytest.raises(TypeError):
        cv_ridge_gp(X, y, "trait")  # type: ignore[misc]


# --------------------------------------------------------------------------- fold construction

def test_grouped_folds_never_split_a_group_across_train_and_test():
    groups = [f"g{i % 11}" for i in range(88)]
    folds = _grouped_fold_indices(groups, 5, seed=0)
    assert sum(len(f) for f in folds) == 88
    assert len({int(i) for f in folds for i in f}) == 88          # a partition, no duplicates
    for f, test_idx in enumerate(folds):
        train_idx = np.concatenate([folds[j] for j in range(len(folds)) if j != f])
        test_groups = {groups[i] for i in test_idx}
        train_groups = {groups[i] for i in train_idx}
        assert not (test_groups & train_groups), f"fold {f} leaks groups {test_groups & train_groups}"


def test_grouped_folds_handle_uneven_group_sizes():
    groups = ["big"] * 20 + ["mid"] * 10 + ["mid2"] * 10 + [f"s{i}" for i in range(20)]
    folds = _grouped_fold_indices(groups, 5, seed=0)
    for f, test_idx in enumerate(folds):
        train_idx = np.concatenate([folds[j] for j in range(len(folds)) if j != f])
        assert not ({groups[i] for i in test_idx} & {groups[i] for i in train_idx})


def test_grouped_folds_are_deterministic_for_a_seed():
    groups = [f"g{i % 9}" for i in range(72)]
    a = _grouped_fold_indices(groups, 4, seed=3)
    b = _grouped_fold_indices(groups, 4, seed=3)
    for x, y_ in zip(a, b):
        assert np.array_equal(x, y_)


def test_grouped_folds_agree_with_GroupKFold_on_equal_sized_groups():
    """Cross-check against sklearn on a benign case.

    GroupKFold is NOT the implementation: it balances group SIZES and silently reduces the
    effective number of splits when one group dominates, where this seam must REFUSE instead
    (see `test_refuses_a_dominant_group`). Agreement on the easy case is a sanity check only.
    """
    from sklearn.model_selection import GroupKFold

    groups = np.array([f"g{i // 8}" for i in range(80)])        # 10 groups x 8 rows
    mine = _grouped_fold_indices(list(groups), 5, seed=0)
    theirs = [te for _, te in GroupKFold(n_splits=5).split(np.zeros((80, 2)), groups=groups)]
    mine_part = sorted(sorted(groups[i] for i in f) for f in mine)
    theirs_part = sorted(sorted(groups[i] for i in f) for f in theirs)
    # Same partition SHAPE (2 groups per fold); which groups land together may differ by seed.
    assert [len(p) for p in mine_part] == [len(p) for p in theirs_part]
    for f in mine:
        assert len({groups[i] for i in f}) == 2


# --------------------------------------------------------------------------- the null shares the folds

def test_the_null_is_computed_on_the_SAME_folds_object_ridge(monkeypatch):
    """The whole design leans on this: assert the identity of the argument, not a number."""
    X, y = _fixture()
    seen: list[int] = []
    real = gp._cv_predictive_r

    def spy(Xa, ya, alphas, folds, groups=None):
        seen.append(id(folds))
        return real(Xa, ya, alphas, folds, groups)

    monkeypatch.setattr(gp, "_cv_predictive_r", spy)
    cv_ridge_gp(X, y, seed=0, n_perm=4, groups=[f"g{i % 10}" for i in range(len(y))])
    assert len(seen) == 5                      # 1 measured + 4 null
    assert len(set(seen)) == 1, "null used a different folds object than the measured run"


def test_the_null_is_computed_on_the_SAME_folds_object_model(monkeypatch):
    X, y = _fixture()
    seen: list[int] = []
    real = gp._cv_r_with

    def spy(Xa, ya, folds, fit_predict):
        seen.append(id(folds))
        return real(Xa, ya, folds, fit_predict)

    monkeypatch.setattr(gp, "_cv_r_with", spy)
    cv_model_gp(X, y, model="ridge", seed=0, n_perm=3,
                groups=[f"g{i % 10}" for i in range(len(y))])
    assert len(seen) == 4
    assert len(set(seen)) == 1


def test_grouped_run_reports_the_scheme_and_group_count():
    X, y = _fixture()
    groups = [f"g{i % 12}" for i in range(len(y))]
    res = cv_ridge_gp(X, y, seed=0, n_perm=3, groups=groups)
    assert res.fold_scheme == "grouped"
    # 90 rows, one dropped for NaN y -> 89 rows, 12 distinct groups
    assert res.n_groups == 12
    assert res.as_dict()["n_groups"] == 12


# --------------------------------------------------------------------------- inner-split leak

def _LEAKY_GROUPS():
    """40 rows in 3-row groups. `_inner_split`'s n_in = 40//5 = 8 falls INSIDE group g2
    (rows 6,7,8), so a positional slice necessarily splits it."""
    return [f"g{i // 3}" for i in range(40)]


def test_inner_split_with_groups_keeps_whole_groups_on_one_side():
    """Regression for the alpha-pick leak: a positional slice can put two members of one group
    on opposite sides of the inner split, so the alpha choice leaks even with a clean OUTER fold.
    """
    groups = _LEAKY_GROUPS()
    train_idx = np.arange(40)
    va, tr2 = _inner_split(train_idx, groups, rng_seed=0)
    va_groups = {groups[i] for i in va}
    tr2_groups = {groups[i] for i in tr2}
    assert va_groups and tr2_groups
    assert not (va_groups & tr2_groups), "inner split leaks a group across the alpha pick"


def test_the_positional_inner_split_DOES_leak_so_the_fix_is_not_vacuous():
    """Pin that the old behaviour was genuinely broken ON THE SAME FIXTURE.

    Group size 3 against n_in=8 is deliberate: with 4-row groups the positional slice lands
    exactly on a group boundary and leaks nothing, so that fixture would make the fix above
    look proven while testing nothing.
    """
    groups = _LEAKY_GROUPS()
    train_idx = np.arange(40)
    n_in = max(1, len(train_idx) // 5)                # the pre-fix positional slice
    va, tr2 = train_idx[:n_in], train_idx[n_in:]
    assert {groups[i] for i in va} & {groups[i] for i in tr2}, "fixture does not exercise the leak"


def test_inner_split_without_groups_is_the_positional_slice():
    train_idx = np.arange(40)
    va, tr2 = _inner_split(train_idx, None, rng_seed=0)
    assert np.array_equal(va, train_idx[:8])
    assert np.array_equal(tr2, train_idx[8:])


def test_inner_split_falls_back_when_one_group_holds_everything():
    groups = ["only"] * 40
    va, tr2 = _inner_split(np.arange(40), groups, rng_seed=0)
    assert len(va) > 0 and len(tr2) > 0


# --------------------------------------------------------------------------- refusals

def test_refuses_length_mismatch_and_names_both_counts():
    X, y = _fixture(with_nan=False)
    with pytest.raises(GroupedFoldError) as e:
        cv_ridge_gp(X, y, seed=0, n_perm=0, groups=["a"] * 5)
    assert "5" in str(e.value) and str(len(y)) in str(e.value)


def test_refuses_fewer_groups_than_splits_and_names_the_count():
    X, y = _fixture(with_nan=False)
    groups = [f"g{i % 3}" for i in range(len(y))]
    with pytest.raises(GroupedFoldError) as e:
        cv_ridge_gp(X, y, n_splits=5, seed=0, n_perm=0, groups=groups)
    assert "3 distinct groups" in str(e.value)


def test_refuses_a_dominant_group_and_names_its_share():
    X, y = _fixture(with_nan=False)
    n = len(y)
    groups = ["dom"] * (n // 2) + [f"g{i}" for i in range(n - n // 2)]
    with pytest.raises(GroupedFoldError) as e:
        cv_ridge_gp(X, y, n_splits=5, seed=0, n_perm=0, groups=groups)
    msg = str(e.value)
    assert "dom" in msg and "leave-" in msg


def test_a_group_exactly_at_the_one_fold_share_is_accepted():
    """The refusal is strictly ABOVE n/n_splits -- a group that is exactly one fold is fine."""
    X, y = _fixture(n=100, with_nan=False)
    groups = ["twenty"] * 20 + [f"g{i % 9}" for i in range(80)]
    res = cv_ridge_gp(X, y, n_splits=5, seed=0, n_perm=0, groups=groups)
    assert res.fold_scheme == "grouped"


def test_filter_rows_applies_the_nan_mask_to_groups_too():
    """A y-NaN filter that does not filter `groups` misaligns labels against rows silently."""
    X, y = _fixture(n=20, m=5, with_nan=False)
    y[0] = np.nan
    y[7] = np.nan
    groups = [f"g{i}" for i in range(20)]
    Xf, yf, gf = gp._filter_rows(X, y, groups)
    assert len(yf) == 18 and len(gf) == 18
    assert "g0" not in gf and "g7" not in gf
    assert gf[0] == "g1"


def test_filter_rows_refuses_a_groups_length_mismatch_before_masking():
    X, y = _fixture(n=20, m=5, with_nan=False)
    with pytest.raises(GroupedFoldError):
        gp._filter_rows(X, y, ["a", "b"])


def test_resolve_folds_passthrough_for_none():
    y = np.zeros(50)
    folds, scheme, ng = _resolve_folds(y, None, 5, 0)
    assert scheme == "random" and ng is None and len(folds) == 5
