"""Guards for the held-out-group transfer harness (F1 Step 3-4).

Offline, pure. The only heavy import is sklearn inside the `ridge` helper, mirroring how
`eval/genomic_prediction.py` keeps it off the module import path.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.transfer import (  # noqa: E402
    MIN_SCORED, POWERED, UNDERPOWERED, WITHIN_GROUP_MIN_N, KShotGridError, UnassignedMemberError,
    derive_k_grid, fold_report, held_out_group_folds, k_shot_curve, power_check, skipped_groups,
)


def ridge(Xtr, ytr, Xte):
    from sklearn.linear_model import Ridge
    return Ridge(alpha=10.0).fit(Xtr, ytr).predict(Xte)


def constant_mean(Xtr, ytr, Xte):
    """A deliberately different method: predicts the training mean, ignoring X."""
    return np.full(len(Xte), float(np.mean(ytr)))


def _panel(n_groups=4, per_group=80, seed=0):
    """The `tests/test_deconfound_package.py:23` fixture shape: big per-group offsets plus a
    within-group signal on feature 0."""
    rng = np.random.default_rng(seed)
    n = n_groups * per_group
    g = np.repeat(np.arange(n_groups), per_group)
    X = rng.integers(0, 2, (n, 3)).astype(float)
    y = np.repeat(rng.normal(0, 5, n_groups), per_group) + 1.5 * X[:, 0] + rng.normal(0, .3, n)
    ids = [f"id{i}" for i in range(n)]
    return X, y, ids, {ids[i]: f"g{g[i]}" for i in range(n)}


# --- import posture -------------------------------------------------------------------------------

def test_importing_the_harness_pulls_no_scipy_or_sklearn():
    """`deconfound` imports both at MODULE scope; this module must import it lazily or every test
    collection touching eval/ pays for it."""
    import subprocess
    code = ("import sys, dna_decode.eval.transfer; "
            "print([m for m in sys.modules if m.startswith(('sklearn','scipy'))])")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(ROOT))
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[]", f"eager heavy import: {out.stdout.strip()}"


# --- folds ----------------------------------------------------------------------------------------

def test_folds_are_disjoint_by_construction_for_every_group():
    _, _, _, group_of = _panel()
    folds = held_out_group_folds(group_of)
    assert len(folds) == 4
    for f in folds:
        assert not (set(f.train_ids) & set(f.test_ids))
        assert set(f.train_ids) | set(f.test_ids) == set(group_of)
        assert all(group_of[i] == f.held_out_group for i in f.test_ids)


def test_an_unassigned_id_raises_by_default_and_buckets_only_when_asked():
    group_of = {"a": "g0", "b": "g0", "c": None, "d": ""}
    with pytest.raises(UnassignedMemberError):
        held_out_group_folds(group_of)
    folds = held_out_group_folds(group_of, allow_unassigned=True)
    assert "__unassigned__" in {f.held_out_group for f in folds}


def test_fold_report_surfaces_STRUCTURE_not_just_a_count():
    """43 clusters once hid 97% of the isolates. A count cannot show domination; a fraction can."""
    group_of = {f"id{i}": ("big" if i < 80 else f"s{i}") for i in range(100)}
    rep = fold_report(held_out_group_folds(group_of))
    assert rep["n_folds"] == 21
    assert rep["largest_group_fraction"] == 0.8, rep
    assert rep["singleton_fraction"] > 0.9, rep


# --- the silent-drop predictor ---------------------------------------------------------------------

def test_skipped_groups_names_the_groups_cv_r2_would_silently_drop():
    """SCOPE CORRECTION, recorded deliberately. The plan said this would be asserted against
    `cv_r2`'s observed per-group NaNs. That is NOT possible through its public API -- `cv_r2` returns
    a single float and its internal `pred` array never escapes. What IS externally observable is the
    BOUNDARY: when every group is skipped it produces no predictions at all and returns `nan`; when
    none are, it returns a finite number. So the prediction is checked against real behaviour at the
    boundary rather than against the docstring, which is weaker than the plan claimed and is the
    honest available check."""
    rng = np.random.default_rng(0)
    from dna_decode.deconfound import cv_r2

    # every group's train side is under 5 -> all skipped
    group_of = {f"id{i}": ("a" if i < 3 else "b") for i in range(6)}
    folds = held_out_group_folds(group_of)
    assert set(skipped_groups(folds)) == {"a", "b"}
    g = np.array(["a"] * 3 + ["b"] * 3)
    assert np.isnan(cv_r2(rng.normal(size=(6, 2)), rng.normal(size=6), groups=g)), \
        "cv_r2 should produce no predictions when every group is skipped"

    # none skipped -> finite
    group_of2 = {f"id{i}": f"g{i // 20}" for i in range(60)}
    assert skipped_groups(held_out_group_folds(group_of2)) == ()
    g2 = np.repeat(["g0", "g1", "g2"], 20)
    assert np.isfinite(cv_r2(rng.normal(size=(60, 2)), rng.normal(size=60), groups=g2))


# --- the derived grid and its two floors ------------------------------------------------------------

def test_the_grid_is_derived_from_the_BINDING_floor_which_is_30_not_10():
    assert WITHIN_GROUP_MIN_N > MIN_SCORED, "the point of max() is that 30 binds, not 10"
    assert derive_k_grid(80) == (0, 1, 2, 4, 8, 16, 32)
    assert derive_k_grid(35) == (0, 1, 2, 4)
    with pytest.raises(KShotGridError):
        derive_k_grid(30)
    # using MIN_SCORED alone would have produced a grid here -- that is the bug being prevented
    assert derive_k_grid(30, within_group_min_n=MIN_SCORED) == (0, 1, 2, 4, 8, 16)


def test_the_30_boundary_matches_within_group_r2s_MEASURED_behaviour():
    """Not asserted against the literal 30: asserted against what `within_group_r2` actually does on
    either side of it. At 20/group it returns (nan, 0 groups used) while the pooled metric returns a
    confident-looking number on the identical data -- which is why 30 binds."""
    from dna_decode.deconfound import cv_r2, within_group_r2

    X, y, ids, group_of = _panel(n_groups=6, per_group=20)
    g = np.array([group_of[i] for i in ids])
    r, used = within_group_r2(X, y, g, min_n=WITHIN_GROUP_MIN_N)
    assert np.isnan(r) and used == 0, (r, used)
    assert np.isfinite(cv_r2(X, y, groups=g)), "the pooled cell still reports a number -- the trap"

    X2, y2, ids2, group_of2 = _panel(n_groups=6, per_group=40)
    g2 = np.array([group_of2[i] for i in ids2])
    r2v, used2 = within_group_r2(X2, y2, g2, min_n=WITHIN_GROUP_MIN_N)
    assert np.isfinite(r2v) and used2 == 6, (r2v, used2)


# --- k-shot bookkeeping ----------------------------------------------------------------------------

def test_shots_are_never_scored_at_any_k():
    X, y, ids, group_of = _panel()
    curve = k_shot_curve(X, y, ids, group_of, fit_predict=ridge)
    assert curve.ks == (0, 1, 2, 4, 8, 16, 32)
    for c in curve.cells:
        assert not (set(c.shot_ids) & set(c.scored_ids)), (c.held_out_group, c.k)
        assert len(c.shot_ids) == c.k
        assert c.n_scored == len(c.scored_ids)


def test_a_cell_with_too_few_scorable_rows_is_unscorable_not_scored():
    X, y, ids, group_of = _panel(n_groups=3, per_group=45)
    curve = k_shot_curve(X, y, ids, group_of, fit_predict=ridge, ks=(0, 36))
    big = [c for c in curve.cells if c.k == 36]
    assert big and all(c.status == "unscorable" and c.score is None for c in big), \
        "45 - 36 = 9 scorable rows is under MIN_SCORED=10; a number here would come from 9 points"


def test_k_at_or_above_the_group_size_raises_rather_than_scoring_nothing():
    X, y, ids, group_of = _panel(n_groups=3, per_group=45)
    with pytest.raises(KShotGridError):
        k_shot_curve(X, y, ids, group_of, fit_predict=ridge, ks=(0, 45))


def test_the_shot_draw_is_PAIRED_across_methods_and_deterministic():
    """Same seed -> identical shots, so a comparison is paired. A difference of medians over
    DIFFERENT items is not a lift."""
    X, y, ids, group_of = _panel()
    a = k_shot_curve(X, y, ids, group_of, fit_predict=ridge, seed=7)
    b = k_shot_curve(X, y, ids, group_of, fit_predict=constant_mean, seed=7)
    for ca in a.cells:
        cb = b.cell(ca.held_out_group, ca.k)
        assert cb is not None and ca.shot_ids == cb.shot_ids
        assert ca.scored_ids == cb.scored_ids
    c = k_shot_curve(X, y, ids, group_of, fit_predict=ridge, seed=8)
    assert any(ca.shot_ids != c.cell(ca.held_out_group, ca.k).shot_ids
               for ca in a.cells if ca.k > 0), "a different seed must move the draw"


# --- the power check comes FIRST --------------------------------------------------------------------

def test_two_identical_methods_are_reported_UNDERPOWERED():
    """A zero-difference run necessarily has a zero gain, so checking the gain first would publish an
    unpowered run as a refutation."""
    X, y, ids, group_of = _panel()
    a = k_shot_curve(X, y, ids, group_of, fit_predict=ridge)
    b = k_shot_curve(X, y, ids, group_of, fit_predict=ridge)
    assert power_check(a, b) == UNDERPOWERED


def test_two_genuinely_different_methods_are_reported_POWERED():
    """Non-vacuity for the guard above: it must be able to say POWERED."""
    X, y, ids, group_of = _panel()
    a = k_shot_curve(X, y, ids, group_of, fit_predict=ridge)
    b = k_shot_curve(X, y, ids, group_of, fit_predict=constant_mean)
    assert power_check(a, b) == POWERED
