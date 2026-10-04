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
    MIN_SCORED, NOTHING_COMPARED, POWERED, UNDERPOWERED, WITHIN_GROUP_MIN_N, KShotCell, KShotCurve, KShotGridError,
    LeakageAuditRefused, UnassignedMemberError, derive_k_grid, fold_report, held_out_group_folds,
    k_shot_curve, leakage_audit, power_check, require_clean_audit, skipped_groups,
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


def test_fold_report_returns_the_SAME_SHAPE_for_an_empty_and_a_populated_partition():
    """FIXED 2026-10-03 (this test previously recorded the gap). The empty branch returned four keys
    while the populated branch returned six, so a consumer indexing `min_group_size` got a KeyError
    on an empty partition instead of the None its siblings return. Both branches now carry six."""
    empty = fold_report([])
    populated = fold_report(held_out_group_folds({"a": "g0", "b": "g1"}))
    assert set(empty) == set(populated), f"shape gap returned: {set(populated) ^ set(empty)}"
    assert empty["n_folds"] == 0 and empty["n_members"] == 0
    for k in ("largest_group_fraction", "singleton_fraction", "min_group_size", "max_group_size"):
        assert empty[k] is None, f"{k} should be None on an empty partition, got {empty[k]!r}"
    assert populated["min_group_size"] == 1 and populated["max_group_size"] == 1


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


def test_the_shot_draw_is_reproducible_ACROSS_processes_from_the_recorded_seed():
    """FIXED 2026-10-03 (this test previously recorded the gap). `_shot_ids` keyed its RNG on
    `abs(hash((seed, tuple(test_ids))))`, and hash() over strings is PYTHONHASHSEED-salted, so four
    fresh interpreters drew four different shot sets. PAIRING always held -- that is what the
    comparison's validity rests on -- but re-derivable-from-the-recorded-seed did not, which is the
    half an artifact's `seed` field promises a later reader. Now keyed on sha256 of the same inputs.

    Deliberately run in SUBPROCESSES with PYTHONHASHSEED unset: an in-process check cannot see the
    salt at all, which is why the defect survived the original same-process determinism test."""
    import os
    import subprocess
    code = ("import sys; sys.path.insert(0, '.');"
            "from dna_decode.eval.transfer import _shot_ids;"
            "print(_shot_ids(['a','b','c','d','e','f','g','h'], 3, 0))")

    def draw():
        env = dict(os.environ)
        env.pop("PYTHONHASHSEED", None)
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             cwd=str(ROOT), env=env)
        assert out.returncode == 0, out.stderr
        return out.stdout.strip()

    draws = {draw() for _ in range(4)}
    assert len(draws) == 1, f"shot draw is not reproducible across processes: {draws}"
    assert draws != {"()"}, "non-vacuity: the draw must be non-empty at k=3"


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


def test_a_cell_SCORED_on_one_side_and_UNSCORABLE_on_the_other_is_POWERED():
    """The third branch, unreached by the two runs above: statuses differing is itself a difference,
    because one method produced a number where the other produced none.

    ISOLATION MATTERS HERE, and the first version of this test did not have it. With the SCORED cell
    as `curve_a` the prediction comparison fires too ((1.0, 2.0) != ()), so deleting the status branch
    entirely left the test green -- it was passing through the wrong branch. Verified by injecting
    exactly that deletion. Driving it from the UNSCORABLE side skips the `a.status == "scored"` guard,
    so the status check is the only thing that can return POWERED.
    """
    scored = KShotCurve(ks=(0,), cells=[
        KShotCell("g0", 0, 0.5, 40, "scored", (), (), (1.0, 2.0))])
    unscorable = KShotCurve(ks=(0,), cells=[
        KShotCell("g0", 0, None, 3, "unscorable", (), (), ())])
    assert power_check(unscorable, scored) == POWERED, \
        "the status-differs branch is the only path to POWERED from the unscorable side"
    assert power_check(scored, unscorable) == POWERED, "and it must be symmetric"


def test_power_check_separates_NOTHING_COMPARED_from_methods_that_never_differed():
    """FIXED 2026-10-03 (this test previously recorded the gap). Two curves whose (group, k) cells
    never line up compared nothing and fell out of the loop as UNDERPOWERED -- the identical token a
    genuine zero-difference result gets. A checker built so an unpowered run is not published as a
    refutation must not itself conflate a plumbing failure with a measured zero."""
    a = KShotCurve(ks=(0,), cells=[KShotCell("g0", 0, 0.5, 40, "scored", (), (), (1.0, 2.0))])
    disjoint = KShotCurve(ks=(0,), cells=[KShotCell("g9", 0, 0.9, 40, "scored", (), (), (9.0, 9.0))])
    assert power_check(a, disjoint) == NOTHING_COMPARED
    assert power_check(KShotCurve(ks=()), KShotCurve(ks=())) == NOTHING_COMPARED

    # non-vacuity: a real comparison still yields the measured tokens, so the new branch has not
    # swallowed the ones it sits beside.
    same = KShotCurve(ks=(0,), cells=[KShotCell("g0", 0, 0.5, 40, "scored", (), (), (1.0, 2.0))])
    assert power_check(a, same) == UNDERPOWERED
    differ = KShotCurve(ks=(0,), cells=[KShotCell("g0", 0, 0.9, 40, "scored", (), (), (3.0, 4.0))])
    assert power_check(a, differ) == POWERED


def test_mismatched_ids_X_y_lengths_raise_rather_than_scoring_a_misaligned_panel():
    """`ids[i]` NAMES row i; a length mismatch silently re-labels rows, so a group vector would point
    at the wrong genotypes while every downstream number still looked well-formed."""
    X, y, ids, group_of = _panel(n_groups=2, per_group=45)
    with pytest.raises(ValueError, match="length mismatch"):
        k_shot_curve(X, y[:-1], ids, group_of, fit_predict=ridge, ks=(0,))
    with pytest.raises(ValueError, match="length mismatch"):
        k_shot_curve(X, y, ids[:-1], group_of, fit_predict=ridge, ks=(0,))


# --- the leakage audit fails CLOSED (Step 4) --------------------------------------------------------

def test_exact_id_overlap_is_an_immediate_failure():
    a = leakage_audit(["x", "y", "z"], ["z", "w"])
    assert a.overlap == ("z",) and not a.ok and not a.incomplete
    with pytest.raises(LeakageAuditRefused):
        require_clean_audit(a)


def test_disjoint_without_a_resolver_passes_but_SAYS_aliasing_was_unchecked():
    """The accession-string check cannot see a GCA assembly and an ERR run naming one isolate. Passing
    is correct; silently implying aliasing was ruled out would not be."""
    a = leakage_audit(["x", "y"], ["w"])
    assert a.ok and not a.incomplete
    assert "NOT checked" in a.reason
    assert require_clean_audit(a)["leakage_audit_degraded"] is False


def test_a_resolver_that_raises_makes_the_audit_INCOMPLETE_not_clean():
    """A membership question we could not answer is not a negative answer."""
    def boom(i):
        raise RuntimeError("entrez down")

    a = leakage_audit(["x", "y"], ["w"], resolver=boom)
    assert a.incomplete and not a.ok
    assert set(a.unresolved) == {"x", "y", "w"}
    assert "UNPROVEN" in a.reason


def test_a_resolver_returning_blank_also_counts_as_unresolved():
    a = leakage_audit(["x"], ["w"], resolver=lambda i: "" if i == "w" else "S1")
    assert a.incomplete and a.unresolved == ("w",)


def test_a_resolver_collapses_aliases_and_catches_a_hidden_overlap():
    """The case the exact-id check is blind to: two different ids, one isolate."""
    alias = {"GCA_1": "SAMN1", "ERR_1": "SAMN1", "GCA_2": "SAMN2"}
    a = leakage_audit(["GCA_1", "GCA_2"], ["ERR_1"], resolver=alias.get)
    assert a.overlap == ("SAMN1",) and not a.ok
    with pytest.raises(LeakageAuditRefused):
        require_clean_audit(a)


def test_a_clean_resolved_audit_passes_and_is_not_degraded():
    alias = {"GCA_1": "SAMN1", "GCA_2": "SAMN2", "ERR_3": "SAMN3"}
    a = leakage_audit(["GCA_1", "GCA_2"], ["ERR_3"], resolver=alias.get)
    assert a.ok and not a.incomplete and a.overlap == ()
    assert require_clean_audit(a)["leakage_audit_degraded"] is False


def test_the_override_is_VISIBLE_in_the_record_rather_than_making_it_clean():
    def boom(i):
        raise RuntimeError("down")

    a = leakage_audit(["x"], ["w"], resolver=boom)
    with pytest.raises(LeakageAuditRefused):
        require_clean_audit(a)                               # refuses by default
    stamp = require_clean_audit(a, allow_incomplete=True)
    assert stamp["leakage_audit_degraded"] is True
    assert stamp["leakage_audit"]["incomplete"] is True, "the override must not rewrite the audit"


def test_the_override_cannot_wave_through_a_REAL_overlap():
    """`allow_incomplete` is about an unproven audit, not a failed one."""
    a = leakage_audit(["z"], ["z"])
    with pytest.raises(LeakageAuditRefused):
        require_clean_audit(a, allow_incomplete=True)


def _imported_names(path: Path) -> set[str]:
    """Every module name this file imports, at ANY scope (module level or inside a function)."""
    import ast
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def test_the_frozen_cohort_manifest_is_never_IMPORTED_by_this_module():
    """cohort_manifest.py is one of the five sha256-pinned frozen files. The concern is IMPORTING it,
    not mentioning it, so this walks the AST rather than grepping -- a substring scan would both
    false-positive on the docstring and miss nothing a rephrase could hide. Non-vacuous: the same
    scan must SEE the imports this module really has."""
    imported = _imported_names(ROOT / "dna_decode" / "eval" / "transfer.py")
    assert not any("cohort_manifest" in m for m in imported), imported
    assert "numpy" in imported and "dna_decode.deconfound" in imported, \
        f"the scanner found nothing recognisable; it is not actually reading imports: {imported}"


def test_the_harness_imports_deconfound_LAZILY_not_at_module_scope():
    """The import-posture guard above proves no heavy module is pulled; this proves WHY -- the
    deconfound import sits inside a function body, not at module level."""
    import ast
    tree = ast.parse((ROOT / "dna_decode" / "eval" / "transfer.py").read_text(encoding="utf-8"))
    top_level = {a.name if isinstance(n, ast.Import) else n.module
                 for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))
                 for a in (n.names if isinstance(n, ast.Import) else [None]) if True}
    assert not any(m and "deconfound" in m for m in top_level), \
        f"deconfound imported at module scope: {top_level}"
