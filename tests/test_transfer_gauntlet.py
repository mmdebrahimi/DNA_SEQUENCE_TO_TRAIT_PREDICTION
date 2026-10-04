"""Guards for the five-control gauntlet and the four-cell verdict matrix (F1 Step 5)."""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.transfer import WITHIN_GROUP_MIN_N  # noqa: E402
from dna_decode.eval.transfer_gauntlet import (  # noqa: E402
    BEATS_ALL_CONTROLS, CONTROL_NAMES, GROUP_OFFSET_LEARNED, GROUP_OFFSET_WITHIN_GROUP_R2,
    WITHIN_GROUP_UNSCORABLE, WITHIN_NOISE, AurocShapedGateRefused, GauntletIncomplete,
    GauntletReport, MetricMatrix, build_matrix, group_offset_baseline, held_out_clade,
    nearest_neighbour_sequence, pca_only_baseline, refuse_auroc_shaped_gate, run_gauntlet,
    shuffled_label_null, structure_only_baseline, verdict,
)

GAUNTLET_SRC = ROOT / "dna_decode" / "eval" / "transfer_gauntlet.py"


def _panel(n_groups=4, per_group=80, signal=1.5, seed=0):
    rng = np.random.default_rng(seed)
    n = n_groups * per_group
    g = np.repeat([f"g{i}" for i in range(n_groups)], per_group)
    X = rng.integers(0, 2, (n, 3)).astype(float)
    y = np.repeat(rng.normal(0, 5, n_groups), per_group) + signal * X[:, 0] + rng.normal(0, .3, n)
    return X, y, g


# --- no statistics live in this module -------------------------------------------------------------

def test_the_gauntlet_imports_no_sklearn_or_scipy_at_ANY_scope():
    """The PCA control goes through numpy SVD precisely so this can hold. Non-vacuous: the same
    scanner must SEE the imports the module really has."""
    imported: set[str] = set()
    for node in ast.walk(ast.parse(GAUNTLET_SRC.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    offenders = [m for m in imported if m.split(".")[0] in {"sklearn", "scipy"}]
    assert not offenders, offenders
    assert "numpy" in imported and any("deconfound" in m for m in imported), imported


def test_the_gauntlet_fits_no_model_of_its_own():
    """The AST import guard alone would be satisfied by reaching sklearn's Ridge THROUGH deconfound's
    namespace -- which would pass the letter of the rule while breaking its purpose. All fitting must
    live inside deconfound's primitives or in the caller's `fit_predict`."""
    src = GAUNTLET_SRC.read_text(encoding="utf-8")
    code_lines = [ln for ln in src.splitlines()
                  if ln.strip() and not ln.strip().startswith("#")]
    code = "\n".join(code_lines)
    tree = ast.parse(src)
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Attribute) and n.func.attr == "fit"]
    assert not calls, "the gauntlet fits a model; move it to deconfound or the caller"
    assert "import Ridge" not in code, "Ridge pulled in via deconfound's namespace -- same violation"


# --- the three discriminating cases ----------------------------------------------------------------

def test_a_real_within_group_signal_beats_every_control():
    X, y, g = _panel(signal=1.5)
    rep = run_gauntlet(X, y, g, n_draws=50)
    assert rep.n_controls_run == 5
    rep.assert_all_controls_ran()
    rep.assert_same_evaluation_basis()
    assert verdict(rep.matrix, rep) == BEATS_ALL_CONTROLS
    assert rep.matrix.kshot_within_group_r2 > 0.5


def test_a_pure_structure_input_does_NOT_pass():
    """Group offsets only, zero within-group signal: the Arabidopsis failure mode in miniature."""
    X, y, g = _panel(signal=0.0, seed=1)
    rep = run_gauntlet(X, y, g, n_draws=50)
    v = verdict(rep.matrix, rep)
    assert v != BEATS_ALL_CONTROLS
    assert v in (WITHIN_NOISE, GROUP_OFFSET_LEARNED) or v.startswith("LOSES_TO_"), v
    assert rep.matrix.kshot_within_group_r2 <= GROUP_OFFSET_WITHIN_GROUP_R2 + 0.01


def test_the_within_group_cell_is_UNSCORABLE_at_20_per_group_and_never_falls_back():
    """THE MEASURED TRAP. `within_group_r2` returns (nan, 0) under min_n=30 while the pooled cell
    still reports a confident-looking number on the identical data. Falling back to pooled there is
    exactly how a group-offset run would be reported as a pass."""
    X, y, g = _panel(n_groups=6, per_group=20)
    rep = run_gauntlet(X, y, g, n_draws=30)
    assert rep.matrix.kshot_within_group_r2 is None
    assert rep.matrix.n_groups_used_within == 0
    assert rep.matrix.zero_shot_group_r2 is not None, "the pooled cell still has a number -- the trap"

    v = verdict(rep.matrix, rep)
    assert v == WITHIN_GROUP_UNSCORABLE
    assert not v.startswith("LOSES_TO_") and v not in (BEATS_ALL_CONTROLS, WITHIN_NOISE,
                                                       GROUP_OFFSET_LEARNED)


def test_GROUP_OFFSET_LEARNED_fires_when_pooled_wins_and_within_group_does_not():
    """A pooled-only win is a NAMED outcome, not a pass. Built as a unit test of the verdict function
    so the exact cell combination is unambiguous."""
    m = MetricMatrix(zero_shot_group_r2=-0.5, kshot_pooled_r2=0.40,
                     kshot_within_group_r2=-0.01, kshot_group_offset_baseline=0.20,
                     n_groups_used_within=4)
    rep = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES})
    assert verdict(m, rep) == GROUP_OFFSET_LEARNED


def test_a_zero_within_group_cell_with_no_pooled_advantage_is_WITHIN_NOISE():
    m = MetricMatrix(zero_shot_group_r2=-0.5, kshot_pooled_r2=0.10,
                     kshot_within_group_r2=0.0, kshot_group_offset_baseline=0.20,
                     n_groups_used_within=4)
    rep = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES})
    assert verdict(m, rep) == WITHIN_NOISE


def test_losing_to_a_named_control_names_it():
    m = MetricMatrix(zero_shot_group_r2=0.1, kshot_pooled_r2=0.1,
                     kshot_within_group_r2=0.20, kshot_group_offset_baseline=0.0,
                     n_groups_used_within=4)
    rep = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                         | {"held_out_clade": 0.50})
    assert verdict(m, rep) == "LOSES_TO_held_out_clade"


def test_structure_only_is_a_POOLED_DIAGNOSTIC_and_NOT_a_verdict_rival():
    """CHANGED 2026-10-03 on a measurement. Under leave-one-group-out the held-out group's one-hot
    column is all-zero in every training row, so the model cannot represent that group's offset at
    all -- it scores WORSE the more structure explains (-0.777 where group explains ~everything vs
    -0.001 where it explains nothing), presenting its easiest bar exactly where it was added to bite.
    The Arabidopsis comparison that motivated it was on POOLED metrics, not the within-group frame the
    verdict reads, where any purely structural predictor is identically 0 and so redundant with the
    offset baseline. It stays in the report as a diagnostic; it is not a rival."""
    m = MetricMatrix(zero_shot_group_r2=0.1, kshot_pooled_r2=0.1,
                     kshot_within_group_r2=0.20, kshot_group_offset_baseline=0.0,
                     n_groups_used_within=4)
    rep = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                         | {"structure_only_baseline": 0.99})
    assert verdict(m, rep) == BEATS_ALL_CONTROLS, \
        "structure_only_baseline must not gate the verdict"
    # but it is still COMPUTED and REPORTED -- demoted, not deleted
    assert "structure_only_baseline" in CONTROL_NAMES
    rep.assert_all_controls_ran()


# --- the two completeness guards, and which one is load-bearing ------------------------------------

def test_assert_all_controls_ran_raises_on_a_short_report():
    rep = GauntletReport(controls={n: 0.1 for n in CONTROL_NAMES[:4]})
    assert rep.n_controls_run == 4
    with pytest.raises(GauntletIncomplete):
        rep.assert_all_controls_ran()


def test_assert_same_evaluation_basis_catches_a_control_scored_on_DIFFERENT_rows():
    """THE GUARD THAT MATTERS. Four of five controls are imported, so five labels existing proves
    almost nothing; what can silently go wrong is one control scored on another id set. Tested by
    MUTATING the basis, not by changing the count."""
    X, y, g = _panel()
    rep = run_gauntlet(X, y, g, n_draws=30)
    rep.assert_same_evaluation_basis()                       # clean to start

    rep.basis["per_control"]["pca_only_baseline"]["scored_ids_hash"] = "deadbeefdeadbeef"
    with pytest.raises(GauntletIncomplete) as exc:
        rep.assert_same_evaluation_basis()
    assert "pca_only_baseline" in str(exc.value)


def test_an_unrecorded_basis_is_refused_rather_than_assumed_consistent():
    rep = GauntletReport(controls={n: 0.1 for n in CONTROL_NAMES})
    with pytest.raises(GauntletIncomplete):
        rep.assert_same_evaluation_basis()


# --- shape safety and provenance --------------------------------------------------------------------

def test_a_continuous_score_aimed_at_an_auroc_gate_is_refused():
    with pytest.raises(AurocShapedGateRefused):
        refuse_auroc_shaped_gate(0.73)
    refuse_auroc_shaped_gate({"cladeA": 0.8, "cladeB": 0.7})     # dict shape is fine


def test_the_auroc_gate_catches_a_numpy_scalar_but_NOT_a_0d_array():
    """MEASURED FINDING about the guard's reach. `np.floating` is in the isinstance tuple, so
    `np.float64(0.73)` is refused -- but a 0-d `ndarray` is neither `float` nor `np.floating` and
    passes straight through, and both shapes arise naturally here: every control returns
    `float(...)` while an unconverted `cv_r2` result would be a numpy value.

    `bool` is caught incidentally (it subclasses `int`), which is correct -- a boolean is not a
    per-clade AUROC mapping either. Tripwire: adding `np.ndarray` to the check fails this test.
    """
    for caught in (0.73, 1, True, np.float64(0.73), np.float32(0.73)):
        with pytest.raises(AurocShapedGateRefused):
            refuse_auroc_shaped_gate(caught)
    refuse_auroc_shaped_gate(np.array(0.73))        # the escape -- the finding
    refuse_auroc_shaped_gate(np.array([0.73, 0.8]))


def test_the_report_records_which_primitives_produced_its_numbers():
    """Drift between the replay arm and the live arms must be detectable FROM THE ARTIFACT."""
    X, y, g = _panel()
    prov = run_gauntlet(X, y, g, n_draws=30).primitive_provenance
    assert len(prov["deconfound_module_sha256"]) == 64
    assert any("within_group_r2" in q for q in prov["qualnames"])
    assert prov["within_group_min_n_passed_explicitly"] == WITHIN_GROUP_MIN_N


def test_n_groups_used_is_reported_in_both_the_scorable_and_unscorable_case():
    X, y, g = _panel(n_groups=4, per_group=80)
    assert run_gauntlet(X, y, g, n_draws=30).matrix.n_groups_used_within == 4
    X2, y2, g2 = _panel(n_groups=6, per_group=20)
    assert run_gauntlet(X2, y2, g2, n_draws=30).matrix.n_groups_used_within == 0


def test_the_null_fails_CLOSED_when_no_binary_feature_qualifies():
    """`univariate_top` needs >=5 units per arm. With no usable feature the null is not computable,
    so the control is nan -> unrun -> assert_all_controls_ran raises. It must not return a number."""
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 2))                             # continuous, never 0/1
    y = rng.normal(size=60)
    g = np.repeat(["a", "b", "c"], 20)
    assert np.isnan(shuffled_label_null(X, y, g, n_draws=20))


def test_a_NaN_control_counts_as_UNRUN_and_makes_assert_all_controls_ran_raise():
    """FIXED 2026-10-03 (this test previously recorded the fail-OPEN). Both the counter and the check
    tested `is None`, and float('nan') is not None -- so a control that could not be computed
    reported as RUN and the report read clean, contradicting the fail-closed behaviour
    `shuffled_label_null`'s own docstring promised. nan is the fail-closed path of the one control the
    module builds itself, so the shape most likely to produce it was least likely to be noticed."""
    rep = GauntletReport(controls={n: 0.1 for n in CONTROL_NAMES}
                         | {"shuffled_label_null": float("nan")})
    assert rep.n_controls_run == 4, "a NaN control must not count as run"
    with pytest.raises(GauntletIncomplete) as exc:
        rep.assert_all_controls_ran()
    assert "shuffled_label_null" in str(exc.value)

    with pytest.raises(GauntletIncomplete):
        GauntletReport(controls={n: 0.1 for n in CONTROL_NAMES}
                       | {"shuffled_label_null": None}).assert_all_controls_ran()
    # non-vacuity: a fully-computed report still passes
    GauntletReport(controls={n: 0.1 for n in CONTROL_NAMES}).assert_all_controls_ran()


def test_the_shuffled_label_null_GATES_the_verdict_on_MATCHED_units():
    """FIXED 2026-10-03, and the dead branch hid a worse defect underneath it.

    The null was consulted behind `w <= 0`, but the branch above already returned on
    `w <= GROUP_OFFSET_WITHIN_GROUP_R2` which IS 0.0 -- so it was unreachable and the control gated
    nothing. Underneath that: `permutation_null` returns a SPEARMAN RHO while the cell it was compared
    against is an R-SQUARED (measured 0.0614 vs 0.7237 on one fixture), so merely making the branch
    reachable would have compared incommensurable quantities -- the same class as the
    AUROC-vs-continuous trap `refuse_auroc_shaped_gate` exists to stop. `shuffled_label_null` now
    returns a within-group r2 null, so the two are the same quantity."""
    m = MetricMatrix(zero_shot_group_r2=0.1, kshot_pooled_r2=0.1, kshot_within_group_r2=0.05,
                     kshot_group_offset_baseline=0.0, n_groups_used_within=4)
    low = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES})
    high = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                          | {"shuffled_label_null": 0.99})
    assert verdict(m, low) == BEATS_ALL_CONTROLS
    assert verdict(m, high) == WITHIN_NOISE, "the null must now gate the verdict"

    # it gates on the MARGIN, not merely on presence: a null just under the candidate still passes
    near = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                          | {"shuffled_label_null": 0.049})
    assert verdict(m, near) == BEATS_ALL_CONTROLS
    rival = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                           | {"pca_only_baseline": 0.99})
    assert verdict(m, rival) == "LOSES_TO_pca_only_baseline"


def test_excluding_SELF_is_what_stops_the_1NN_control_being_trivially_perfect():
    """The load-bearing default, tested directly because `run_gauntlet` never passes the flag.

    With `exclude_self=False` the nearest neighbour of every row IS that row, so the control predicts
    y from y and scores EXACTLY 1.0 -- it would then beat any candidate and every verdict would come
    back LOSES_TO_nearest_neighbour_sequence regardless of the data. The diagonal fill is the whole
    control, so it gets its own assertion rather than being implied by an end-to-end pass."""
    y = np.array([1.0, 2.0, 3.0, 10.0])
    dist = np.abs(y[:, None] - y[None, :])
    assert nearest_neighbour_sequence(y, dist, exclude_self=False) == 1.0
    honest = nearest_neighbour_sequence(y, dist)
    assert honest < 1.0, "self-exclusion must cost the control its free perfect score"

    # and the degenerate control really would capture the verdict
    m = MetricMatrix(zero_shot_group_r2=0.1, kshot_pooled_r2=0.1, kshot_within_group_r2=0.80,
                     kshot_group_offset_baseline=0.0, n_groups_used_within=4)
    rep = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                         | {"nearest_neighbour_sequence": 1.0})
    assert verdict(m, rep) == "LOSES_TO_nearest_neighbour_sequence"


def test_over_asking_the_pca_control_for_components_is_safe_and_the_explicit_min_is_REDUNDANT():
    """Behavioural pin plus a correction to my own first draft of it.

    The panels here carry 3 features against `n_components=5`, so every end-to-end test above already
    over-asks. That is safe -- but NOT because of the explicit `min(n_components, min(Xc.shape))`:
    numpy slice bounds already clamp, so `vt[:50]` on a rank-3 decomposition is `vt[:3]`. Verified by
    deleting the `min()` and re-running: nothing failed. So the line is defensive and inert, and a test
    named for it would have been a guard over a mechanism that cannot break.

    What IS worth pinning is the contract a caller depends on: over-asking returns a finite number
    equal to asking for exactly the available rank, rather than raising or silently changing the score.
    """
    X, y, g = _panel()
    assert X.shape[1] == 3, "this test is about asking for more components than exist"
    at_rank = pca_only_baseline(X, y, g, n_components=3)
    assert np.isfinite(at_rank)
    assert pca_only_baseline(X, y, g, n_components=50) == at_rank
    assert pca_only_baseline(X, y, g, n_components=5) == at_rank, \
        "the default over-asks on these panels; it must agree with the at-rank score"
    # and the control is not degenerate: fewer components than rank DOES change the answer
    assert pca_only_baseline(X, y, g, n_components=1) != at_rank


def test_the_structure_only_control_is_ANTI_informative_under_the_grouped_CV_it_is_called_with():
    """MEASURED FINDING, and it inverts the control's stated purpose.

    `structure_only_baseline` builds one-hot GROUP membership and scores it with
    `cv_r2(..., groups=groups)` -- leave-one-group-out. But the held-out group's own indicator column
    is all-zero across every training row (measured below: 0 rows), so ridge has no coefficient for it
    and predicts roughly the training mean -- a mean computed with that group's offset EXCLUDED.

    The consequence is backwards. Measured on identical group structure:

        phenotype IS the group offset (one-hot explains ~0.98 under K-fold) -> -0.777
        phenotype independent of group                                      -> -0.001

    It scores WORSE the MORE population structure explains. And in `verdict` this control is a RIVAL
    the candidate must exceed (`score <= rival -> LOSES_TO_...`), so it presents its EASIEST bar in
    exactly the situation it was added to catch. The Arabidopsis comparison that motivated it -- 'the
    control that BEAT the embedding' -- was a POOLED/global one; under leave-one-group-out the
    comparison it encodes cannot be made at all.

    Tripwire, not an endorsement: scoring this control on a group-blind split (or dropping the
    held-out indicator) would make the first assertion fail, which is the fix.
    """
    _, y_offset, g = _panel(signal=0.0, seed=3)
    rng = np.random.default_rng(3)
    y_noise = rng.normal(size=len(g))

    explains_everything = structure_only_baseline(y_offset, g)
    explains_nothing = structure_only_baseline(y_noise, g)
    assert explains_everything < 0.0, (
        f"the structure-only control now returns {explains_everything:.4f} on a pure group offset -- "
        f"it is no longer anti-informative under grouped CV. That is the fix; update this test.")
    assert explains_everything < explains_nothing, (
        "the inversion is the finding: a phenotype fully explained by group must not score BELOW a "
        "phenotype independent of it")

    # the mechanism, asserted rather than asserted-about: nothing trains the held-out indicator
    levels = sorted(set(g.tolist()))
    onehot = np.zeros((len(g), len(levels)))
    for j, lv in enumerate(levels):
        onehot[g == lv, j] = 1.0
    for lv in levels:
        assert onehot[g != lv, levels.index(lv)].sum() == 0.0

    # non-vacuity: the SAME features on a group-blind split recover the offset, so the features are
    # informative and it is the splitting convention that destroys the control.
    from dna_decode.deconfound import cv_r2
    assert cv_r2(onehot, y_offset) > 0.9


def test_the_held_out_clade_control_derives_its_folds_from_the_DISTANCE_matrix():
    """Non-vacuity by construction: three well-separated blocks in the distance matrix must produce a
    different number than one in which the same rows are mutually equidistant."""
    X, y, g = _panel(n_groups=3, per_group=40, signal=1.5, seed=5)
    d2 = ((X[:, None, :] - X[None, :, :]) ** 2).sum(-1)
    euclid = np.sqrt(np.maximum(d2, 0.0))
    blocked = (g[:, None] != g[None, :]).astype(float)      # 0 within group, 1 across
    assert np.isfinite(held_out_clade(X, y, euclid))
    assert held_out_clade(X, y, blocked) != held_out_clade(X, y, euclid), \
        "the clade folds must depend on the supplied distances, not only on X"


def test_the_offset_baseline_with_NO_shots_degenerates_to_the_global_mean():
    """The k=0 cell. With no shots there is nothing group-specific to predict, so every row gets the
    global mean and pooled r2 is exactly 0 -- which is what makes the k>0 value above it meaningful."""
    y = np.array([1.0, 2.0, 3.0, 10.0])
    assert group_offset_baseline(y, np.array(["a", "a", "b", "b"])) == 0.0
    assert group_offset_baseline(y, np.array(["a", "a", "b", "b"]), {}) == 0.0


def test_the_offset_baselines_within_group_r2_is_zero_BY_CONSTRUCTION():
    """A predictor constant inside each group explains none of the within-group variance, which is
    why 'beat the offset baseline on the within-group cell' reduces to 'within-group r2 > 0'."""
    assert GROUP_OFFSET_WITHIN_GROUP_R2 == 0.0
    X, y, g = _panel()
    pooled = group_offset_baseline(y, g, {gg: list(np.where(g == gg)[0][:4]) for gg in set(g)})
    assert pooled > 0.5, "the offset baseline should score well POOLED -- that is the whole danger"


def test_build_matrix_takes_the_kshot_pooled_cell_from_the_CALLER():
    """Computing it here would mean fitting a model in this module; the caller already owns the
    estimator via `fit_predict`."""
    X, y, g = _panel()
    m0 = build_matrix(X, y, g)
    assert m0.kshot_pooled_r2 == m0.zero_shot_group_r2, "with no k-shot run, k=0 is the honest value"
    m1 = build_matrix(X, y, g, kshot_pooled=0.123)
    assert m1.kshot_pooled_r2 == 0.123 and m1.zero_shot_group_r2 == m0.zero_shot_group_r2
