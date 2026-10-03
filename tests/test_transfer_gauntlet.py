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
    GauntletReport, MetricMatrix, build_matrix, group_offset_baseline, refuse_auroc_shaped_gate,
    run_gauntlet, shuffled_label_null, verdict,
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
                         | {"structure_only_baseline": 0.50})
    assert verdict(m, rep) == "LOSES_TO_structure_only_baseline"


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
