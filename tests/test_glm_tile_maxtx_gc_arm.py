"""The frozen bar for the max_tx+gc source-diversity arm, and the refusals that outrank a headline.

Every test here is offline and synthetic. The prediction tables live on a volume that is detached
from this host, so the arm cannot be RUN -- which is exactly why the bar is testable now: nothing
here can have been tuned to a result nobody can compute.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = _ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(name, mod)
    spec.loader.exec_module(mod)
    return mod


arm = _load("glm_tile_maxtx_gc_arm")


def _paired(median: float, wins: int, n: int = 10) -> dict:
    return {"median": median, "mean": median, "wins": wins, "n": n, "ci95": [median, median]}


def _arms(gain: float, wins: int, ctrl_gain: float = -0.12, ctrl_wins: int = 0) -> dict:
    return {
        "_paired_hybrid_vs_comparator": _paired(gain, wins),
        "_paired_control_vs_comparator": _paired(ctrl_gain, ctrl_wins),
    }


# --- the bar itself ----------------------------------------------------------------------------


def test_the_bar_mirrors_the_hybrid_arm_so_the_two_results_are_comparable():
    """If either margin or seed count drifts, the two combine results stop being comparable."""
    hybrid = _load("glm_tile_hybrid_arm")
    assert arm.PREREGISTERED["margin"] == hybrid.PREREGISTERED["margin"]
    assert arm.PREREGISTERED["min_seeds_winning"] == hybrid.PREREGISTERED["min_seeds_winning"]
    assert arm.PREREGISTERED["n_seeds_expected"] == hybrid.PREREGISTERED["n_seeds_expected"]
    assert arm.PREREGISTERED["comparator"] == hybrid.PREREGISTERED["comparator"] == "max_tx"


def test_the_combo_contains_no_learned_component():
    """The whole point: conv absent. A conv in here answers the hybrid arm's question, not this one."""
    assert arm.COMBO == ("max_tx", "gc")
    assert "conv" not in arm.COMBO


def test_the_prior_strike_is_declared_so_the_bar_cannot_read_as_expecting_a_lift():
    strike = arm.PREREGISTERED["prior_strike"]
    assert "+0.0302" in strike["fact"] and "+0.0339" in strike["fact"]
    assert "NO_GAIN" in strike["consequence"]


def test_it_declares_that_it_does_not_discharge_the_planned_stacking_arm():
    dnd = arm.PREREGISTERED["does_not_discharge"]
    assert dnd["arm"] == "orthogonal-modality-stack-promcalc-plus-N-channels"
    assert dnd["status_when_frozen"] == "planned"


def test_orthogonality_is_reported_and_never_gates_the_run():
    pre = arm.PREREGISTERED["precondition"]
    assert pre["measured_first"] is True
    assert "not_a_gate" in pre


# --- verdict: the refusals outrank the headline ---------------------------------------------------


def test_missing_scorer_refuses_before_anything_else():
    a = _arms(0.05, 10)
    a["_missing_scorer"] = ["gc"]
    a["_reconcile_problems"] = ["anchor drift"]
    v, detail = arm.verdict(a)
    assert v == "INDETERMINATE_MISSING_SCORER"
    assert "gc" in detail


def test_reconcile_failure_refuses_even_on_a_clean_lift():
    a = _arms(0.05, 10)
    a["_reconcile_problems"] = ["gc seed 3 = 0.9, not 0.3112"]
    v, _ = arm.verdict(a)
    assert v == "INDETERMINATE_RECONCILE_FAILED"


def test_a_lifting_control_outranks_a_passing_headline():
    """If shuffled gc lifts too, the number came from rank-averaging, not from gc."""
    v, detail = arm.verdict(_arms(0.05, 10, ctrl_gain=0.04, ctrl_wins=9))
    assert v == "INDETERMINATE_CONTROL_ALSO_LIFTS"
    assert "artifact" in detail


# --- verdict: the graded outcomes -----------------------------------------------------------------


def test_a_clean_lift_says_source_diversity_suffices():
    v, _ = arm.verdict(_arms(0.0339, 10))
    assert v == "SOURCE_DIVERSITY_SUFFICES"


def test_the_expected_outcome_is_a_real_verdict_not_a_failed_run():
    """The three-way already hints against a lift. NO_GAIN must be reachable and informative."""
    v, detail = arm.verdict(_arms(-0.03, 1))
    assert v == "NO_GAIN"
    assert "CONTENT" in detail


@pytest.mark.parametrize(
    "gain,wins",
    [(0.019, 10),   # clears the seed bar, misses the margin
     (0.05, 7)],    # clears the margin, misses the seed bar
)
def test_each_half_of_the_bar_bites_on_its_own(gain, wins):
    """Non-vacuity: neither margin nor seed count is decorative."""
    assert arm.verdict(_arms(gain, wins))[0] == "WEAK"


def test_exactly_at_the_bar_passes():
    assert arm.verdict(_arms(0.02, 8))[0] == "SOURCE_DIVERSITY_SUFFICES"


# --- mechanics ------------------------------------------------------------------------------------


def test_as_table_keys_are_positional_so_a_length_mismatch_cannot_be_papered_over():
    assert arm.as_table([0.5, 0.1, 0.9]) == {"0": 0.5, "1": 0.1, "2": 0.9}


def test_combine_returns_one_value_per_tile_in_the_original_order():
    preds = {"max_tx": [0.1, 0.5, 0.9, 0.3], "gc": [0.9, 0.5, 0.1, 0.7]}
    out = arm.combine(preds, arm.COMBO)
    assert len(out) == 4
    # Opposed rankings rank-average to a flat-ish middle; what matters is it is defined everywhere.
    assert all(isinstance(v, float) for v in out)


def test_the_control_permutes_and_never_invents_values():
    vec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    out = arm.shuffled(vec, 3)
    assert sorted(out) == sorted(vec)
    assert out != vec, "a control that returns the input unchanged is not a control"


def test_the_control_is_deterministic_per_seed():
    vec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    assert arm.shuffled(vec, 3) == arm.shuffled(vec, 3)


def test_the_control_differs_from_the_hybrid_arms_control_stream():
    """Distinct seed offsets, so the two arms' controls are not the same permutation."""
    hybrid = _load("glm_tile_hybrid_arm")
    vec = list(range(40))
    assert arm.shuffled(vec, 1) != hybrid.shuffled(vec, 1)


def test_orthogonality_of_a_vector_with_itself_is_one():
    vec = [0.1, 0.4, 0.2, 0.9, 0.5]
    assert arm.orthogonality(vec, vec) == pytest.approx(1.0)


# --- the absent-substrate refusal ------------------------------------------------------------------


def test_an_unreachable_substrate_refuses_with_exit_2_and_never_a_verdict(monkeypatch, capsys, tmp_path):
    """The volume holding the prediction tables is detached. That must read as INSUFFICIENT_DATA,
    not as a result -- a broken extraction wearing the costume of a finding."""
    missing = tmp_path / "not-attached"
    monkeypatch.setattr(sys, "argv", ["glm_tile_maxtx_gc_arm.py", "--in-dir", str(missing)])
    assert arm.main() == 2
    out = capsys.readouterr().out
    assert "INSUFFICIENT_DATA" in out
    assert "NO_GAIN" not in out and "SOURCE_DIVERSITY_SUFFICES" not in out


def test_an_attached_but_empty_substrate_also_refuses(monkeypatch, capsys, tmp_path):
    empty = tmp_path / "attached-empty"
    empty.mkdir()
    monkeypatch.setattr(sys, "argv", ["glm_tile_maxtx_gc_arm.py", "--in-dir", str(empty)])
    assert arm.main() == 2
    assert "INSUFFICIENT_DATA" in capsys.readouterr().out
