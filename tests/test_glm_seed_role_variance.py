"""The variance decomposition must have three genuinely different arms, and must refuse a degenerate one.

The experiment is three lines of seed wiring, so that is where it can silently stop measuring anything:
two arms that share a seed plan would produce two identical SDs and a plausible-looking decomposition. The
arm wiring is therefore a pure function, pinned here directly, and the script refuses rather than reports
when an arm collapses to a single repeated value.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_seed_role_variance.py"
OUT = ROOT / "wiki" / "glm_seed_role_variance_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_seedvar", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the arm wiring

def test_arm_A_moves_all_three_seeds_together(mod):
    assert mod.seed_plan("A_all_coupled", 3) == (3, 3, 3)
    assert mod.seed_plan("A_all_coupled", 9) == (9, 9, 9)


def test_arm_B_holds_corpus_and_split_constant_while_generation_varies(mod):
    a = mod.seed_plan("B_generation_only", 1, anchor=7)
    b = mod.seed_plan("B_generation_only", 5, anchor=7)
    assert a == (7, 1, 7) and b == (7, 5, 7)
    assert a[0] == b[0] and a[2] == b[2], "the corpus or split moved in the generation-only arm"
    assert a[1] != b[1], "the generation seed did not move"


def test_arm_C_is_the_converse_of_arm_B(mod):
    a = mod.seed_plan("C_corpus_only", 1, anchor=7)
    b = mod.seed_plan("C_corpus_only", 5, anchor=7)
    assert a == (1, 7, 1) and b == (5, 7, 5)
    assert a[1] == b[1], "the generation seed moved in the corpus-only arm"
    assert a[0] != b[0] and a[2] != b[2], "the corpus or split did not move"


def test_an_unknown_arm_is_refused_not_defaulted(mod):
    with pytest.raises(ValueError):
        mod.seed_plan("D_made_up", 3)


def test_the_three_arms_are_not_accidentally_identical(mod):
    """A decomposition whose arms share one plan measures nothing and must fail loudly."""
    assert mod.arms_are_distinct(list(range(10))) is True


def test_arms_collapse_detectably_when_there_is_nothing_to_vary(mod):
    """Non-vacuity for the distinctness guard: on a single-seed list equal to the anchor, all three
    arms reduce to the same plan, and the guard must SAY so rather than pass."""
    assert mod.arms_are_distinct([7], anchor=7) is False


# ------------------------------------------------------------------ the mechanism helpers

def test_independent_components_give_an_additivity_ratio_near_one(mod):
    """Control for the mechanism check: uncorrelated arms should NOT look like cancellation."""
    a = [0.50, 0.52, 0.48, 0.51, 0.49, 0.53, 0.47, 0.50, 0.52, 0.48]
    out = mod.variance_additivity(a, a, [0.5] * 9 + [0.5000001])
    assert out["sum_over_coupled_ratio"] == pytest.approx(1.0, abs=0.05), out


def test_a_coupled_variance_far_below_the_sum_is_flagged_as_cancellation(mod):
    """The real shape: tight coupled arm, wide decoupled arms -> a large ratio."""
    coupled = [0.50, 0.505, 0.495, 0.50, 0.505, 0.495, 0.50, 0.505, 0.495, 0.50]
    wide_b = [0.46, 0.54, 0.47, 0.53, 0.48, 0.52, 0.49, 0.51, 0.45, 0.55]
    wide_c = [0.47, 0.53, 0.48, 0.52, 0.46, 0.54, 0.50, 0.50, 0.49, 0.51]
    out = mod.variance_additivity(coupled, wide_b, wide_c)
    assert out["sum_over_coupled_ratio"] > 3, out


def test_the_verdict_bands_are_stated_and_ordered(mod):
    assert mod.verdict(0.0266, 0.0100) == "DECOUPLING_WOULD_HALVE_THE_SPREAD_OR_BETTER"
    assert mod.verdict(0.0266, 0.0180) == "DECOUPLING_WOULD_HELP_MODESTLY"
    assert mod.verdict(0.0140, 0.0280) == "DECOUPLING_LARGELY_INERT"
    assert mod.verdict(None, 0.02) == "INDETERMINATE"
    assert mod.verdict(0.0, 0.02) == "INDETERMINATE"


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_every_arm_produced_more_than_one_distinct_value():
    """An SD computed over a repeated number is not a measurement, and must never be reported as one."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for arm, modes in rec["summary"].items():
        for mode, s in modes.items():
            assert s["n_distinct"] >= 2, (arm, mode, s["n_distinct"])


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_artifact_carries_the_mechanism_and_the_refuted_rival():
    """The surprising direction of this result needs a mechanism and a named rival, not a shrug."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["mechanism_variance_additivity"]["both"]["sum_over_coupled_ratio"] > 1
    sens = rec["mechanism_anchor_sensitivity"]
    assert len(sens["by_anchor"]) >= 3, "one anchor cannot refute an anchor-specificity rival"
    assert sens["rival_status"].startswith("REFUTED")


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_artifact_scopes_itself_to_the_markov_probe():
    """It bounds the INSTRUMENT. Reading it as a statement about the HF generator is the live misreading."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert "NOT the HF generator" in rec["probe_candidate"]
    assert any("HF generator" in lim for lim in rec["honest_limits"])
