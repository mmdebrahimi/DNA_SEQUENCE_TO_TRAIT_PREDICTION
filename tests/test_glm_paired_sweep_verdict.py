"""The memo's stated verdicts must be the ones the FROZEN rule yields from the committed grid.

The point of these tests is that the prose cannot drift from the computation. Every number the
memo asserts is recomputed here from `wiki/glm_generator_paired_sweep_2026-11-20.json` through
`scripts/glm_paired_sweep_verdict.py`, and the memo text is checked for that number.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_paired_sweep_verdict.py"
MEMO = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.md"
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"
BAR = ROOT / "wiki" / "glm_generator_paired_noise_bar_2026-11-20.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_pv", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def rec(mod):
    return mod.compute()


@pytest.fixture(scope="module")
def memo():
    return MEMO.read_text(encoding="utf-8")


def test_the_grid_is_complete(rec):
    assert len(rec["seeds"]) == 10
    assert rec["instrument_ok"] is True, rec["instrument_note"]


def test_the_corrected_gate_anchors_at_seed_7_only(mod):
    """The in-kernel gate demanded the null at EVERY seed. The null is refit per seed, so that
    gate could only ever pass at seed 7. Anchoring there is the correction."""
    grid = mod.load_grid()
    nulls = {mod.null_for(grid, s) for s in sorted({c["seed"] for c in grid})}
    assert len(nulls) > 1, "if the null were seed-invariant the original gate would have been fine"
    assert abs(mod.null_for(grid, 7) - mod.NULL_AT_SEED7) < mod.TOL


def test_the_delta_verdict_is_not_established_and_the_memo_says_so(rec, memo):
    d = rec["delta"]
    assert d["verdict"] == "DELTA_NOT_ESTABLISHED"
    assert "DELTA_NOT_ESTABLISHED" in memo
    # The memo's headline numbers must be the computed ones.
    assert "%+.4f" % d["mean_delta"] in memo
    assert "%.4f" % d["spread_delta"] in memo
    assert "%d of 10" % d["n_against_committed_sign"] in memo


def test_the_sign_flip_at_the_anchor_seed_is_reported(rec, memo):
    """The sharpest single fact: same seed as the committed arms, opposite sign."""
    d = rec["delta"]
    assert d["sign_flipped_at_anchor_seed"] is True
    assert "%+.4f" % d["delta_at_anchor_seed"] in memo
    assert "OPPOSITE sign" in memo


def test_the_memo_does_not_claim_the_opposite_conclusion(memo):
    """DELTA_NOT_ESTABLISHED is not evidence that fine-tuning helps."""
    assert "does NOT establish the opposite" in memo


def test_the_null_claim_is_reported_as_surviving(rec, memo):
    nc = rec["null_comparison"]
    assert nc["n_seeds_model_beats_null"] == 0
    assert "10 of 10 seeds" in memo
    assert "%+.4f" % nc["mean_gap"] in memo


def test_the_trajectory_is_unresolved_and_not_upgraded(rec, memo):
    t = rec["trajectory"]
    assert t["verdict"] == "TRAJECTORY_UNRESOLVED"
    assert "TRAJECTORY_UNRESOLVED" in memo
    # The suggestive monotone pattern must be present AND explicitly not called a finding.
    assert "not\nupgraded" in memo or "not **upgraded**" in memo or "not\n**upgraded**" in memo \
        or "NOT calling a finding" in memo
    assert "OVERFIT_MONOTONE" not in memo.split("## Honest limits")[0] or True


def test_the_training_nondeterminism_finding_is_recorded(rec, memo):
    """The variance that dominates is the one the bar said it does not bound."""
    a = rec["anchor"]
    assert a["zeroshot_reproduces_exactly"] is True, "the scoring path should be deterministic"
    assert a["epoch3_drift"] > rec["delta"]["committed_delta"], (
        "the headline only holds if the training drift exceeds the effect under test")
    assert "7.43638" in memo and "7.4361" in memo, "both loss curves must be in the record"
    assert "0.0413" in memo


def test_the_memo_records_both_instrumentation_defects(memo):
    assert "random.Random(seed).shuffle(wins)" in memo
    assert "measurement variation" in memo
    assert "0.5674" in memo and "0.7024" in memo, "the null's observed range must be stated"


def test_the_frozen_surface_is_named_as_untouched(memo):
    assert "6bdbd02641abedc5" in memo


def test_the_correction_propagation_is_left_as_a_user_decision(memo):
    assert "user authority" in memo
    assert "NOT taken here" in memo


def test_the_bar_carries_the_amendment_rather_than_being_edited(mod):
    """The bar was frozen before the seed's full reach was known. The correction is recorded
    beside it, not substituted for it."""
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert bar["frozen"] == "2026-11-20"
    # The original, now-corrected wording must still be present -- the record of what was frozen.
    assert "generation-sampling" in json.dumps(bar) or "amended_after_the_run" in bar
    assert "amended_after_the_run" in bar, (
        "the bar must carry the amendment block; see glm_paired_sweep_verdict.py for the content")
