"""Label-threshold sweep — pure logic only (no solver, no feba.db, no model)."""
from __future__ import annotations

from dna_decode.fba.conditional_essentiality import GeneRecord
from scripts.fba_label_threshold_sweep import (
    call_mechanism,
    score_setting,
    verdict_for_sweep,
)

KEYS = ("a", "b", "c")


def test_a_nan_growth_is_an_infeasible_call_not_a_threshold_crossing():
    """The distinction the whole mechanism question turns on: an infeasible solve involves NO threshold,
    so no amount of threshold retuning can move it."""
    assert call_mechanism(None, 1.0) == "infeasible"
    assert call_mechanism(float("nan"), 1.0) == "infeasible"


def test_a_finite_ratio_below_the_bar_is_the_only_retunable_call():
    assert call_mechanism(0.001, 1.0) == "sub_threshold"        # 0.001 < 0.01 * 1.0
    assert call_mechanism(0.5, 1.0) == "not_essential"


def test_the_bar_scales_with_wildtype_growth():
    """1% of wild type, not 1% absolute -- a low-growth condition has a proportionally lower bar."""
    assert call_mechanism(0.004, 1.0) == "sub_threshold"
    assert call_mechanism(0.004, 0.2) == "not_essential"        # 0.004 > 0.01 * 0.2


def _rec(gid, exp):
    return GeneRecord(gid, gid, dict(zip(KEYS, exp, strict=True)), {}, True)


def test_a_constant_prediction_is_not_counted_as_a_commitment():
    """All-dispensable and all-essential are refusals to commit, not predictions."""
    subset = [_rec("gAll", (True, False, False)), _rec("gNone", (False, True, False))]
    calls = {"a": {"gAll": True, "gNone": False}, "b": {"gAll": True, "gNone": False},
             "c": {"gAll": True, "gNone": False}}          # gAll all-essential, gNone all-dispensable
    growth = {c: {"gAll": None, "gNone": 0.9} for c in KEYS}
    got = score_setting(subset, calls, KEYS, growth, {c: 1.0 for c in KEYS})
    assert got["n_committed"] == 0


def test_commitment_mechanism_is_split_by_infeasible_vs_sub_threshold():
    subset = [_rec("gInf", (True, False, False)), _rec("gSub", (False, True, False))]
    calls = {"a": {"gInf": True, "gSub": False},
             "b": {"gInf": False, "gSub": True},
             "c": {"gInf": False, "gSub": False}}
    growth = {"a": {"gInf": None, "gSub": 0.9}, "b": {"gInf": 0.9, "gSub": 0.001},
              "c": {"gInf": 0.9, "gSub": 0.9}}
    got = score_setting(subset, calls, KEYS, growth, {c: 1.0 for c in KEYS})
    assert got["n_committed"] == 2
    assert got["essential_calls_by_mechanism"] == {"infeasible": 1, "sub_threshold": 1}
    assert got["n_committed_all_infeasible"] == 1
    assert got["infeasible_share_of_commitment_calls"] == 0.5


def _cell(n=50, lift=0.05, share=0.9, tbar=None, mode=None, commit=5):
    return {"n_conditionally_essential": n, "lift_over_null": lift,
            "infeasible_share_of_commitment_calls": share, "min_abs_t": tbar,
            "min_abs_t_mode": mode, "n_committed": commit}


def test_the_mechanism_verdict_fires_only_when_feasibility_carries_nearly_everything():
    """If ~all commitment calls are infeasible, a threshold retune cannot help this metric at all."""
    assert verdict_for_sweep([_cell(share=0.97)])["mechanism"] ==         "CONDITIONAL_SIGNAL_IS_BINARY_FEASIBILITY"
    assert verdict_for_sweep([_cell(share=0.30)])["mechanism"] ==         "THRESHOLD_CROSSING_CARRIES_REAL_SIGNAL"


def test_the_fit_axis_verdict_reports_cutoff_dependence_honestly():
    assert verdict_for_sweep([_cell() for _ in range(4)])["headline_label_sensitivity"] ==         "SURVIVES_EVERY_POWERED_SETTING"
    mostly_bad = [_cell(lift=-0.01) for _ in range(3)] + [_cell(lift=0.05)]
    assert verdict_for_sweep(mostly_bad)["headline_label_sensitivity"] == "CUTOFF_DEPENDENT"


def test_a_degenerate_instrument_is_named_as_such_not_as_cutoff_dependence():
    """THE axis-pooling trap. The all_conditions t-bar removes every switcher, so its settings show
    zero commitments and a negative lift. Pooled with the sound axis that reads as 'the headline is
    cutoff-dependent' -- a false negative manufactured by a broken instrument."""
    cells = [_cell(lift=0.05, tbar=None)]
    cells += [_cell(lift=-0.2, tbar=2.0, mode="all_conditions", commit=0) for _ in range(3)]
    got = verdict_for_sweep(cells)
    assert got["headline_label_sensitivity"] == "SURVIVES_EVERY_POWERED_SETTING"   # fit axis is clean
    assert got["t_bar_axis_all_conditions"]["status"] == "INSTRUMENT_DEGENERATE_NO_COMMITMENTS"
    assert "ANTI-SELECTIVE" in got["note_on_all_conditions_mode"]


def test_an_underpowered_axis_refuses_to_render_a_verdict():
    """A cutoff that shrinks the set to a handful must not produce a confident-looking headline."""
    got = verdict_for_sweep([_cell(n=4, lift=0.9, share=1.0)])
    assert got["headline_label_sensitivity"] == "UNDERPOWERED"
