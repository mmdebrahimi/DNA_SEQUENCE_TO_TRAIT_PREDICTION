"""Offline tests for the nitrogen conditional-essentiality axis.

Pure/contract tests only -- no feba.db, no cobra solves, so these stay green without D: attached.
"""
from __future__ import annotations

import pytest

from dna_decode.fba.nitrogen import (
    DEFAULT_CARBON,
    NITROGEN_EXCHANGES,
    NITROGEN_UNMAPPABLE,
    REDACTED_ON_NONDETERMINISM,
    apply_nitrogen_condition,
    determinism_verdict,
    nitrogen_conditions,
    redact_unverified,
)

_KEYS = ("Glycine", "L-Arginine")
_GENES = ("b0001", "b0002")


def _panel(vals):
    """{cond: {gene: ratio}} from a flat 2x2 list, in _KEYS x _GENES order."""
    it = iter(vals)
    return {c: {g: next(it) for g in _GENES} for c in _KEYS}


class _Rxn:
    def __init__(self, rid):
        self.id = rid


class _FakeModel:
    """Minimal stand-in exposing the two attributes apply_nitrogen_condition touches."""

    def __init__(self, exchanges, medium=None):
        self.id = "fake"
        self.exchanges = [_Rxn(e) for e in exchanges]
        self.medium = dict(medium or {})


_EX = set(NITROGEN_EXCHANGES.values()) | {DEFAULT_CARBON, "EX_o2_e", "EX_pi_e"}


def test_curated_map_and_exclusions_are_disjoint_and_documented():
    assert not set(NITROGEN_EXCHANGES) & set(NITROGEN_UNMAPPABLE)
    # every exclusion carries a REASON, so the exclusion is auditable rather than silent
    assert all(isinstance(v, str) and v for v in NITROGEN_UNMAPPABLE.values())
    assert len(NITROGEN_EXCHANGES) == 13


def test_ammonium_is_closed_when_another_source_is_applied():
    """The load-bearing guard: a residual NH4 uptake would make every condition score as ammonium."""
    m = _FakeModel(_EX, medium={"EX_nh4_e": 10.0, DEFAULT_CARBON: 10.0})
    apply_nitrogen_condition(m, "EX_gly_e", all_nitrogen=tuple(NITROGEN_EXCHANGES.values()))
    assert "EX_nh4_e" not in m.medium
    assert m.medium["EX_gly_e"] == 10.0


def test_every_other_nitrogen_candidate_is_closed():
    m = _FakeModel(_EX, medium={e: 10.0 for e in NITROGEN_EXCHANGES.values()})
    apply_nitrogen_condition(m, "EX_arg__L_e", all_nitrogen=tuple(NITROGEN_EXCHANGES.values()))
    open_n = {e for e in NITROGEN_EXCHANGES.values() if e in m.medium}
    assert open_n == {"EX_arg__L_e"}


def test_carbon_source_is_held_open():
    """Nitrogen varies; carbon must NOT be removed or the condition is a starvation, not an N swap."""
    m = _FakeModel(_EX, medium={"EX_nh4_e": 10.0})
    apply_nitrogen_condition(m, "EX_ser__L_e")
    assert m.medium[DEFAULT_CARBON] == 10.0


def test_applying_ammonium_itself_leaves_it_open():
    m = _FakeModel(_EX, medium={})
    apply_nitrogen_condition(m, "EX_nh4_e", all_nitrogen=tuple(NITROGEN_EXCHANGES.values()))
    assert m.medium["EX_nh4_e"] == 10.0


def test_missing_exchange_raises_rather_than_silently_scoring():
    m = _FakeModel({DEFAULT_CARBON, "EX_nh4_e"})
    with pytest.raises(KeyError):
        apply_nitrogen_condition(m, "EX_ptrc_e")


def test_missing_carbon_exchange_raises():
    m = _FakeModel(set(NITROGEN_EXCHANGES.values()))
    with pytest.raises(KeyError):
        apply_nitrogen_condition(m, "EX_gly_e")


def test_determinism_gate_passes_on_pure_float_noise():
    """The real case: drifts ~1e-11 on ratios that sit nowhere near the 0.01 line."""
    a = _panel([1.0, 0.0, 0.9987, 0.0])
    b = _panel([1.0 - 3e-11, 0.0, 0.9987 + 2e-11, 0.0])
    v = determinism_verdict(a, b, _KEYS, _GENES, metric_a=0.6799, metric_b=0.6799)
    assert v["deterministic_at_claim_level"]
    assert v["n_call_flips"] == 0
    # and it passes with room to spare, not marginally -- the bar is not tuned to this panel
    assert v["safety_factor"] > 1e6


def test_determinism_gate_fails_on_a_call_flip():
    """A cell landing on opposite sides of FRAC changes a claim, so no amount of smallness saves it."""
    a = _panel([1.0, 0.00999, 1.0, 1.0])
    b = _panel([1.0, 0.01001, 1.0, 1.0])
    v = determinism_verdict(a, b, _KEYS, _GENES)
    assert not v["deterministic_at_claim_level"]
    assert v["n_call_flips"] == 1


def test_determinism_gate_fails_when_a_cell_sits_near_the_line():
    """The behaviour a FIXED tolerance cannot have, and the reason the margin is DERIVED.

    No call flips here and the drift is a tiny 1e-9 -- a fixed `max_delta < 1e-9`-style bar would pass
    this. But a cell sits 1e-7 from the decision line, so noise is only two orders of magnitude below
    the distance that would flip a call. The gate must fail: reproducibility of the CLAIM is not
    established when noise is that close to the boundary.
    """
    a = _panel([1.0, 0.0100001, 1.0, 0.0])
    b = _panel([1.0, 0.0100001 + 1e-9, 1.0, 0.0])
    v = determinism_verdict(a, b, _KEYS, _GENES)
    assert v["n_call_flips"] == 0, "precondition: this case must NOT be caught by the call-flip rule"
    assert v["max_abs_delta"] < 1e-8, "precondition: the drift is tiny, so only the margin can catch it"
    assert not v["deterministic_at_claim_level"]
    assert v["safety_factor"] < v["min_safety_factor"]


def test_determinism_gate_fails_when_the_headline_metric_moves():
    """The check that would have caught the retracted graded result: same-looking solves, moving number."""
    a = _panel([1.0, 0.0, 1.0, 0.0])
    v = determinism_verdict(a, a, _KEYS, _GENES, metric_a=0.6428, metric_b=0.6216)
    assert not v["deterministic_at_claim_level"]
    assert not v["headline_metric_matches"]


def test_redact_unverified_removes_every_solve_derived_claim():
    """Reproduces the ACTUAL defect: a `deterministic: false` artifact that still carried quotable numbers."""
    payload = {"deterministic": False, "per_cell_agreement": 0.6799, "best_constant_null": 0.5355,
               "per_condition": {"Glycine": {"agreement": 0.72}},
               "predictions": {"P1_bimodality": {"result": "REPLICATES"}},
               "verdict": "ALL_THREE_REPLICATE"}
    out = redact_unverified(payload, deterministic=False)
    for f in REDACTED_ON_NONDETERMINISM:
        assert out[f] is None, f"{f} must be REMOVED, not merely flagged"
    assert out["verdict"] == "NON_DETERMINISTIC_NO_VERDICT"
    assert "redaction_note" in out
    # label-derived numbers survive: they do not depend on any solve
    assert out["best_constant_null"] == 0.5355
    # the caller's dict is not mutated in place
    assert payload["per_cell_agreement"] == 0.6799


def test_redact_unverified_is_a_noop_when_the_gate_passes():
    """A control that always fires gets ignored, so it must be satisfiable."""
    payload = {"per_cell_agreement": 0.6799, "predictions": {"P1": 1}, "per_condition": {},
               "verdict": "ALL_THREE_REPLICATE"}
    assert redact_unverified(payload, deterministic=True) == payload


def test_nitrogen_conditions_intersects_assay_and_model():
    class _Conn:
        def execute(self, q, params):
            return [("Glycine",), ("L-Arginine",), ("casamino acids",), ("Gly-Glu",)]

    m = _FakeModel({"EX_gly_e", DEFAULT_CARBON})  # model lacks arginine exchange
    got = nitrogen_conditions(_Conn(), m)
    # present in assay AND mappable AND in the model -> glycine only
    assert got == {"Glycine": "EX_gly_e"}


# --------------------------------------------------------------- determinism gate (claim-level)
def test_identical_passes_are_deterministic():
    p = _panel([1.0, 0.0, 1.0, 0.0])
    v = determinism_verdict(p, p, _KEYS, _GENES, metric_a=0.68, metric_b=0.68)
    assert v["deterministic_at_claim_level"] is True
    assert v["n_call_flips"] == 0
    assert v["max_abs_delta"] == 0.0
    assert v["safety_factor"] == float("inf")


def test_float_noise_far_from_the_line_still_passes():
    """The whole point: ~1e-10 drift on cells sitting at 0.0/1.0 cannot change a call."""
    a = _panel([1.0, 0.0, 1.0, 0.0])
    b = _panel([1.0 + 3e-11, 0.0, 1.0 - 2e-11, 0.0])
    v = determinism_verdict(a, b, _KEYS, _GENES, metric_a=0.68, metric_b=0.68)
    assert v["deterministic_at_claim_level"] is True
    assert v["n_call_flips"] == 0


def test_determinism_gate_fails_when_a_cell_sits_near_the_line():
    """The anti-tuning guard the gate's docstring promises.

    Same tiny drift as the passing case, but one cell sits just above the 0.01 call line. The margin
    collapses, so the DERIVED safety factor drops below the bar and the gate fails on its own -- which a
    fixed tolerance could never do.
    """
    a = _panel([1.0, 0.0, 1.0, 0.010000001])
    b = _panel([1.0, 0.0, 1.0, 0.010000001 + 3e-11])
    v = determinism_verdict(a, b, _KEYS, _GENES, metric_a=0.68, metric_b=0.68)
    assert v["deterministic_at_claim_level"] is False
    assert v["safety_factor"] < v["min_safety_factor"]
    assert v["n_call_flips"] == 0, "fails on MARGIN, not on a flip -- that is the point"


def test_a_call_flip_fails_even_with_a_huge_margin_elsewhere():
    a = _panel([1.0, 0.0, 1.0, 0.005])
    b = _panel([1.0, 0.0, 1.0, 0.5])
    v = determinism_verdict(a, b, _KEYS, _GENES, metric_a=0.68, metric_b=0.68)
    assert v["deterministic_at_claim_level"] is False
    assert v["n_call_flips"] == 1
    assert v["call_flip_examples"]


def test_headline_metric_mismatch_fails_even_when_cells_agree():
    p = _panel([1.0, 0.0, 1.0, 0.0])
    v = determinism_verdict(p, p, _KEYS, _GENES, metric_a=0.6799, metric_b=0.6801)
    assert v["deterministic_at_claim_level"] is False
    assert v["headline_metric_matches"] is False


# --------------------------------------------------------------- redaction control
def test_redaction_removes_every_solve_derived_claim():
    payload = {"per_cell_agreement": 0.68, "per_condition": {"x": 1}, "predictions": {"P1": "REPLICATES"},
               "best_constant_null": 0.5355, "verdict": "ALL_THREE_REPLICATE"}
    out = redact_unverified(payload, deterministic=False)
    for f in REDACTED_ON_NONDETERMINISM:
        assert out[f] is None, f"{f} must not survive a determinism failure"
    assert out["verdict"] == "NON_DETERMINISTIC_NO_VERDICT"
    assert "redaction_note" in out
    # label-derived, NOT solve-derived -> must survive
    assert out["best_constant_null"] == 0.5355


def test_redaction_is_a_noop_when_deterministic():
    payload = {"per_cell_agreement": 0.68, "per_condition": {}, "predictions": {},
               "verdict": "ALL_THREE_REPLICATE"}
    assert redact_unverified(payload, deterministic=True) == payload
