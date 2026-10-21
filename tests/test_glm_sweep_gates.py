"""The post-run grid gates.

WHY (2026-10-20/21): the sweep plan's Step 4 excludes any point whose `knob_verified_end_to_end`
is false. Verification is gated on a non-default knob having been REQUESTED, so at T=1.0 none
was, and the real anchor artifact records all three knob fields false. A literal reading of that
rule therefore excludes the ANCHOR -- the comparator the sharpened points are measured against.

A default point is not UNVERIFIED; it is NOTHING-TO-VERIFY. These tests pin that distinction
and pin the cross-point null invariant beside it.
"""
from __future__ import annotations

import json
import pathlib

import pytest

from dna_decode.glm.sweep_gates import (
    ELIGIBLE,
    ELIGIBLE_DEFAULT,
    ELIGIBLE_VERIFIED,
    EXCLUDED_KNOB_DEAD,
    EXCLUDED_UNVERIFIED,
    GridGateError,
    assemble_grid,
    null_reproduced,
    verdict_eligibility,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent

_DEFAULT = {"requested_non_default": False, "knob_wired": False, "knob_verified_end_to_end": False}
_VERIFIED = {"requested_non_default": True, "knob_wired": True, "knob_verified_end_to_end": True}
_UNVERIFIED = {"requested_non_default": True, "knob_wired": True, "knob_verified_end_to_end": False}
_DEAD = {"requested_non_default": True, "knob_wired": False, "knob_verified_end_to_end": False}


# ------------------------------------------------- the defect this module exists to prevent
def test_the_DEFAULT_point_is_ELIGIBLE_not_excluded():
    """The whole point. A bool-returning rule would have excluded the comparator."""
    assert verdict_eligibility(_DEFAULT) == ELIGIBLE_DEFAULT
    assert verdict_eligibility(_DEFAULT) in ELIGIBLE


def test_the_REAL_anchor_artifact_is_ELIGIBLE():
    """Pinned against the artifact on disk, not a fixture written to agree with the rule."""
    art = ROOT / "wiki" / "glm_sweep_anchor_2026-10-20.json"
    if not art.is_file():
        pytest.skip("anchor artifact absent")
    rec = json.loads(art.read_text(encoding="utf-8")).get("knob_contract") or {}
    # all three fields are false in the real file -- exactly the shape that reads as "unverified"
    assert rec.get("knob_verified_end_to_end") is False
    assert verdict_eligibility(rec) == ELIGIBLE_DEFAULT


def test_a_naive_boolean_rule_WOULD_have_excluded_it():
    """Shows the counterfactual rather than asserting it in prose.

    If this ever fails, the artifact's shape changed and the module's reason for existing
    should be re-read before the rule is simplified.
    """
    naive = bool(_DEFAULT["knob_verified_end_to_end"])
    assert naive is False                                   # the naive rule excludes
    assert verdict_eligibility(_DEFAULT) in ELIGIBLE        # the real rule does not


# ------------------------------------------------------------------ the full ladder
def test_a_verified_nondefault_point_is_eligible():
    assert verdict_eligibility(_VERIFIED) == ELIGIBLE_VERIFIED


def test_a_requested_but_unverified_point_is_EXCLUDED():
    assert verdict_eligibility(_UNVERIFIED) == EXCLUDED_UNVERIFIED


def test_a_DEAD_knob_is_excluded_as_dead_not_as_unverified():
    """Order matters: a knob that installed no warper is inert whatever verification says.

    This is the 2026-10-20 regression's own signature, so it gets its own reason string
    rather than being folded into the unverified bucket.
    """
    assert verdict_eligibility(_DEAD) == EXCLUDED_KNOB_DEAD


def test_all_four_outcomes_are_reachable_and_distinct():
    got = {verdict_eligibility(r) for r in (_DEFAULT, _VERIFIED, _UNVERIFIED, _DEAD)}
    assert got == {ELIGIBLE_DEFAULT, ELIGIBLE_VERIFIED, EXCLUDED_UNVERIFIED, EXCLUDED_KNOB_DEAD}


def test_an_OLD_artifact_without_the_stamp_is_eligible_as_default():
    """A pre-stamp artifact records default sampling; that is case 1, honestly eligible."""
    assert verdict_eligibility({}) == ELIGIBLE_DEFAULT


def test_a_non_dict_is_REFUSED_rather_than_coerced():
    with pytest.raises(GridGateError, match="must be a dict"):
        verdict_eligibility(None)


# ------------------------------------------------------------------ the null invariant
def test_an_identical_null_reproduces():
    r = null_reproduced(0.5990, 0.5990, 0.5990)
    assert r["ok"] and r["reason"] == "null_reproduced"


def test_a_null_that_MOVED_vs_the_anchor_fails():
    r = null_reproduced(0.7475, 0.5990, 0.5990)
    assert not r["ok"] and r["reason"] == "null_moved_vs_anchor"


def test_a_null_matching_the_anchor_but_NOT_the_record_fails_separately():
    """Two different failures, distinguishable: within-sweep drift vs drift from the record."""
    r = null_reproduced(0.7475, 0.7475, 0.5990)
    assert not r["ok"] and r["reason"] == "null_moved_vs_recorded"


def test_a_MISSING_null_is_reported_not_treated_as_a_match():
    assert null_reproduced(None, 0.5990)["reason"] == "null_missing"
    assert null_reproduced(0.5990, None)["reason"] == "null_missing"


def test_tol_defaults_to_ZERO_so_strictness_cannot_be_lost_by_omission():
    assert null_reproduced(0.5991, 0.5990)["ok"] is False
    assert null_reproduced(0.5991, 0.5990, tol=0.001)["ok"] is True


# ------------------------------------------------------------------ assembly
def _grid():
    return [
        {"label": "T=1.0", "knob_contract": _DEFAULT, "null": 0.5990},
        {"label": "T=0.8", "knob_contract": _VERIFIED, "null": 0.5990},
        {"label": "T=0.6", "knob_contract": _VERIFIED, "null": 0.5990},
    ]


def test_assemble_grid_keeps_the_anchor_AND_the_sharpened_points():
    g = assemble_grid(_grid(), recorded_null=0.5990)
    assert [e["label"] for e in g["eligible"]] == ["T=1.0", "T=0.8", "T=0.6"]
    assert g["excluded"] == []
    assert g["anchor_null"] == 0.5990        # inferred from the default-sampling point


def test_assemble_grid_EXCLUDES_an_unverified_point_WITH_a_reason():
    pts = _grid()
    pts[1]["knob_contract"] = _UNVERIFIED
    g = assemble_grid(pts, recorded_null=0.5990)
    assert [e["label"] for e in g["eligible"]] == ["T=1.0", "T=0.6"]
    assert g["excluded"] == [{"label": "T=0.8", "reason": EXCLUDED_UNVERIFIED}]


def test_assemble_grid_EXCLUDES_a_point_whose_NULL_MOVED():
    pts = _grid()
    pts[2]["null"] = 0.7475
    g = assemble_grid(pts, recorded_null=0.5990)
    assert [e["label"] for e in g["eligible"]] == ["T=1.0", "T=0.8"]
    bad = g["excluded"][0]
    assert bad["label"] == "T=0.6"
    assert bad["reason"] == "excluded:null_moved_vs_anchor"
    assert bad["null_check"]["delta_vs_anchor"] > 0


def test_NO_point_is_ever_dropped_silently():
    """An unexplained absence reads as 'not run', which is the wrong story."""
    pts = _grid()
    pts[1]["knob_contract"] = _DEAD
    pts[2]["null"] = 0.7475
    g = assemble_grid(pts, recorded_null=0.5990)
    assert len(g["eligible"]) + len(g["excluded"]) == g["n_points"] == 3
    assert all(e.get("reason") for e in g["excluded"])
def test_the_plan_prose_cites_THIS_function_as_the_single_definition():
    """Step 4 of the sweep plan aligns its wording with this module; pin the citation.

    The rule previously lived in three places and disagreed with itself, so the prose must
    point at the code rather than restating the rule a fourth time.
    """
    plan = (ROOT / "plans" / "Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan"
            / "technical-plan.md")
    if not plan.is_file():
        pytest.skip("sweep plan not at this path (archived?)")
    t = plan.read_text(encoding="utf-8")
    assert "verdict_eligibility" in t, "the plan must cite the code-owned rule"

    # The real pin is the ABSENCE of the unscoped sentence. Measured: before the amendment
    # `verdict_eligibility` appeared 0 times (so that assertion bites) while "non-default"
    # already appeared once at gate 2 -- so asserting its PRESENCE was vacuous, and is dropped.
    flat = " ".join(t.split())
    assert "a point whose knob is unverified is excluded from the verdict" not in flat, (
        "the unscoped wording is back; it excludes the T=1.0 anchor, which is the comparator")
    assert "requested a** non-default knob" in flat or "requested a non-default knob" in flat \
        or "**requested a non-default knob**" in flat, (
        "the exclusion must be scoped to points that REQUESTED a non-default knob")
