# -*- coding: utf-8 -*-
"""Amendment 5 of the tile-chain substrate wall memo: the promcalc half is closed.

A separate file from `test_tile_chain_substrate_wall.py` on purpose -- that file pins
amendments 1-3 and is already long. Amendment 4 turned out NOT to be pinned anywhere
(the saved plan's test strategy says "as amendment 4 is pinned", which is false), so the
both-or-neither guard below covers 4 and 5 together rather than leaving the older one
unguarded while adding a guard for the newer.

The load-bearing test here is the LAST one. Amendment 5 says the substrate blocker is
gone, and the single most expensive way to misread that is as "the ceiling arm is
settled". R6 says an arm dies only by its own executed evidence; this pins that the
sibling NO_GAIN verdict did not retire it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MEMO = REPO / "wiki" / "tile_chain_substrate_wall_2026-10-21.md"
SIDECAR = REPO / "wiki" / "tile_chain_substrate_wall_2026-10-21.json"
ARMS = REPO / "soraya_runs" / "arms-ledger.json"

CEILING = ("ceiling: best achievable prediction of bacterial promoter activity from "
           "sequence, reported against the replicate-derived 0.941 ceiling")
ARM = "orthogonal-modality-stack-promcalc-plus-N-channels"


def _side() -> dict:
    return json.loads(SIDECAR.read_text(encoding="utf-8"))


def _a5() -> dict:
    return _side()["amendment_5"]


def test_amendment_5_pins_a_wall_class_from_the_memos_OWN_closed_vocabulary():
    """Amendment 3's first draft invented a fourth spelling. The vocabulary is the guard.

    Read the allowed set from the sidecar rather than restating it here: a test that
    hardcodes the vocabulary cannot notice the vocabulary changing underneath it."""
    side = _side()
    vocab = side["wall_class_vocabulary"]
    assert _a5()["wall_class_held"] in vocab, (
        f"amendment 5 names a class outside {vocab}")
    assert _a5()["wall_class_vocabulary_checked"] is True


def test_the_closure_route_is_named_and_explains_why_it_was_never_external_wall():
    """`resolved: true` alone is not a record. A future reader needs to know HOW, because
    the route is what distinguishes a wall this host could close from one it could not."""
    a5 = _a5()
    assert a5["resolved"] is True
    assert "external-wall" in a5["closure_route"], (
        "say why this was never the class that needs someone else")


def test_the_FROZEN_dependencies_array_is_declared_untouched_and_IS_untouched():
    """Amendments 3 and 4 amend; they do not rewrite the frozen record. Check the claim
    AND the fact -- a self-report that the array was not edited is not evidence."""
    side = _side()
    assert _a5()["frozen_dependencies_array_edited"] is False
    deps = side["dependencies"]
    assert len(deps) == 2
    assert {d["wall_class"] for d in deps} == {"code-closable", "data-absent"}, (
        "the frozen array still records the classes as they were at the time of writing")


def test_both_the_MEMO_PROSE_and_the_SIDECAR_carry_amendments_4_and_5():
    """The sidecar is machine-readable; the memo is what a human opens. Both or neither."""
    md = MEMO.read_text(encoding="utf-8")
    side = _side()
    for n in (4, 5):
        assert f"AMENDMENT {n}" in md, f"memo prose is missing amendment {n}"
        assert f"amendment_{n}" in side, f"sidecar is missing amendment_{n}"
    assert md.count("## AMENDMENT 5") == 1, "appended twice"


def test_risk_flag_1s_mitigation_is_recorded_UNMET_rather_than_quietly_dropped():
    """The saved plan promised "Step 6 records the version so this risk cannot recur".
    It does not. A plan's mitigation silently not happening is the failure mode a wall
    memo exists to catch, so the memo must say UNMET in both surfaces."""
    r = _a5()["risk_flag_1_mitigation"]
    assert r["status"] == "UNMET"
    assert "Step 6 records the version" in r["plan_claim"], (
        "quote the claim being contradicted, or the contradiction is not checkable")
    assert "UNMET" in MEMO.read_text(encoding="utf-8")


def test_the_unrecoverability_rests_on_several_INDEPENDENT_checks_not_one():
    """One negative lookup is a lookup; four agreeing negatives plus an untouched sklearn
    are a finding. The sklearn check is the one that is positive evidence rather than
    absence -- it confirms amendment 2's refutation held."""
    c = _a5()["risk_flag_1_mitigation"]["checks_run"]
    assert len(c) >= 4
    assert c["grep_count_in_uv_lock"] == 0
    assert c["sklearn_version_now"] == "1.8.0", (
        "sklearn being untouched is what rules out the refuted 1.9.1 story")


def test_the_uv_run_mechanism_is_labelled_a_LEADING_account_not_a_finding():
    """I did not watch uv remove the package. This project prices asserted-then-refuted
    causes at three consecutive instances in one session, so the hedge is the point."""
    m = _a5()["risk_flag_1_mitigation"]["leading_mechanism"]
    assert "NOT watched" in m
    assert "inferred" in m.lower()


def test_provenance_by_reproduction_records_what_it_does_NOT_support():
    """Exact agreement with the original run is strong and is NOT a version pin. Without
    this clause the next session reads a closed debt where an open one exists."""
    p = _a5()["provenance_by_reproduction"]
    assert "what_this_does_not_support" in p, (
        "the does-not-support clause must be its OWN key, not buried in prose")
    assert "pinned" in p["what_this_does_not_support"]
    assert "OPEN" in p["what_this_does_not_support"]
    for anchor in ("gc_seed0", "max_tx_median"):
        a = p[anchor]
        assert a["exact"] is True
        assert abs(a["measured"] - a["expected"]) <= a["tolerance"]


def test_the_reconcile_outcome_is_declared_absent_as_a_FIELD_and_re_derivable():
    """A reader who goes looking for a `reconcile: passed` key will not find one. Saying
    so, and saying where the evidence IS, is cheaper than them concluding it never ran."""
    p = _a5()["provenance_by_reproduction"]["reconcile_outcome_is_not_a_field"]
    assert "INDETERMINATE_RECONCILE_FAILED" in p
    assert "promcalc_decomposition_arm.py" in p


def test_the_naming_defect_is_recorded_as_found_NOT_fixed_with_its_verdict_impact():
    """`control_max_abs` tests a median. The defect matters to a reader, not to the
    verdict, and conflating those two would either fake a problem or hide one."""
    d = _a5()["naming_defect_found_not_fixed"]
    assert d["verdict_impact"].startswith("none")
    assert d["median_measured"] == pytest.approx(0.0342)
    assert set(d["seeds_exceeding_bar_in_abs_value"]) == {"seed_5", "seed_7"}
    assert all(abs(v) > 0.05 for v in d["seeds_exceeding_bar_in_abs_value"].values()), (
        "the seeds named must actually exceed the bar, or the defect is misdescribed")
    assert "frozen PREREGISTERED" in d["why_not_renamed"]


def test_the_verdict_carries_the_clause_that_stops_it_being_over_read():
    """NO_GAIN is about the DECOMPOSITION vs the scalar reduced from it -- not about
    promoter thermodynamics, which beats GC 10/10 on both arms. The compressed reading
    is the error this project prices highest."""
    v = _a5()["verdict"]
    assert v["value"] == "NO_GAIN"
    assert v["includes_zero"] is True and v["wins"] == 2 and v["n"] == 10
    assert v["paired_median"] < 0
    dn = v["do_not_over_read"]
    assert "does NOT say" in dn
    assert "+0.0802" in dn and "+0.0865" in dn, "carry both GC margins or the rebuttal is bare"


def test_the_arm_RESOLVED_BY_ITS_OWN_RUN_and_the_ceiling_is_STILL_blocked():
    """THE load-bearing test, updated 2026-10-26 when the arm was actually run.

    It previously pinned `status == "planned"` / `metric is None` / `evidence == ""`, which was the
    correct guard while the arm was UNRUN: it stopped a SIBLING question's NO_GAIN verdict from
    retiring an arm that had measured nothing (R6 -- an arm dies only by its OWN executed evidence,
    the 2026-07-09 defect this rail exists to prevent).

    The arm has now been run (`wiki/glm_tile_stack_arm_2026-10-26.md`), so those three assertions
    would pin the OPPOSITE failure -- hiding executed evidence. The rail is kept in a stronger form:
    the arm may only leave `planned` by naming its own artifact, and `blocking()` must STILL name the
    criterion, now through the remaining `planned` siblings.
    """
    arm = json.loads(ARMS.read_text(encoding="utf-8"))["criteria"][CEILING][ARM]

    # R6, enforced at the seam that matters: it left `planned` only with its OWN evidence
    assert arm["status"] == "measured"
    assert arm["metric"] is not None
    assert "glm_tile_stack_arm_2026-10-26" in arm["evidence"], (
        "the arm resolved without naming the artifact that resolved it")
    assert arm["resolved_at"], "a resolved arm with no resolved_at is not auditable"

    # the verdict is recorded honestly -- it did NOT clear the bar
    v = arm["verdict_2026_10_26"]
    assert "WEAK_STACK" in v
    assert "did not clear the bar" in v

    # the 2026-10-23 note is SUPERSEDED, not rewritten -- the memo's own convention
    note = arm["substrate_status_2026_10_23"]
    assert "RUNNABLE, NOT RUN" in note, "the superseded note was edited rather than left in place"

    # the untested form is recorded as a GAP, and deliberately not declared as a new arm
    assert "NEW modality family" in arm["narrowed_successor_NOT_declared_as_an_arm"]

    # and the ceiling is still blocked -- by siblings, verified rather than asserted
    crit = json.loads(ARMS.read_text(encoding="utf-8"))["criteria"][CEILING]
    planned = [k for k, a in crit.items() if a.get("status") == "planned"]
    assert planned, "no `planned` arm remains, so the ceiling would read as settled"
    assert ARM not in planned

    sys.path.insert(0, str(Path.home() / ".claude" / "skills" / "soraya" / "scripts"))
    try:
        import arms_ledger
    except ImportError:  # pragma: no cover - the skill is not installed in CI
        pytest.skip("soraya arms_ledger not importable here")
    blocked = arms_ledger.ArmsLedger(REPO).blocking()
    assert any(CEILING in str(b) for b in blocked), (
        f"blocking() no longer names the ceiling criterion; got {blocked!r}")


def test_the_raw_JSON_edit_PRESERVED_every_unmodelled_key_on_the_arm():
    """`Arm.to_dict()` models five fields; this arm carries eleven. A `ArmsLedger.save()`
    round-trip would have silently dropped six, so the Step-7 edit was raw JSON -- and that
    choice is only safe if the preservation is checked."""
    arm = json.loads(ARMS.read_text(encoding="utf-8"))["criteria"][CEILING][ARM]
    for key in ("name", "cost_tier", "status", "metric", "evidence", "declared_at", "resolved_at",
                "rationale", "kill_option_withdrawn", "verdict_consequence",
                "substrate_status_2026_10_23"):
        assert key in arm, f"raw-JSON edit dropped {key}"
    # the three fields a to_dict() round-trip would NOT model, still carrying their text
    assert "WITHDRAWN" in arm["kill_option_withdrawn"]
    assert "blocked:untried-arm" in arm["verdict_consequence"]


def test_superseded_statements_are_RECORDED_not_rewritten():
    """The memo still says "No arm was run"; an arm has since run. Amendment 1 established
    the convention for exactly this -- supersede in the amendment, leave the old text so the
    record of what was believed when stays readable. Amendment 5 follows it rather than
    quietly editing a claim that was true on 2026-10-21."""
    sup = _a5()["statements_above_now_false_left_in_place"]
    assert len(sup) >= 2
    joined = " ".join(sup)
    assert "No arm was run" in joined, "name the superseded claim verbatim or it is not checkable"
    assert "NOT done and was not started" in joined, "amendment 4's claim is superseded too"
    # The old text must still BE there -- the convention is supersede, not delete.
    assert "not a result" in MEMO.read_text(encoding="utf-8").lower()
