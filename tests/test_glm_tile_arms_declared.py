"""The arms ledger is DURABLE STATE, and nothing used to prove the writes reached disk.

**The measured failure this file exists to stop (2026-10-14).** `ArmsLedger.declare_arm` and
`measure_arm` mutate in memory only; `save()` is a separate explicit call. Three ledger writes in one
session reported success and the file's mtime never moved, because every verification called
`ceiling_verdict()` on the *same in-memory object* that had just been mutated -- the check and the thing
checked were one object. The results were recoverable only because they also lived in committed wiki
artifacts.

So this file reads the JSON **from disk** and never constructs a ledger to ask it about itself.

**Why a repo test rather than a procedural rule:** an arm that silently fails to persist looks exactly
like an arm that was never declared, and R6's whole point is that a declared arm cannot evaporate -- so
a cheap arm passing can never cancel an expensive one that was never run. A rule in prose cannot catch
that; a test that reads the file can.

`arms_ledger` itself lives in the Soraya skills tree OUTSIDE this repo, so the one assertion that needs
it is guarded and skips with a reason rather than failing on a machine that has the repo but not the
skill. Everything else is computed from the JSON, which is committed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

LEDGER = Path(__file__).resolve().parents[1] / "soraya_runs" / "arms-ledger.json"

#: The exact criterion key. Spelled out rather than substring-matched: two criteria could share a
#: prefix, and this repo has a recorded trap where a lexical presence was a statement about the search
#: string rather than about the data.
CRITERION = ("ceiling: best achievable prediction of bacterial promoter activity from sequence, "
             "reported against the replicate-derived 0.941 ceiling")

#: The arm this plan declared, BEFORE running it (R6).
NEW_ARM = "promcalc-plus-conv-rank-average-hybrid-150bp"

#: Every arm that existed before this plan. Pinned so a later edit cannot quietly DROP one -- which is
#: the other half of R6: an arm dies only by its own executed evidence, never as a side effect.
PRE_EXISTING = (
    "conv-3width-13k-params-150bp",
    "fragment-rung-321k-rows-size-lever-LB-M9",
    "gc-single-feature",
    "multiresolution-pyramid-4-decades-UNRUN",
    "promcalc-full-decomposition-ridge-150bp",
    "promoter-calculator-biophysical-346params-fwd_max_tx",
    "promoteratlas-base-embeddings-frozen-150bp-padded-to-200",
    "pyramid-4scale-1.70decades-150bp",
    "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN",
)

SKILL_SCRIPTS = Path.home() / ".claude" / "skills" / "soraya" / "scripts"


def _arms() -> dict:
    """Every arm under the ceiling criterion, read FROM DISK."""
    assert LEDGER.is_file(), f"arms ledger missing at {LEDGER}"
    raw = json.loads(LEDGER.read_text(encoding="utf-8"))
    crits = raw.get("criteria", {})
    assert CRITERION in crits, (
        f"the ceiling criterion is absent from the ledger. Present keys: {sorted(crits)}"
    )
    return crits[CRITERION]


def test_the_hybrid_arm_is_on_disk_not_merely_in_memory():
    """The declaration survived the process that made it."""
    arms = _arms()
    assert NEW_ARM in arms, (
        f"{NEW_ARM!r} is not in the ledger on disk. If it was just declared, `save()` was probably "
        f"not called -- declare_arm mutates in memory only. Present: {sorted(arms)}"
    )


def test_the_hybrid_arm_status_is_declared_or_measured():
    """`planned` before the run, `measured` after it -- both are valid; anything else is not.

    Deliberately accepts BOTH so this test survives the arm actually being measured. A test that
    pinned `planned` would have to be edited the moment the plan succeeded, and a test edited to make
    it pass is not a guard.
    """
    arm = _arms()[NEW_ARM]
    assert arm.get("status") in {"planned", "measured"}, (
        f"unexpected status {arm.get('status')!r} -- a killed or malformed arm here would mean the "
        f"R6 no-evaporation property was bypassed"
    )


def test_no_pre_existing_arm_was_dropped():
    """Adding an arm must not remove one. R6: an arm dies only by its own executed evidence."""
    arms = _arms()
    missing = [n for n in PRE_EXISTING if n not in arms]
    assert not missing, f"arms disappeared from the ledger: {missing}"


def test_a_planned_arm_blocks_the_criterion_in_the_json_itself():
    """Computed from the committed JSON, with no ledger object involved.

    This is the property R5 rests on: while any arm is still `planned` the ceiling question is OPEN,
    so a cheap arm clearing a threshold cannot settle it. Asserting it from the raw file means the
    check holds even where the Soraya skill is not installed.
    """
    arms = _arms()
    planned = [n for n, a in arms.items() if a.get("status") == "planned"]
    assert planned, (
        "no arm is planned, so this test can no longer observe the blocking property. If every arm "
        "really is measured, replace this with an assertion about the settled verdict instead of "
        "deleting it -- a test whose subject vanished should be rewritten, not bumped."
    )


def test_ceiling_verdict_agrees_with_the_json_while_an_arm_is_untried():
    """The same property through the real `ceiling_verdict`, when the skill is available.

    Guarded + self-retiring: it runs wherever the Soraya scripts resolve and skips with a reason
    elsewhere, rather than coupling this repo's suite to a path outside it.
    """
    if not (SKILL_SCRIPTS / "arms_ledger.py").is_file():
        pytest.skip(f"arms_ledger.py not resolvable at {SKILL_SCRIPTS}; JSON-level check covers it")
    sys.path.insert(0, str(SKILL_SCRIPTS))
    try:
        from arms_ledger import ArmsLedger
    finally:
        sys.path.pop(0)
    verdict = ArmsLedger(str(LEDGER.parents[1])).ceiling_verdict(CRITERION).verdict
    assert verdict.startswith("blocked:"), (
        f"a planned arm is present in the JSON but ceiling_verdict says {verdict!r}; the ledger and "
        f"its verdict function disagree about whether the ceiling question is settled"
    )
