"""Guards for two KNOWN arms-ledger hazards that nothing could previously check.

Neither hazard is new. Both were recorded in the ledger itself and left deliberately unfixed; what was
missing is a check, so a future edit could re-introduce them silently.

HAZARD 1 -- A METRIC THAT IS NOT IN THE CRITERION'S UNIT.
`ceiling_verdict` computes `best` as a raw `max(measured, key=metric)` with no unit check, and the
criterion is measured in bacterial tile Spearman. Two arms under it are NOT in that unit:
  * `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` -- 0.4853 is a HUMAN lentiMPRA
    Spearman (different assay, different noise floor), and its own long-standing
    `metric_not_commensurable_with_siblings` note says a future human arm above 0.5076 "WOULD silently
    become `best`".
  * `fragment-rung-321k-rows-size-lever-LB-M9` -- 0.185 is not a correlation at all. It is a FRACTION
    of the 0.7928 FRAGMENT ceiling, while the criterion quotes the 0.941 TILE ceiling.
**This file does NOT change `best` semantics and does NOT touch the shared `arms_ledger` module.**
Restructuring the guard was explicitly surfaced to the user and deliberately parked; overriding that
here would be taking an authority call. The test closes only the named hazard: such an arm must never
BE `best`. Today it cannot be (0.4853 and 0.185 both sit below the bacterial best), so this is a
tripwire for the future, and the non-vacuity test below proves it would actually fire.

HAZARD 2 -- A NAME THAT CONTRADICTS ITS OWN STATUS.
Two arm names end in `-UNRUN`; one of them has `status: measured`, `metric` 0.4853 and evidence opening
"RAN 2026-11-01". **Renaming was considered and withdrawn:** the name is the ledger KEY and 21 files
reference these two names, three of them test files, so renaming rewrites history and breaks tests to
fix a cosmetic read. The fix is that the row DECLARES the contradiction.

Why a required field on the offending row rather than an allowlist here: this repo has a measured
lesson that a hand-enumerated exclusion list beside the data drifts and fails by PASSING. A required
acknowledgement ON the row inverts that -- the data carries its own declaration, and a NEW offender
fails until someone writes one. There is no list in this file to go stale.

Pure structural checks against the live ledger: no network, no torch, no Kaggle, no D: drive.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")

import arms_ledger as al  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
LEDGER = REPO / "soraya_runs" / "arms-ledger.json"
CRITERION_PREFIX = "ceiling: best achievable prediction of bacterial promoter"

NON_COMMENSURABLE_KEY = "metric_not_commensurable_with_siblings"
STALE_NAME_KEY = "name_is_a_stale_declaration_label"


def _live_ledger() -> al.ArmsLedger:
    return al.ArmsLedger(cwd=str(REPO))


def _criterion(ledger: al.ArmsLedger) -> str:
    return next(k for k in ledger._data["criteria"] if k.startswith(CRITERION_PREFIX))


@pytest.fixture
def sandbox(tmp_path: Path) -> Path:
    """A throwaway copy -- the mutating non-vacuity tests must never touch the real ledger."""
    dest = tmp_path / "soraya_runs"
    dest.mkdir(parents=True)
    shutil.copy2(LEDGER, dest / "arms-ledger.json")
    return tmp_path


# ---------------------------------------------------------------------------------------------------
# the two checks, as ONE function each, so the live test and its non-vacuity twin share an implementation
# ---------------------------------------------------------------------------------------------------
def _best_is_commensurable(ledger: al.ArmsLedger) -> tuple[bool, str]:
    crit = _criterion(ledger)
    verdict = ledger.ceiling_verdict(crit)
    best = verdict.best
    if best is None:
        return True, "no best arm (nothing measured)"
    arm = ledger._data["criteria"][crit][best]
    if NON_COMMENSURABLE_KEY in arm:
        return False, (
            "best=%r carries %s, so `best` is a max across different units"
            % (best, NON_COMMENSURABLE_KEY))
    return True, "best=%r is in the criterion's unit" % best


def _unacknowledged_stale_names(ledger: al.ArmsLedger) -> list[str]:
    """An arm whose NAME asserts UNRUN while its STATUS says it ran must declare the contradiction."""
    crit = _criterion(ledger)
    bad = []
    for name, arm in ledger._data["criteria"][crit].items():
        if not name.endswith("-UNRUN"):
            continue
        if STALE_NAME_KEY not in arm:
            bad.append(name)
    return bad


# ---------------------------------------------------------------------------------------------------
# hazard 1
# ---------------------------------------------------------------------------------------------------
def test_a_non_commensurable_arm_is_never_best():
    ok, why = _best_is_commensurable(_live_ledger())
    assert ok, why


def test_the_non_commensurable_marker_is_actually_USED():
    """Non-vacuity for the guard above: if nothing were marked, it could never fail."""
    ledger = _live_ledger()
    crit = _criterion(ledger)
    marked = [n for n, a in ledger._data["criteria"][crit].items() if NON_COMMENSURABLE_KEY in a]
    assert len(marked) >= 2, (
        "expected at least the human-panel arm and the fragment-rung arm to be marked "
        "non-commensurable; got %r" % marked)
    for n in marked:
        note = ledger._data["criteria"][crit][n][NON_COMMENSURABLE_KEY]
        assert len(note) > 120, "%s's marker must say WHY, not just flag" % n


def test_the_commensurability_guard_FIRES_when_such_an_arm_becomes_best(sandbox: Path):
    """Proves the tripwire is live rather than merely true-by-luck.

    Today both non-commensurable arms sit below the bacterial best, so the guard passes for free. Here
    one is pushed above it in a THROWAWAY copy; `best` must move to it and the guard must then fail.
    """
    path = sandbox / "soraya_runs" / "arms-ledger.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    crit = next(k for k in data["criteria"] if k.startswith(CRITERION_PREFIX))
    target = "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN"
    assert NON_COMMENSURABLE_KEY in data["criteria"][crit][target]
    data["criteria"][crit][target]["metric"] = 0.99
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")

    ledger = al.ArmsLedger(cwd=str(sandbox))
    assert ledger.ceiling_verdict(crit).best == target, "the mutation did not move `best`"
    ok, why = _best_is_commensurable(ledger)
    assert not ok, "the guard did NOT fire on a non-commensurable best: %s" % why
    assert NON_COMMENSURABLE_KEY in why


# ---------------------------------------------------------------------------------------------------
# hazard 2
# ---------------------------------------------------------------------------------------------------
def test_every_UNRUN_named_arm_declares_that_its_name_is_stale():
    bad = _unacknowledged_stale_names(_live_ledger())
    assert not bad, (
        "these arms carry `-UNRUN` in their name without a %s field declaring it: %r. The name is the "
        "ledger key and is referenced across the repo, so renaming is not the fix -- declare it."
        % (STALE_NAME_KEY, bad))


def test_the_stale_name_check_is_non_vacuous(sandbox: Path):
    """Remove one acknowledgement and the check must name that arm."""
    path = sandbox / "soraya_runs" / "arms-ledger.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    crit = next(k for k in data["criteria"] if k.startswith(CRITERION_PREFIX))
    target = "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN"
    del data["criteria"][crit][target][STALE_NAME_KEY]
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")

    bad = _unacknowledged_stale_names(al.ArmsLedger(cwd=str(sandbox)))
    assert bad == [target], "expected the check to name %r, got %r" % (target, bad)


def test_the_measured_UNRUN_arm_really_is_the_contradiction_this_guards():
    """Anchors the guard to the real case, so it is not an abstract rule about suffixes.

    Reads status from the LIVE ledger rather than pinning a metric -- this repo converted its arms
    guards from exact-state pins to invariants, and a pinned number would fail on any legitimate
    re-measurement.
    """
    ledger = _live_ledger()
    crit = _criterion(ledger)
    arm = ledger._data["criteria"][crit]["substrate-change-constructed-eukaryotic-or-human-panel-UNRUN"]
    assert arm["status"] == "measured"
    assert isinstance(arm["metric"], float)
    assert arm["evidence"].startswith("RAN "), (
        "the contradiction this guard exists for is a name saying UNRUN on an arm whose evidence "
        "opens by saying it ran")
