"""Guard: a SUPERSEDED plan may not keep a resumable execution-state file.

WHY THIS IS A TEST AND NOT A THIRD HAND-FIX. The hazard has now occurred twice. A plan gets
partially executed, a successor plan takes over its remaining steps, and its
`.claude/execute-plan-state/<name>.json` is left behind still listing those steps as `pending` --
so `/execute-plan` on the old path offers to resume a design that measurement has already
refuted. It was caught by hand for
`Declare_And_Run_All_Three_Successor_Arms_…` and again for
`Finish_The_Three_Successor_Arms_…`. Catching it a third time by hand would be the wrong
response.

THE REMEDY IS A RENAME, NOT A DELETE, AND THE GUARD ENFORCES THAT DISTINCTION. Renaming to
`.json.superseded` keeps the record of which steps DID execute; deleting would destroy it. So the
guard additionally requires that a superseded plan's record still be readable, which means it
cannot be satisfied by `rm`.

SCOPE, STATED BECAUSE IT IS NARROWER THAN IT FIRST LOOKS. The guard enforces ONE invariant --
`superseded` status versus a live state file. It deliberately does NOT require every entry to
carry a `**Status:**` field, because that was MEASURED to fail on 39 pre-existing entries, and
satisfying it would mean inventing a status for each from its state file without reading the plan.
Instead the coverage gap is RATCHETED: the existing 39 are grandfathered and a NEW statusless
entry fails. A ratchet that cannot be met is a ratchet that gets deleted.

`as-built` is deliberately COMPATIBLE with a live state file, and that is the load-bearing
distinction. `as-built` means the plan ran and execution diverged -- its state file is the record
of what happened and should stay. `superseded` means another plan took its remaining steps over,
which is what makes resuming it actively wrong.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
INDEX = REPO / "wiki/plans-index.md"
STATE_DIR = REPO / ".claude/execute-plan-state"

#: Statuses that make resuming a plan actively wrong.
INCOMPATIBLE_WITH_A_LIVE_STATE_FILE = ("superseded",)

#: Statuses that legitimately coexist with a live state file.
COMPATIBLE_WITH_A_LIVE_STATE_FILE = ("candidate", "active", "as-built", "executed")

#: MEASURED 2026-11-03, after this plan's two annotations: 83 entries, 39 without a Status field.
#: Grandfathered, not endorsed. A NEW statusless entry pushes the count to 40 and fails.
MAX_ENTRIES_WITHOUT_STATUS = 39


# --------------------------------------------------------------------- pure helpers


def index_statuses(text: str) -> dict:
    """`{plan_key: status_or_None}` for every `## [plan_file: …]` entry.

    Pure, so the invariant can be driven on fixtures as well as on the real index. Real entries
    write both `**Status:** superseded` and `**Status:** as-built — long prose…`, so the status is
    the leading run of letters and HYPHENS, which stops at the space before the em-dash.

    The hyphen must be INSIDE the class. The first version excluded it to stop at a separator and
    truncated `as-built` to `as`, which silently dropped `as-built` from the compatible set — a
    latent bug, because `superseded` has no hyphen and looked correct throughout.
    """
    out: dict = {}
    for block in re.split(r"(?m)^## \[plan_file: ", text)[1:]:
        key = block.split("]", 1)[0].strip()
        m = re.search(r"(?m)^\*\*Status:\*\*\s*([A-Za-z][A-Za-z-]*)", block)
        out[key] = m.group(1).strip() if m else None
    return out


def live_state_plans(names) -> set:
    """Plan names with a RESUMABLE state file. `.json.superseded` is not resumable."""
    return {n[: -len(".json")] for n in names if n.endswith(".json")}


def violations(statuses: dict, live: set) -> list:
    """Plans whose status forbids resuming while a resumable state file still exists."""
    bad = []
    for key, status in sorted(statuses.items()):
        if status in INCOMPATIBLE_WITH_A_LIVE_STATE_FILE and key.rstrip("/") in live:
            bad.append(key.rstrip("/"))
    return bad


# --------------------------------------------------------------------- non-vacuity


def test_the_checker_FAILS_on_the_pre_fix_condition():
    """Non-vacuity, by constructing the exact state the tree was in before this plan's Step 1.

    Without this, the guard could be trivially satisfied by a parser that finds nothing, and would
    pass on a tree that still carries the hazard.
    """
    text = ("## [plan_file: Some_Plan/] 2026-11-02\n"
            "**Status:** superseded\n"
            "**Summary:** x\n\n---\n")
    live = live_state_plans(["Some_Plan.json"])
    assert violations(index_statuses(text), live) == ["Some_Plan"]


def test_the_checker_PASSES_once_the_state_file_is_renamed():
    """The other half of the proof: the SAME index text passes once the file is neutralised."""
    text = ("## [plan_file: Some_Plan/] 2026-11-02\n"
            "**Status:** superseded\n"
            "**Summary:** x\n\n---\n")
    live = live_state_plans(["Some_Plan.json.superseded"])
    assert violations(index_statuses(text), live) == []


@pytest.mark.parametrize("status", COMPATIBLE_WITH_A_LIVE_STATE_FILE)
def test_a_compatible_status_may_keep_its_state_file(status):
    """`as-built` especially: the plan RAN and diverged, so its state file is the record of what
    happened and must not be treated as a hazard."""
    text = ("## [plan_file: Some_Plan/] 2026-11-02\n"
            "**Status:** %s\n**Summary:** x\n\n---\n" % status)
    assert violations(index_statuses(text), live_state_plans(["Some_Plan.json"])) == []


def test_an_as_built_status_with_trailing_prose_still_parses():
    """Real entries write `**Status:** as-built — Steps 1-7 EXECUTED …`. A parser that swallowed
    the prose into the status would match nothing and the guard would pass vacuously."""
    text = ("## [plan_file: Some_Plan/] 2026-11-02\n"
            "**Status:** as-built — Steps 1-7 EXECUTED and committed\n**Summary:** x\n\n---\n")
    assert index_statuses(text)["Some_Plan/"] == "as-built"


def test_a_missing_status_reads_as_None_not_as_a_violation():
    text = "## [plan_file: Some_Plan/] 2026-11-02\n**Summary:** x\n\n---\n"
    assert index_statuses(text)["Some_Plan/"] is None
    assert violations(index_statuses(text), live_state_plans(["Some_Plan.json"])) == []


def test_only_json_counts_as_resumable():
    names = ["A.json", "B.json.superseded", "C.txt"]
    assert live_state_plans(names) == {"A"}


# --------------------------------------------------------------------- the real tree


def test_no_superseded_plan_keeps_a_resumable_state_file():
    """The invariant, on the live tree. This is the one that would have caught both recurrences."""
    statuses = index_statuses(INDEX.read_text(encoding="utf-8"))
    live = live_state_plans([p.name for p in STATE_DIR.iterdir()]) if STATE_DIR.is_dir() else set()
    bad = violations(statuses, live)
    assert bad == [], (
        "these plans are marked superseded but still carry a resumable state file, so "
        "/execute-plan would offer to resume work another plan took over: %r. Rename each to "
        "<name>.json.superseded -- do NOT delete, the file is the record of which steps ran." % bad)


def test_a_superseded_plans_record_is_PRESERVED_not_deleted():
    """The remedy is a rename. A guard satisfiable by `rm` would destroy the execution record it
    exists to protect, so every superseded entry must still have a readable `.json.superseded`."""
    statuses = index_statuses(INDEX.read_text(encoding="utf-8"))
    superseded = [k.rstrip("/") for k, v in statuses.items()
                  if v in INCOMPATIBLE_WITH_A_LIVE_STATE_FILE]
    assert superseded, "non-vacuity: at least one entry must be marked superseded"
    missing = []
    for name in superseded:
        p = STATE_DIR / (name + ".json.superseded")
        if not p.is_file():
            missing.append(name)
            continue
        rec = json.loads(p.read_text(encoding="utf-8"))
        assert rec.get("steps"), "%s records no steps; the execution record was lost" % name
        assert any(s.get("status") == "done" for s in rec["steps"]), (
            "%s records no completed step, so the rename preserved nothing worth keeping" % name)
    assert missing == [], (
        "superseded plans with no preserved record: %r. The rename is what keeps the record; "
        "a delete loses which steps actually executed." % missing)


def test_the_statusless_entry_count_does_not_REGRESS():
    """A ratchet, not a requirement.

    MEASURED: 83 entries, 39 without a `**Status:**` field. Requiring all 83 would fail on 39
    pre-existing entries and satisfying it would mean inventing a status for each without reading
    the plan. Grandfathering them while failing a NEW one is the version of this check that can
    actually hold.
    """
    text = INDEX.read_text(encoding="utf-8")
    blocks = re.split(r"(?m)^## \[plan_file: ", text)[1:]
    n_without = sum(1 for b in blocks if not re.search(r"(?m)^\*\*Status:\*\*", b))
    assert n_without <= MAX_ENTRIES_WITHOUT_STATUS, (
        "%d entries now lack a **Status:** field, up from the grandfathered %d. A new plan entry "
        "should declare one of %r." % (n_without, MAX_ENTRIES_WITHOUT_STATUS,
                                       COMPATIBLE_WITH_A_LIVE_STATE_FILE))


def test_every_live_state_file_has_an_index_entry():
    """A state file with no entry at all is a plan nobody can look up. Measured: all 9 have one."""
    if not STATE_DIR.is_dir():
        pytest.skip("no execution-state directory on this checkout")
    statuses = index_statuses(INDEX.read_text(encoding="utf-8"))
    keys = {k.rstrip("/") for k in statuses}
    live = live_state_plans([p.name for p in STATE_DIR.iterdir()])
    orphans = sorted(live - keys)
    assert orphans == [], "state files with no plans-index entry: %r" % orphans


def test_this_plans_own_two_annotations_landed():
    """Pins Step 1's actual edits, so a later index rewrite cannot quietly drop them."""
    text = INDEX.read_text(encoding="utf-8")
    statuses = index_statuses(text)
    assert statuses.get(
        "Finish_The_Three_Successor_Arms_The_Expression_Heads_And_Fine_Tune_Halves_Re_Grounded_"
        "In_What_The_First_Seven_Steps_Measured_Plan/") == "superseded"
    assert statuses.get(
        "Two_Unattended_Windows_Tonight_S_10_12_Hours_And_The_10_Day_Vacation_Run_Plan/"
    ) == "as-built"
    assert "Do not re-execute:" in text
