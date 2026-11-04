"""Tests for the fine-tune run's relaunch watchdog.

THE ORDER OF THE CHECKS IS THE SAFETY PROPERTY, so it is asserted directly rather than inferred
from happy paths. Two mistakes would each cost the ten-day window in a different direction: racing
a healthy run would stack a second torch job (this host's documented access-violation hazard), and
treating any stop as final would leave a resumable run dead for days. So `decide` is tested as a
total function over the cross-product of its inputs, not just its intended branches.

THE MALFORMED-STATE DIRECTION IS THE ONE THAT MATTERS. A damaged watchdog state file must read as
EMPTY, never as "budget spent" -- the latter silently disables the watchdog, which is exactly the
failure that costs the window. That is the same malformed-versus-absent distinction that cost 3.8
hours when a NUL-filled sidecar read as a run that never happened.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import scripts.watchdog_finetune_arm as W
from scripts.watchdog_finetune_arm import (
    ALL_DECISIONS,
    DONE,
    EXHAUSTED,
    NOT_OWNER,
    NO_LAUNCHER,
    UNREACHABLE,
    RELAUNCH,
    RUNNING,
    decide,
    read_state,
    write_state,
)

HOLDER = {"pid": 4242, "label": "finetune_arm"}

#: The owning host, for the calls that are allowed to reach a RELAUNCH decision.
OWNED = {"owner": "THE-OWNER", "this_host": "THE-OWNER"}


def _artifact(tmp_path, present: bool) -> Path:
    """The run artifact, created OR REMOVED to match `present`.

    The remove branch is load-bearing: two calls in one test share this path, so a helper that only
    ever created the file left it behind and a later `present=False` call silently described the
    opposite state. That made one test assert RUNNING against an artifact that was still on disk.
    """
    p = tmp_path / "run.json"
    if present:
        p.write_text('{"schema": "finetune-arm-run-v1"}', encoding="utf-8")
    elif p.exists():
        p.unlink()
    return p


def _launcher(tmp_path, present: bool = True) -> Path:
    p = tmp_path / "run_finetune_arm_detached.bat"
    if present:
        p.write_text("@echo off\n", encoding="ascii")
    return p


# ------------------------------------------------------------------ the four decisions


def test_an_existing_artifact_means_DONE(tmp_path):
    d = decide(artifact=_artifact(tmp_path, True), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path))
    assert d["decision"] == DONE


def test_a_LIVE_lock_holder_means_RUNNING(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=HOLDER, relaunches=0,
               launcher=_launcher(tmp_path))
    assert d["decision"] == RUNNING
    assert "4242" in d["reason"]


def test_no_artifact_and_no_holder_means_RELAUNCH(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path), **OWNED)
    assert d["decision"] == RELAUNCH
    assert "resumes from the last per-epoch checkpoint" in d["reason"]


def test_the_relaunch_budget_is_BOUNDED(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=None, relaunches=12,
               launcher=_launcher(tmp_path), max_relaunches=12, **OWNED)
    assert d["decision"] == EXHAUSTED
    assert "needs a human" in d["reason"]


def test_an_absent_launcher_REFUSES_rather_than_improvising(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path, present=False))
    assert d["decision"] == NO_LAUNCHER


# ------------------------------------------------------------------ the ORDER of the checks


def test_a_FINISHED_run_is_never_relaunched_even_with_budget_left(tmp_path):
    """The first check must be the artifact: a finished run has no lock holder, so any later check
    would conclude RELAUNCH and overwrite a good result."""
    d = decide(artifact=_artifact(tmp_path, True), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path))
    assert d["decision"] == DONE


def test_a_FINISHED_run_outranks_an_exhausted_budget(tmp_path):
    d = decide(artifact=_artifact(tmp_path, True), lock_holder=None, relaunches=999,
               launcher=_launcher(tmp_path), max_relaunches=12)
    assert d["decision"] == DONE


def test_a_LIVE_run_outranks_the_budget_so_it_is_never_RACED(tmp_path):
    """Stacking a second torch job beside a live one is the hazard the run lock exists for; the
    watchdog must not be the thing that causes it."""
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=HOLDER, relaunches=999,
               launcher=_launcher(tmp_path), max_relaunches=12)
    assert d["decision"] == RUNNING


def test_a_LIVE_run_outranks_an_absent_launcher(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=HOLDER, relaunches=0,
               launcher=_launcher(tmp_path, present=False))
    assert d["decision"] == RUNNING


@pytest.mark.parametrize("has_art", [True, False])
@pytest.mark.parametrize("has_holder", [True, False])
@pytest.mark.parametrize("relaunches", [0, 5, 12, 99])
@pytest.mark.parametrize("has_launcher", [True, False])
@pytest.mark.parametrize("owns", [True, False])
def test_decide_is_TOTAL_and_only_RELAUNCH_ever_acts(tmp_path, has_art, has_holder,
                                                     relaunches, has_launcher, owns):
    """Across the whole cross-product: always one of the declared decisions, and RELAUNCH only
    when there is genuinely nothing running and nothing finished."""
    d = decide(artifact=_artifact(tmp_path, has_art),
               lock_holder=(HOLDER if has_holder else None),
               relaunches=relaunches, launcher=_launcher(tmp_path, has_launcher),
               max_relaunches=12,
               owner=("THE-OWNER" if owns else "SOMEONE-ELSE"), this_host="THE-OWNER")
    assert d["decision"] in ALL_DECISIONS
    if d["decision"] == RELAUNCH:
        assert not has_art and not has_holder and has_launcher and relaunches < 12 and owns


# ------------------------------------------------------------------ the state file


def test_an_absent_state_reads_as_zero_relaunches(tmp_path):
    st = read_state(tmp_path)
    assert st["relaunches"] == 0 and st["history"] == []


def test_a_MALFORMED_state_reads_as_EMPTY_not_as_EXHAUSTED(tmp_path):
    """The load-bearing direction. Reading damage as 'budget spent' would silently disable the
    watchdog, which is the failure that costs the window."""
    W.state_path(tmp_path).parent.mkdir(parents=True, exist_ok=True)
    W.state_path(tmp_path).write_text('{"relaunches": 7, "history"', encoding="utf-8")
    assert read_state(tmp_path)["relaunches"] == 0


def test_a_NUL_filled_state_reads_as_EMPTY(tmp_path):
    """The signature a power cut leaves on this host, twice measured."""
    W.state_path(tmp_path).parent.mkdir(parents=True, exist_ok=True)
    W.state_path(tmp_path).write_bytes(b"\x00" * 64)
    assert read_state(tmp_path)["relaunches"] == 0


def test_a_state_that_is_not_an_object_reads_as_EMPTY(tmp_path):
    W.state_path(tmp_path).parent.mkdir(parents=True, exist_ok=True)
    W.state_path(tmp_path).write_text("[1, 2, 3]", encoding="utf-8")
    assert read_state(tmp_path)["relaunches"] == 0


def test_write_state_round_trips_and_leaves_no_tmp(tmp_path):
    write_state({"relaunches": 3, "history": []}, tmp_path)
    assert read_state(tmp_path)["relaunches"] == 3
    assert not list(tmp_path.glob("*.tmp"))


def test_write_state_uses_fsync_and_replace():
    """Asserted over the function's ast CALL set: its docstring must name the durable-write pair,
    so a text search for the absence of `write_text` is unsatisfiable."""
    import ast

    src = (Path(W.__file__)).read_text(encoding="utf-8")
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "write_state")
    called = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Call):
            f = node.func
            called.add(f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", ""))
    assert {"fsync", "replace"} <= called
    assert "write_text" not in called


# ------------------------------------------------------------------ main, and what it records


def test_a_dry_run_starts_nothing_and_writes_no_state(tmp_path, monkeypatch):
    started = []
    monkeypatch.setattr(W, "launch", lambda *a, **k: started.append(1))
    monkeypatch.setattr(W, "ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(W, "LAUNCHER", _launcher(tmp_path))
    assert W.main(["--state-dir", str(tmp_path)]) == 0
    assert started == []
    assert not W.state_path(tmp_path).is_file(), "a dry run wrote state"


def test_apply_relaunches_once_and_increments_the_counter(tmp_path, monkeypatch):
    monkeypatch.setattr(W, "launch", lambda *a, **k: 31337)
    monkeypatch.setattr(W, "ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(W, "LAUNCHER", _launcher(tmp_path))
    W.main(["--claim", "--state-dir", str(tmp_path)])
    W.main(["--apply", "--state-dir", str(tmp_path)])
    st = read_state(tmp_path)
    assert st["relaunches"] == 1
    assert st["history"][-1]["launched_pid"] == 31337
    assert st["last_decision"] == RELAUNCH


def test_apply_stops_relaunching_once_the_budget_is_spent(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(W, "launch", lambda *a, **k: calls.append(1) or 1)
    monkeypatch.setattr(W, "ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(W, "LAUNCHER", _launcher(tmp_path))
    W.main(["--claim", "--state-dir", str(tmp_path)])
    for _ in range(5):
        W.main(["--apply", "--state-dir", str(tmp_path), "--max-relaunches", "2"])
    assert len(calls) == 2, "the bound did not hold: %d launches" % len(calls)
    assert read_state(tmp_path)["last_decision"] == EXHAUSTED


def test_every_decision_is_recorded_not_only_the_actions(tmp_path, monkeypatch):
    """A history of `running` entries is how a reader later confirms the watchdog was alive the
    whole time, rather than having died on day two itself."""
    monkeypatch.setattr(W, "ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(W, "LAUNCHER", _launcher(tmp_path))
    monkeypatch.setattr("scripts.glm_run_lock.holder", lambda p: HOLDER)
    W.main(["--apply", "--state-dir", str(tmp_path)])
    st = read_state(tmp_path)
    assert st["relaunches"] == 0
    assert st["history"][-1]["decision"] == RUNNING


def test_the_history_is_capped_so_it_cannot_grow_without_bound(tmp_path, monkeypatch):
    monkeypatch.setattr(W, "ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(W, "LAUNCHER", _launcher(tmp_path))
    monkeypatch.setattr("scripts.glm_run_lock.holder", lambda p: HOLDER)
    write_state({"relaunches": 0, "history": [{"at": "x"}] * 500}, tmp_path)
    W.main(["--apply", "--state-dir", str(tmp_path)])
    assert len(read_state(tmp_path)["history"]) <= 200


def test_apply_does_nothing_when_the_artifact_is_already_there(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(W, "launch", lambda *a, **k: calls.append(1) or 1)
    monkeypatch.setattr(W, "ARTIFACT", _artifact(tmp_path, True))
    monkeypatch.setattr(W, "LAUNCHER", _launcher(tmp_path))
    W.main(["--apply", "--state-dir", str(tmp_path)])
    assert calls == []
    assert read_state(tmp_path)["last_decision"] == DONE


# ------------------------------------------------------------------ the real repo wiring


def test_the_watchdog_points_at_the_real_launcher_and_artifact_names():
    assert W.LAUNCHER.name == "run_finetune_arm_detached.bat"
    assert W.ARTIFACT.name == "finetune_arm_run_2026-11-13.json"
    assert W.state_path().name == "finetune_arm.watchdog.json"


def test_the_lock_path_matches_the_launchers_own():
    """The watchdog and the launcher must agree on which lock file means "running"; disagreeing
    would make the watchdog blind to a live run."""
    bat = (W.REPO / "run_finetune_arm_detached.bat")
    if not bat.is_file():
        pytest.skip("launcher absent on this checkout")
    text = bat.read_text(encoding="ascii", errors="replace")
    assert "finetune_arm.lock" in text

# ------------------------------------------------------------------ the ownership guard


def test_an_UNCLAIMED_run_is_never_relaunched(tmp_path):
    """The footgun this guard exists for. `run_finetune_arm_detached.bat` is COMMITTED, so it is
    present in every checkout, and its cohort cache resolves on more than one machine -- so without
    ownership every other check passes and the watchdog would start a four-thread torch job on
    whichever host happened to run it, including the one that hard-powers-off under this load."""
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path), owner=None, this_host="SOME-HOST")
    assert d["decision"] == NOT_OWNER
    assert "unclaimed" in d["reason"]


def test_a_run_owned_by_ANOTHER_host_is_never_relaunched_here(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path), owner="SERVER2", this_host="THE-LAPTOP")
    assert d["decision"] == NOT_OWNER
    assert "SERVER2" in d["reason"] and "THE-LAPTOP" in d["reason"]


def test_ownership_does_NOT_suppress_the_read_only_facts(tmp_path):
    """DONE and RUNNING are facts about the run, true on any host, and a reader on the wrong
    machine still needs them -- only the ACTING branches are gated."""
    assert decide(artifact=_artifact(tmp_path, True), lock_holder=None, relaunches=0,
                  launcher=_launcher(tmp_path), owner="OTHER",
                  this_host="HERE")["decision"] == DONE
    assert decide(artifact=_artifact(tmp_path, False), lock_holder=HOLDER, relaunches=0,
                  launcher=_launcher(tmp_path), owner="OTHER",
                  this_host="HERE")["decision"] == RUNNING


def test_claim_records_THIS_host_and_then_relaunching_is_permitted(tmp_path):
    import platform

    assert W.read_owner(tmp_path) is None
    assert W.main(["--claim", "--state-dir", str(tmp_path)]) == 0
    assert W.read_owner(tmp_path) == platform.node()


def test_the_LIVE_repo_is_NOT_claimed_so_this_laptop_cannot_relaunch():
    """Guard on the real state dir: the run lives on server2, and this checkout must never be able
    to start it. If someone deliberately claims it here, this fails and says so."""
    assert W.read_owner() is None, (
        "this host has claimed the fine-tune run; it hard-powered-off three times in eighteen "
        "hours under exactly this load")


# ------------------------------------------------------------------ the precondition check


def _launcher_with_env(tmp_path, pa: str, mpra: str) -> Path:
    """A launcher in the shape the generator emits: quoted `set` lines, outputs plus paths."""
    p = tmp_path / "run_finetune_arm_detached.bat"
    rows = [
        "@echo off",
        'set "LOG=%~dp0data/unattended/finetune_arm.log"',
        'set "SENTINEL=%~dp0data/unattended/finetune_arm.sentinel"',
        'set "LOCK=%~dp0data/unattended/finetune_arm.lock"',
        'set "MPRA_CACHE=%s"' % mpra,
        'set "PA_REPO=%s"' % pa,
    ]
    p.write_text("\n".join(rows) + "\n", encoding="ascii")
    return p


def test_declared_env_paths_reads_the_PATH_variables_from_the_launcher(tmp_path):
    lp = _launcher_with_env(tmp_path, str(tmp_path), str(tmp_path))
    assert W.declared_env_paths(lp) == {"PA_REPO": str(tmp_path), "MPRA_CACHE": str(tmp_path)}


def test_declared_env_paths_EXCLUDES_the_launchers_own_outputs(tmp_path):
    """LOG, SENTINEL and LOCK are created BY the launcher. Requiring them to pre-exist would make
    the watchdog refuse every first launch."""
    got = W.declared_env_paths(_launcher_with_env(tmp_path, str(tmp_path), str(tmp_path)))
    for out in ("LOG", "SENTINEL", "LOCK"):
        assert out not in got


def test_unreachable_env_names_only_the_paths_that_do_not_resolve(tmp_path):
    lp = _launcher_with_env(tmp_path, str(tmp_path), str(tmp_path / "gone"))
    bad = W.unreachable_env(lp)
    assert set(bad) == {"MPRA_CACHE"}


def test_an_unreachable_declared_path_REFUSES_and_names_it(tmp_path):
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path), unreachable={"MPRA_CACHE": "D:/gone"}, **OWNED)
    assert d["decision"] == UNREACHABLE
    assert "MPRA_CACHE" in d["reason"]


def test_an_unreachable_volume_does_NOT_spend_the_relaunch_budget(tmp_path, monkeypatch):
    """The point of the check. A detached volume cannot be fixed by relaunching, so burning 12
    attempts on it would leave nothing for the crash the budget exists for."""
    calls = []
    monkeypatch.setattr(W, "launch", lambda *a, **k: calls.append(1) or 1)
    monkeypatch.setattr(W, "ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(W, "LAUNCHER", _launcher_with_env(tmp_path, str(tmp_path / "gone"),
                                                          str(tmp_path)))
    W.main(["--claim", "--state-dir", str(tmp_path)])
    for _ in range(4):
        W.main(["--apply", "--state-dir", str(tmp_path)])
    assert calls == []
    st = read_state(tmp_path)
    assert st["relaunches"] == 0, "an unreachable volume spent the budget"
    assert st["last_decision"] == UNREACHABLE


def test_a_LIVE_run_outranks_an_unreachable_path(tmp_path):
    """Never disturb a running job over a path check -- the job already has what it opened."""
    d = decide(artifact=_artifact(tmp_path, False), lock_holder=HOLDER, relaunches=0,
               launcher=_launcher(tmp_path), unreachable={"MPRA_CACHE": "D:/gone"}, **OWNED)
    assert d["decision"] == RUNNING


def test_a_FINISHED_run_outranks_an_unreachable_path(tmp_path):
    d = decide(artifact=_artifact(tmp_path, True), lock_holder=None, relaunches=0,
               launcher=_launcher(tmp_path), unreachable={"MPRA_CACHE": "D:/gone"}, **OWNED)
    assert d["decision"] == DONE


def test_the_REAL_committed_launcher_declares_both_path_variables():
    """Real-data check: if the generator stops emitting these, the precondition check would go
    silently vacuous -- an empty unreachable set reads as 'everything resolves'."""
    bat = W.REPO / "run_finetune_arm_detached.bat"
    if not bat.is_file():
        pytest.skip("launcher absent on this checkout")
    got = W.declared_env_paths(bat)
    assert set(got) == {"PA_REPO", "MPRA_CACHE"}, got
