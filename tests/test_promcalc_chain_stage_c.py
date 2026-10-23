"""Guards for the unattended stage-B -> validate -> stage-C chain (`scripts/promcalc_chain_stage_c.py`).

The chain's whole value is a GATE ORDER, so that is what these tests pin: the validator runs BEFORE
stage C, and a non-zero validator means stage C is never invoked at all. Both halves matter. If the
order inverted, or if a blocked validator still let the scorer fire, the morning artifact would be
`INDETERMINATE_RECONCILE_FAILED` on a defective table -- which is indistinguishable from genuine
model drift, so the reader would wake up to a plausible wrong diagnosis instead of a defect report.

A chain only ever observed on the happy path proves nothing (the lesson the smoke gate taught), so
every refusal branch -- timed out, workers still alive, worker count UNKNOWN, validator blocked,
stage C non-zero, lock already held -- is exercised, and each asserts WHAT WAS NOT RUN rather than
only what was.

`run()` and `workers_alive()` are the only impure seams; both are injected. Nothing here touches the
real `data/promcalc/` table (a precious 12 h run), spawns PowerShell, or sleeps meaningfully.
"""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CHAIN = REPO / "scripts" / "promcalc_chain_stage_c.py"
LAUNCHER = REPO / "run_promcalc_chain_detached.bat"


def _load():
    spec = importlib.util.spec_from_file_location("_promcalc_chain", CHAIN)
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture
def chain(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """The module with every file path and every wall-clock constant redirected somewhere harmless."""
    m = _load()
    monkeypatch.setattr(m, "LOG", tmp_path / "chain.log")
    monkeypatch.setattr(m, "LOCK", tmp_path / "chain.lock")
    monkeypatch.setattr(m, "VERDICT", tmp_path / "verdict.txt")
    monkeypatch.setattr(m, "SENTINEL", tmp_path / "sentinel.txt")
    monkeypatch.setattr(m, "POLL_S", 0.01)
    monkeypatch.setattr(m, "SETTLE_S", 0)
    monkeypatch.setattr(m, "MAX_WAIT_S", 10)
    return m


@pytest.fixture
def calls(chain, monkeypatch: pytest.MonkeyPatch):
    """Record every subprocess the chain would have launched, in order, and run none of them."""
    seen: list[tuple[str, list[str]]] = []
    rc_by_tag = {"validate": 0, "stage-C": 0}

    def fake_run(cmd: list[str], tag: str) -> int:
        seen.append((tag, list(cmd)))
        return rc_by_tag[tag]

    monkeypatch.setattr(chain, "run", fake_run)
    monkeypatch.setattr(chain, "stage_b_done", lambda: (True, "test: complete"))
    monkeypatch.setattr(chain, "workers_alive", lambda: 0)
    return seen, rc_by_tag


def _tags(calls) -> list[str]:
    return [tag for tag, _ in calls[0]]


# --- THE gate order: validate, THEN score -------------------------------------------------------

def test_a_nonzero_validator_means_stage_C_is_NEVER_invoked(chain, calls):
    """The non-vacuity test of the whole script. A defective table must not be scored at all --
    not scored-and-flagged, not scored-and-ignored. Never invoked."""
    calls[1]["validate"] = 2
    assert chain.main([]) == 1
    assert _tags(calls) == ["validate"], "stage C fired on a table the validator rejected"
    assert "BLOCKED_VALIDATOR" in chain.VERDICT.read_text(encoding="utf-8")


def test_an_indeterminate_validator_also_blocks_rather_than_guessing(chain, calls):
    """Exit 3 is 'unknown', and an unknown table is not a scoreable one -- the fail direction."""
    calls[1]["validate"] = 3
    assert chain.main([]) == 1
    assert _tags(calls) == ["validate"]


def test_the_validator_runs_BEFORE_stage_C_and_not_the_reverse(chain, calls):
    assert chain.main([]) == 0
    assert _tags(calls) == ["validate", "stage-C"]


def test_the_validator_is_invoked_WITHOUT_expect_partial(chain, calls):
    """`--expect-partial` downgrades the coverage check. Passing it here would let the chain score a
    half-finished table, in which missing tiles enter stage C as no-promoter observations."""
    chain.main([])
    validate_cmd = next(cmd for tag, cmd in calls[0] if tag == "validate")
    assert "--expect-partial" not in validate_cmd
    assert any("promcalc_validate_output.py" in part for part in validate_cmd)


def test_stage_C_runs_under_uv_because_it_needs_the_project_env(chain, calls):
    """The chain itself is stdlib-only and runs under the plain interpreter; stage C is not."""
    chain.main([])
    stage_c = next(cmd for tag, cmd in calls[0] if tag == "stage-C")
    assert stage_c[:2] == ["uv", "run"]
    assert any("promcalc_decomposition_arm.py" in part for part in stage_c)


def test_a_clean_run_scores_and_exits_zero(chain, calls):
    assert chain.main([]) == 0
    assert chain.VERDICT.read_text(encoding="utf-8").splitlines()[1] == "SCORED"


def test_a_nonzero_stage_C_is_its_own_verdict_not_a_validator_block(chain, calls):
    """A stage-C refusal is a legitimate outcome, and conflating it with a blocked validator would
    send the morning reader to the wrong half of the pipeline."""
    calls[1]["stage-C"] = 1
    assert chain.main([]) == 1
    assert _tags(calls) == ["validate", "stage-C"]
    body = chain.VERDICT.read_text(encoding="utf-8")
    assert "STAGE_C_NONZERO" in body and "BLOCKED_VALIDATOR" not in body


# --- refusal branches, each asserting what was NOT run ------------------------------------------

def test_workers_that_never_clear_block_validation_entirely(chain, calls, monkeypatch):
    """Decompressing the table during a live sweep is the recorded trap that turned a ~1 h job into
    a 5.0 h projection -- observing the sweep slowed the sweep."""
    monkeypatch.setattr(chain, "workers_alive", lambda: 2)
    monkeypatch.setattr(chain, "MAX_WAIT_S", 0.05)
    assert chain.main([]) == 1
    assert calls[0] == [], "validation started while the sweep was still writing"
    assert "TIMED_OUT" in chain.VERDICT.read_text(encoding="utf-8")


def test_an_UNKNOWN_worker_count_is_treated_as_alive_and_never_proceeds(chain, calls, monkeypatch):
    """`workers_alive()` returns -1 when it cannot tell. Over-counting only makes the chain wait,
    which is the safe direction; treating -1 as 'none alive' would validate mid-run."""
    monkeypatch.setattr(chain, "workers_alive", lambda: -1)
    monkeypatch.setattr(chain, "MAX_WAIT_S", 0.05)
    assert chain.main([]) == 1
    assert calls[0] == []


def test_stage_b_never_completing_times_out_and_runs_nothing(chain, calls, monkeypatch):
    """A stage-B crash must leave the chain waiting -- the watchdog owns relaunching, and a partial
    table must never be scored just because the night ran out."""
    monkeypatch.setattr(chain, "stage_b_done", lambda: (False, "coverage 7000/44106"))
    monkeypatch.setattr(chain, "MAX_WAIT_S", -1)
    assert chain.main([]) == 1
    assert calls[0] == []
    assert "stage B not complete" in chain.VERDICT.read_text(encoding="utf-8")


# --- single-instance lock -----------------------------------------------------------------------

def test_an_existing_lock_refuses_to_double_run(chain, calls):
    """Two chains would validate and score the same table concurrently."""
    chain.LOCK.write_text("1234", encoding="utf-8")
    assert chain.main([]) == 3
    assert calls[0] == []
    assert not chain.VERDICT.exists(), "a refused double-run must not overwrite a real verdict"
    assert chain.LOCK.read_text(encoding="utf-8") == "1234", (
        "the refused run must leave the HOLDER's lock intact -- releasing it in a finally would free "
        "the lock out from under a chain that is still running, which is worse than not starting"
    )


def test_the_lock_is_released_even_when_the_validator_blocks(chain, calls):
    """Leaving it behind would make every later night refuse to start."""
    calls[1]["validate"] = 2
    chain.main([])
    assert not chain.LOCK.exists()


def test_the_lock_is_released_on_the_happy_path_too(chain, calls):
    chain.main([])
    assert not chain.LOCK.exists()


def test_dry_run_decides_nothing_runs_nothing_and_takes_no_lock(chain, calls):
    assert chain.main(["--dry-run"]) == 0
    assert calls[0] == []
    assert not chain.LOCK.exists()
    assert not chain.VERDICT.exists()


# --- completion detection -----------------------------------------------------------------------

def test_sentinel_exit_reads_the_LAST_exit_line(chain):
    """The sentinel is appended to across watchdog relaunches, so the first line is stale."""
    chain.SENTINEL.write_text("EXIT=1\nrelaunching\nEXIT=0\n", encoding="utf-8")
    assert chain.sentinel_exit() == "EXIT=0"


def test_sentinel_exit_is_None_when_the_file_is_absent(chain):
    assert chain.sentinel_exit() is None


def test_sentinel_exit_is_None_when_stage_b_was_killed_before_writing_one(chain):
    chain.SENTINEL.write_text("started\nworking\n", encoding="utf-8")
    assert chain.sentinel_exit() is None


def test_sentinel_EXIT_0_alone_is_enough(chain, monkeypatch):
    """Stage B's own clean exit is authoritative; no coverage read is needed to trust it."""
    chain.SENTINEL.write_text("EXIT=0\n", encoding="utf-8")

    def boom() -> int:
        raise AssertionError("coverage must not be read when the sentinel already says EXIT=0")

    monkeypatch.setattr(chain, "completed_tiles", boom)
    done, why = chain.stage_b_done()
    assert done and "EXIT=0" in why


def test_a_crashed_stage_b_with_full_coverage_is_still_done(chain, monkeypatch):
    """The second route exists precisely because a killed worker never writes the exit line."""
    chain.SENTINEL.write_text("EXIT=1\n", encoding="utf-8")
    monkeypatch.setattr(chain, "completed_tiles", lambda: chain.EXPECT_TILES)
    done, why = chain.stage_b_done()
    assert done and str(chain.EXPECT_TILES) in why


def test_short_coverage_is_not_done(chain, monkeypatch):
    monkeypatch.setattr(chain, "completed_tiles", lambda: chain.EXPECT_TILES - 1)
    done, _ = chain.stage_b_done()
    assert not done, "one tile short is still a partial table"


def test_an_unreadable_coverage_is_not_done_rather_than_a_crash(chain, monkeypatch):
    """A transient read failure must keep the chain waiting, not end the night."""
    monkeypatch.setattr(chain, "completed_tiles", lambda: None)
    done, why = chain.stage_b_done()
    assert not done
    assert "None" in why


def test_completed_tiles_delegates_to_stage_Bs_own_resume_logic(chain, monkeypatch):
    """Re-implementing the resume rule here would let the two drift apart silently."""
    stub = types.ModuleType("promcalc_decompose")
    stub.completed_tiles = lambda d: ["t1", "t2", "t3"]  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "promcalc_decompose", stub)
    assert chain.completed_tiles() == 3


def test_completed_tiles_returns_None_rather_than_raising(chain, monkeypatch):
    stub = types.ModuleType("promcalc_decompose")

    def raiser(d):
        raise RuntimeError("shard half-written")

    stub.completed_tiles = raiser  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "promcalc_decompose", stub)
    assert chain.completed_tiles() is None


def test_the_worker_query_stays_deliberately_OVER_broad(chain, monkeypatch):
    """`uv venv`'s python.exe is a 46 KB trampoline, so multiprocessing spawn launches the BASE
    interpreter -- a venv-path grep would miss every worker and the chain would decompress the table
    mid-sweep. Tightening this query is the failure mode, so the breadth is pinned rather than the
    count. No PowerShell is spawned: the subprocess seam is refused outright."""
    captured: list[list[str]] = []

    def no_subprocess(cmd, **kw):
        captured.append(list(cmd))
        raise AssertionError("no PowerShell in tests")

    monkeypatch.setattr(chain.subprocess, "run", no_subprocess)
    assert chain.workers_alive() == -1, "an unanswerable query must report UNKNOWN, not zero"
    query = " ".join(captured[0])
    assert "promcalc_decompose" in query
    assert "multiprocessing" in query, (
        "spawned workers carry no script name on their command line; dropping this clause is the "
        "regression that would let validation start during a live sweep"
    )


# --- verdict reporting --------------------------------------------------------------------------

def test_finish_exit_code_is_zero_only_for_SCORED(chain):
    assert chain.finish("SCORED", "ok") == 0
    for verdict in ("BLOCKED_VALIDATOR", "TIMED_OUT", "STAGE_C_NONZERO"):
        assert chain.finish(verdict, "why") == 1, f"{verdict} must not exit 0"


def test_the_verdict_file_records_the_verdict_and_its_detail(chain):
    chain.finish("BLOCKED_VALIDATOR", "validator exit=2")
    lines = chain.VERDICT.read_text(encoding="utf-8").splitlines()
    assert lines[1] == "BLOCKED_VALIDATOR"
    assert "exit=2" in lines[2]


# --- the launcher's cmd quirk -------------------------------------------------------------------

def test_the_launcher_puts_the_redirect_BEFORE_the_echo():
    """Written `echo CHAIN_EXIT=%RC%>>file`, cmd reads the digit before `>>` as a file handle and
    the line vanishes -- the sentinel would be silently empty, which reads as 'never ran'."""
    body = LAUNCHER.read_text(encoding="utf-8", errors="replace")
    sentinel_lines = [ln.strip() for ln in body.splitlines()
                      if "CHAIN_EXIT" in ln and not ln.strip().upper().startswith("REM")]
    assert sentinel_lines, "the launcher no longer writes a chain sentinel"
    for ln in sentinel_lines:
        assert ln.startswith(">>"), f"redirect must precede the echo: {ln!r}"
        assert ln.index(">>") < ln.lower().index("echo")
