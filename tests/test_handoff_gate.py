"""Tests for scripts/handoff_gate.py.

Covers the 5 pure-logic verdict helpers + the accept/reject aggregator +
basic markdown rendering. Subprocess-touching paths (git, pytest) are
exercised via integration tests where it's safe to do so.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from scripts.handoff_gate import (
    DEFAULT_CONTRACT_TESTS,
    GateCheck,
    GateReport,
    evaluate_locked_parameter_coverage,
    evaluate_pytest_result,
    evaluate_receiver_pull,
    evaluate_sender_push,
    evaluate_sync_verdict,
    is_handoff_accepted,
    render_markdown,
    run_all_checks,
)


# ---------------------------------------------------------------------------
# Pure-logic helpers
# ---------------------------------------------------------------------------


class TestEvaluateSenderPush:
    def test_pass_when_ahead_zero(self):
        c = evaluate_sender_push(ahead=0, behind=0)
        assert c.status == "PASS"
        assert c.name == "sender_push_status"
        assert c.extra["ahead"] == 0

    def test_pass_even_when_behind_nonzero(self):
        """The sender-push gate doesn't care about 'behind' — that's gate 2's job."""
        c = evaluate_sender_push(ahead=0, behind=5)
        assert c.status == "PASS"

    def test_fail_when_ahead_one(self):
        c = evaluate_sender_push(ahead=1, behind=0, sample=["abc123 commit msg"])
        assert c.status == "FAIL"
        assert "git push" in c.detail

    def test_fail_includes_sample_in_extra(self):
        sample = [f"hash{i} msg{i}" for i in range(10)]
        c = evaluate_sender_push(ahead=10, behind=0, sample=sample)
        assert c.status == "FAIL"
        assert len(c.extra["sample"]) == 5  # capped to 5

    def test_fail_handles_none_sample(self):
        c = evaluate_sender_push(ahead=2, behind=0, sample=None)
        assert c.status == "FAIL"
        assert c.extra["sample"] == []


class TestEvaluateReceiverPull:
    def test_pass_when_behind_zero(self):
        c = evaluate_receiver_pull(ahead=0, behind=0)
        assert c.status == "PASS"

    def test_fail_when_behind_nonzero(self):
        c = evaluate_receiver_pull(ahead=0, behind=3)
        assert c.status == "FAIL"
        assert "git pull --ff-only" in c.detail
        assert c.extra["behind"] == 3


class TestEvaluateSyncVerdict:
    def test_pass_when_all_markers_present(self):
        divergences = [
            {"file": "a.py", "marker": "x", "status": "PRESENT"},
            {"file": "b.py", "marker": "y", "status": "PRESENT"},
        ]
        c = evaluate_sync_verdict(divergences)
        assert c.status == "PASS"
        assert c.extra["missing"] == 0

    def test_fail_when_marker_missing(self):
        divergences = [
            {"file": "a.py", "marker": "x", "status": "PRESENT"},
            {"file": "b.py", "marker": "y", "status": "MISSING_LIKELY_STALE"},
        ]
        c = evaluate_sync_verdict(divergences)
        assert c.status == "FAIL"
        assert c.extra["missing"] == 1

    def test_fail_when_file_missing(self):
        divergences = [{"file": "c.py", "marker": "z", "status": "FILE_MISSING"}]
        c = evaluate_sync_verdict(divergences)
        assert c.status == "FAIL"

    def test_fail_when_read_error(self):
        divergences = [{"file": "c.py", "marker": "z", "status": "READ_ERROR: foo"}]
        c = evaluate_sync_verdict(divergences)
        assert c.status == "FAIL"

    def test_pass_when_empty_divergence_list(self):
        """Edge case: no markers tracked = vacuously in-sync."""
        c = evaluate_sync_verdict([])
        assert c.status == "PASS"
        assert c.extra["markers_total"] == 0


class TestEvaluateLockedParameterCoverage:
    def test_pass_when_all_test_files_exist(self, tmp_path: Path):
        # Build a fake repo with all required test files
        (tmp_path / "tests").mkdir()
        for t in DEFAULT_CONTRACT_TESTS:
            (tmp_path / t).parent.mkdir(parents=True, exist_ok=True)
            (tmp_path / t).write_text("# stub\n")
        c = evaluate_locked_parameter_coverage(tmp_path, DEFAULT_CONTRACT_TESTS)
        assert c.status == "PASS"

    def test_fail_when_test_file_missing(self, tmp_path: Path):
        # Empty tmp_path = all test files missing
        c = evaluate_locked_parameter_coverage(tmp_path, DEFAULT_CONTRACT_TESTS)
        assert c.status == "FAIL"
        assert len(c.extra["missing"]) == len(DEFAULT_CONTRACT_TESTS)

    def test_fail_when_one_of_three_missing(self, tmp_path: Path):
        (tmp_path / "tests").mkdir()
        # Write only the first two; third stays missing
        for t in DEFAULT_CONTRACT_TESTS[:2]:
            (tmp_path / t).write_text("# stub\n")
        c = evaluate_locked_parameter_coverage(tmp_path, DEFAULT_CONTRACT_TESTS)
        assert c.status == "FAIL"
        assert len(c.extra["missing"]) == 1


class TestEvaluatePytestResult:
    def test_pass_when_exit_zero(self):
        c = evaluate_pytest_result(0, ["", "20 passed in 0.30s"], DEFAULT_CONTRACT_TESTS)
        assert c.status == "PASS"
        assert "20 passed" in c.detail

    def test_fail_when_exit_nonzero(self):
        c = evaluate_pytest_result(
            1, ["FAILED tests/test_x.py::test_y", "1 failed, 19 passed"],
            DEFAULT_CONTRACT_TESTS,
        )
        assert c.status == "FAIL"
        assert "exit 1" in c.detail

    def test_fail_with_empty_stdout(self):
        c = evaluate_pytest_result(2, [], DEFAULT_CONTRACT_TESTS)
        assert c.status == "FAIL"


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------


class TestIsHandoffAccepted:
    def test_accepts_when_all_pass(self):
        checks = [
            GateCheck(name=str(i), description="", status="PASS")
            for i in range(5)
        ]
        assert is_handoff_accepted(checks) is True

    def test_rejects_when_one_fails(self):
        checks = [
            GateCheck(name=str(i), description="", status="PASS")
            for i in range(4)
        ]
        checks.append(GateCheck(name="bad", description="", status="FAIL"))
        assert is_handoff_accepted(checks) is False

    def test_rejects_when_error(self):
        checks = [GateCheck(name="a", description="", status="ERROR")]
        assert is_handoff_accepted(checks) is False

    def test_rejects_when_empty(self):
        """No checks = no evidence of acceptance."""
        assert is_handoff_accepted([]) is False


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------


class TestRenderMarkdown:
    def _build_report(self, accepted: bool) -> GateReport:
        status = "PASS" if accepted else "FAIL"
        checks = [
            GateCheck(name=f"check_{i}", description="", status=status, detail=f"d{i}")
            for i in range(5)
        ]
        return GateReport(
            run_date="2026-05-26",
            repo_path="/tmp/foo",
            handoff_accepted=accepted,
            checks=checks,
        )

    def test_accepted_report_says_accepted(self):
        md = render_markdown(self._build_report(accepted=True))
        assert "ACCEPTED" in md
        assert "no longer provisional" in md

    def test_provisional_report_says_provisional(self):
        md = render_markdown(self._build_report(accepted=False))
        assert "PROVISIONAL" in md
        assert "provisional" in md.lower()

    def test_table_lists_all_five_checks(self):
        md = render_markdown(self._build_report(accepted=True))
        for i in range(5):
            assert f"check_{i}" in md

    def test_suggested_actions_appear_on_failure(self):
        """Failing checks should produce next-step suggestions."""
        checks = [
            GateCheck(name="sender_push_status", description="", status="FAIL", detail=""),
            GateCheck(name="receiver_pull", description="", status="PASS"),
            GateCheck(name="sync_verdict", description="", status="PASS"),
            GateCheck(name="locked_parameter_coverage", description="", status="PASS"),
            GateCheck(name="contract_tests", description="", status="PASS"),
        ]
        report = GateReport(
            run_date="2026-05-26", repo_path="/tmp", handoff_accepted=False,
            checks=checks,
        )
        md = render_markdown(report)
        assert "git push origin main" in md


# ---------------------------------------------------------------------------
# Integration — actually run against the live repo
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_run_all_checks_on_live_repo_does_not_crash():
    """Smoke test: run_all_checks on the real repo returns a 5-check report.

    We do NOT assert PASS/FAIL — local state can vary between runs. We assert
    structure only.
    """
    repo = Path.cwd()
    if not (repo / ".git").exists():
        pytest.skip("not a git repo")
    report = run_all_checks(repo)
    assert len(report.checks) == 5
    names = [c.name for c in report.checks]
    assert names == [
        "sender_push_status",
        "receiver_pull",
        "sync_verdict",
        "locked_parameter_coverage",
        "contract_tests",
    ]
    # Every check must declare a status
    for c in report.checks:
        assert c.status in ("PASS", "FAIL", "ERROR")
