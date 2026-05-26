"""Handoff gate runner — per plans/Two_Machine_Operating_Contract.md §2.1.

Runs all 5 mandatory handoff-gate checks in sequence:

  1. Sender push status   — local repo has zero unpushed commits ahead of origin
  2. Receiver pull        — local repo is not behind origin (after `git fetch`)
  3. Sync verdict         — all KNOWN_DIVERGENCE_TARGETS markers PRESENT
  4. Locked-parameter coverage — required contract regression test files exist
  5. Contract tests pass  — those test files pass on receiver side

Returns a single accept/reject verdict. A handoff is **accepted** only when
all 5 checks PASS. Any failure marks the handoff **provisional** — per the
contract, work depending on a provisional artifact is itself provisional.

Output:
  - Console summary (per-check PASS/FAIL + bottom-line verdict)
  - reports/handoff_gate_<DATE>.md
  - reports/handoff_gate_<DATE>.json

Exit code:
  0 = handoff accepted (all 5 PASS)
  1 = handoff provisional (≥1 check FAIL)
  2 = tool-runtime error

Run from repo root:

    uv run python -m scripts.handoff_gate
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import date as _date
from pathlib import Path

# Reuse existing primitives from cross_machine_sync_check.
from scripts.cross_machine_sync_check import (
    check_commit_gap,
    check_spec_divergence,
)


# Default required contract regression test files (gate check 5).
# Same-commit rule (contract §4): when a new locked parameter is added, the
# corresponding regression test should land in the same commit AND be added
# to this tuple if a new test file is introduced.
DEFAULT_CONTRACT_TESTS: tuple[str, ...] = (
    "tests/test_drug_mechanism_phenotype_merge_contract.py",
    "tests/test_cross_machine_sync_check.py",
    "tests/test_two_machine_operating_contract.py",
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class GateCheck:
    """Result of a single gate check."""

    name: str
    description: str
    status: str  # "PASS" / "FAIL" / "ERROR"
    detail: str = ""
    extra: dict = field(default_factory=dict)


@dataclass
class GateReport:
    """Aggregated 5-check handoff-gate report."""

    run_date: str
    repo_path: str
    handoff_accepted: bool
    checks: list[GateCheck] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Pure-logic verdict helpers (tested without subprocess)
# ---------------------------------------------------------------------------


def evaluate_sender_push(ahead: int, behind: int, sample: list[str] | None = None) -> GateCheck:
    """Gate check 1: sender push status. PASS iff ahead == 0."""
    if ahead == 0:
        return GateCheck(
            name="sender_push_status",
            description="Sender pushed; no local-ahead unpushed commits",
            status="PASS",
            detail=f"ahead=0; behind={behind}",
            extra={"ahead": ahead, "behind": behind},
        )
    return GateCheck(
        name="sender_push_status",
        description="Sender pushed; no local-ahead unpushed commits",
        status="FAIL",
        detail=f"{ahead} commit(s) unpushed — sender must `git push origin main`",
        extra={"ahead": ahead, "behind": behind, "sample": (sample or [])[:5]},
    )


def evaluate_receiver_pull(ahead: int, behind: int) -> GateCheck:
    """Gate check 2: receiver pulled. PASS iff behind == 0 after fetch."""
    if behind == 0:
        return GateCheck(
            name="receiver_pull",
            description="Receiver is up to date with origin",
            status="PASS",
            detail="behind=0",
            extra={"behind": behind},
        )
    return GateCheck(
        name="receiver_pull",
        description="Receiver is up to date with origin",
        status="FAIL",
        detail=f"{behind} commit(s) behind — run `git pull --ff-only origin main`",
        extra={"behind": behind},
    )


def evaluate_sync_verdict(divergences: list[dict]) -> GateCheck:
    """Gate check 3: all KNOWN_DIVERGENCE_TARGETS markers PRESENT.

    A marker is missing if its status is anything other than PRESENT — including
    `READ_ERROR: <exception>` which carries the exception repr suffix.
    """
    missing = [
        d for d in divergences
        if d.get("status", "") != "PRESENT"
    ]
    total = len(divergences)
    if not missing:
        return GateCheck(
            name="sync_verdict",
            description="All KNOWN_DIVERGENCE_TARGETS markers PRESENT",
            status="PASS",
            detail=f"all {total} markers PRESENT",
            extra={"markers_total": total, "missing": 0},
        )
    return GateCheck(
        name="sync_verdict",
        description="All KNOWN_DIVERGENCE_TARGETS markers PRESENT",
        status="FAIL",
        detail=f"{len(missing)} of {total} markers missing or unreadable",
        extra={
            "markers_total": total,
            "missing": len(missing),
            "missing_detail": missing,
        },
    )


def evaluate_locked_parameter_coverage(
    repo: Path, contract_tests: tuple[str, ...]
) -> GateCheck:
    """Gate check 4: required contract test files exist on disk."""
    missing = [t for t in contract_tests if not (repo / t).exists()]
    if not missing:
        return GateCheck(
            name="locked_parameter_coverage",
            description="All required contract regression test files present",
            status="PASS",
            detail=f"{len(contract_tests)} contract test files present",
            extra={"test_files": list(contract_tests)},
        )
    return GateCheck(
        name="locked_parameter_coverage",
        description="All required contract regression test files present",
        status="FAIL",
        detail=f"missing: {missing}",
        extra={"missing": missing, "expected": list(contract_tests)},
    )


def evaluate_pytest_result(
    exit_code: int, stdout_tail: list[str], tests_run: tuple[str, ...]
) -> GateCheck:
    """Gate check 5: contract tests pass."""
    if exit_code == 0:
        return GateCheck(
            name="contract_tests",
            description="All contract regression tests PASS on receiver",
            status="PASS",
            detail=stdout_tail[-1] if stdout_tail else "exit 0",
            extra={"exit_code": exit_code, "tests_run": list(tests_run)},
        )
    return GateCheck(
        name="contract_tests",
        description="All contract regression tests PASS on receiver",
        status="FAIL",
        detail=f"pytest exit {exit_code}; tail: {stdout_tail[-3:] if stdout_tail else []}",
        extra={
            "exit_code": exit_code,
            "stdout_tail": stdout_tail[-10:],
            "tests_run": list(tests_run),
        },
    )


def is_handoff_accepted(checks: list[GateCheck]) -> bool:
    """A handoff is accepted iff every check returns PASS."""
    return bool(checks) and all(c.status == "PASS" for c in checks)


# ---------------------------------------------------------------------------
# Side-effectful runners (subprocess-touching)
# ---------------------------------------------------------------------------


def check_sender_push_status(repo: Path) -> GateCheck:
    ahead, behind, sample = check_commit_gap(repo)
    return evaluate_sender_push(ahead, behind, sample)


def check_receiver_pull(repo: Path) -> GateCheck:
    # check_commit_gap already runs `git fetch origin --quiet` so 'behind' is fresh.
    ahead, behind, _ = check_commit_gap(repo)
    return evaluate_receiver_pull(ahead, behind)


def check_sync_verdict(repo: Path) -> GateCheck:
    divergences = check_spec_divergence(repo)
    return evaluate_sync_verdict(divergences)


def check_locked_parameter_coverage(
    repo: Path, contract_tests: tuple[str, ...] = DEFAULT_CONTRACT_TESTS
) -> GateCheck:
    return evaluate_locked_parameter_coverage(repo, contract_tests)


def check_contract_tests_pass(
    repo: Path, contract_tests: tuple[str, ...] = DEFAULT_CONTRACT_TESTS
) -> GateCheck:
    """Run pytest on the contract test files."""
    existing = [t for t in contract_tests if (repo / t).exists()]
    if not existing:
        return GateCheck(
            name="contract_tests",
            description="All contract regression tests PASS on receiver",
            status="ERROR",
            detail="no contract test files found to run",
            extra={"expected": list(contract_tests)},
        )
    cmd = ["uv", "run", "pytest", *existing, "-q", "--tb=line"]
    try:
        p = subprocess.run(cmd, cwd=repo, capture_output=True, text=True)
    except FileNotFoundError as e:
        return GateCheck(
            name="contract_tests",
            description="All contract regression tests PASS on receiver",
            status="ERROR",
            detail=f"could not run pytest: {e!r}",
        )
    stdout_lines = [l for l in p.stdout.splitlines() if l.strip()]
    return evaluate_pytest_result(p.returncode, stdout_lines, contract_tests)


def run_all_checks(
    repo: Path, contract_tests: tuple[str, ...] = DEFAULT_CONTRACT_TESTS
) -> GateReport:
    """Run all 5 gate checks in sequence; return aggregated report."""
    checks = [
        check_sender_push_status(repo),
        check_receiver_pull(repo),
        check_sync_verdict(repo),
        check_locked_parameter_coverage(repo, contract_tests),
        check_contract_tests_pass(repo, contract_tests),
    ]
    return GateReport(
        run_date=_date.today().isoformat(),
        repo_path=str(repo),
        handoff_accepted=is_handoff_accepted(checks),
        checks=checks,
        notes=[],
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def render_markdown(report: GateReport) -> str:
    verdict_word = "ACCEPTED" if report.handoff_accepted else "PROVISIONAL"
    lines = [
        f"# Handoff Gate Report — {report.run_date}",
        "",
        f"**Verdict:** `{verdict_word}`",
        f"**Repo:** `{report.repo_path}`",
        "",
        "Per `plans/Two_Machine_Operating_Contract.md` §2.1.",
        "",
        "| # | Check | Status | Detail |",
        "|---|---|---|---|",
    ]
    for i, c in enumerate(report.checks, start=1):
        detail = c.detail.replace("|", "\\|")
        lines.append(f"| {i} | {c.name} | `{c.status}` | {detail} |")
    if not report.handoff_accepted:
        lines.extend([
            "",
            "## Bottom line",
            "",
            "Handoff is **provisional**. Per contract §2.1: work depending on this",
            "artifact is itself provisional until the failing check(s) resolve.",
        ])
        # Suggest next action per failing check
        lines.append("")
        lines.append("### Suggested next actions")
        for c in report.checks:
            if c.status == "PASS":
                continue
            if c.name == "sender_push_status":
                lines.append("- Sender: `git push origin main`")
            elif c.name == "receiver_pull":
                lines.append("- Receiver: `git pull --ff-only origin main`")
            elif c.name == "sync_verdict":
                lines.append("- Verify the missing markers; sender pushes the lock; receiver re-runs gate")
            elif c.name == "locked_parameter_coverage":
                lines.append("- Add the missing contract regression test files in the SAME commit that locked the parameter")
            elif c.name == "contract_tests":
                lines.append("- Investigate failing test(s); contract tests must pass for handoff to be accepted")
    else:
        lines.extend([
            "",
            "## Bottom line",
            "",
            "Handoff is **accepted**. Work depending on this artifact is no longer provisional.",
        ])
    return "\n".join(lines) + "\n"


def write_report(
    report: GateReport, output_dir: Path
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    md_path = output_dir / f"handoff_gate_{report.run_date}.md"
    json_path = output_dir / f"handoff_gate_{report.run_date}.json"
    md_path.write_text(render_markdown(report), encoding="utf-8")
    json_path.write_text(json.dumps(asdict(report), indent=2), encoding="utf-8")
    return md_path, json_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run all 5 handoff-gate checks per plans/Two_Machine_Operating_Contract.md §2.1"
    )
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Where to write the report. Defaults to <repo>/reports/.",
    )
    parser.add_argument(
        "--json-only",
        action="store_true",
        help="Skip writing the markdown report.",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Print to console only; no report files written.",
    )
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    if not (repo / ".git").exists():
        print(f"[handoff-gate] error: {repo} is not a git repo", file=sys.stderr)
        return 2

    try:
        report = run_all_checks(repo)
    except Exception as e:
        print(f"[handoff-gate] tool-runtime error: {e!r}", file=sys.stderr)
        return 2

    # Console summary — ASCII-only so Windows cp1252 console doesn't choke.
    verdict = "ACCEPTED" if report.handoff_accepted else "PROVISIONAL"
    print(f"=== Handoff Gate ({report.run_date}) -- {verdict} ===")
    for i, c in enumerate(report.checks, start=1):
        marker = {"PASS": "[ OK ]", "FAIL": "[FAIL]", "ERROR": "[ERR ]"}.get(c.status, "[ ?? ]")
        print(f"  {marker} [{i}] {c.name:30s} {c.status:6s}  {c.detail}")
    print()

    if not args.no_write:
        out_dir = args.output_dir or (repo / "reports")
        md_path, json_path = write_report(report, out_dir)
        if args.json_only:
            md_path.unlink(missing_ok=True)
        else:
            print(f"[handoff-gate] wrote {md_path}")
        print(f"[handoff-gate] wrote {json_path}")

    return 0 if report.handoff_accepted else 1


if __name__ == "__main__":
    sys.exit(main())
