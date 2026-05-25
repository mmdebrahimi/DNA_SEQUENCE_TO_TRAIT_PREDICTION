"""Cross-machine sync diagnostic.

Symptom this addresses: 2026-05-24 — Codex on the Precision 7780 shipped
v0 + relocked the v0 spec but didn't push to origin. The GTX 860M laptop
+ Precision 7780 silently diverged. The user shared two markdown files via
Downloads as the only transfer signal; meanwhile the actual model retrain
+ pipeline.py changes + release packet + reports/ artifacts stayed on
Precision 7780 only.

This tool reports drift on five axes:
  1. Commit gap: origin vs local commit-hash diff
  2. Working-tree dirtiness: modified/untracked files vs HEAD
  3. Untracked artifact landings: Downloads/ files referencing recent dates
     that haven't made it into repo
  4. Spec-divergence spot-check: known-cross-machine docs (e.g.,
     wiki/decoder_v0_ux_and_success_criterion.md) diffed against a manifest
     of expected fields/headlines
  5. Test count snapshot: pytest --collect-only for regression-vs-prior

Output:
  - Console summary
  - reports/cross_machine_sync_<DATE>.md (narrative)
  - reports/cross_machine_sync_<DATE>.json (machine-readable)

Exit code:
  0 = in sync
  1 = drift detected (one or more axes)
  2 = tool-runtime error
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import date as _date
from pathlib import Path


# Known cross-machine artifacts that historically diverge.
KNOWN_DIVERGENCE_TARGETS: tuple[tuple[str, str], ...] = (
    ("wiki/decoder_v0_ux_and_success_criterion.md", "RELOCKED"),  # Codex relocked 2026-05-23
    ("scripts/pipeline.py", "cv_strategy"),                       # Codex added cv_strategy field
    ("scripts/pipeline.py", "leave_one_accession_out"),           # Codex's leakage-safe CV
)


@dataclass
class DriftReport:
    """Aggregated drift signals across axes."""

    run_date: str
    repo_path: str
    in_sync: bool = True
    commit_gap_ahead: int = 0           # local commits NOT on origin
    commit_gap_behind: int = 0          # origin commits NOT local
    working_tree_dirty: bool = False
    modified_files: list[str] = field(default_factory=list)
    untracked_files: list[str] = field(default_factory=list)
    downloads_recent: list[str] = field(default_factory=list)
    spec_divergences: list[dict[str, str]] = field(default_factory=list)
    pytest_collected: int | None = None
    notes: list[str] = field(default_factory=list)


def _run(cmd: list[str], cwd: Path | None = None, check: bool = False) -> tuple[int, str, str]:
    """Run a subprocess; return (exit, stdout, stderr)."""
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(f"{cmd!r} failed: {p.stderr}")
    return p.returncode, p.stdout.strip(), p.stderr.strip()


def check_commit_gap(repo: Path) -> tuple[int, int, list[str]]:
    """Returns (ahead, behind, sample_ahead_messages)."""
    _run(["git", "fetch", "origin", "--quiet"], cwd=repo)
    rc_a, out_a, _ = _run(
        ["git", "rev-list", "--count", "origin/main..HEAD"], cwd=repo
    )
    rc_b, out_b, _ = _run(
        ["git", "rev-list", "--count", "HEAD..origin/main"], cwd=repo
    )
    ahead = int(out_a) if rc_a == 0 and out_a.isdigit() else 0
    behind = int(out_b) if rc_b == 0 and out_b.isdigit() else 0
    sample: list[str] = []
    if ahead > 0:
        _, msgs, _ = _run(
            ["git", "log", "--oneline", "origin/main..HEAD"], cwd=repo
        )
        sample = msgs.splitlines()[:5]
    return ahead, behind, sample


def check_working_tree(repo: Path) -> tuple[bool, list[str], list[str]]:
    """Returns (dirty, modified, untracked)."""
    _, out, _ = _run(["git", "status", "--short"], cwd=repo)
    if not out:
        return False, [], []
    modified: list[str] = []
    untracked: list[str] = []
    for line in out.splitlines():
        if not line.strip():
            continue
        status = line[:2]
        path = line[3:].strip()
        if status.startswith("??"):
            untracked.append(path)
        else:
            modified.append(path)
    return True, modified, untracked


def check_downloads_recent(downloads_dir: Path, lookback_days: int = 3) -> list[str]:
    """List Downloads/ files that look like project artifacts + are recent."""
    if not downloads_dir.exists():
        return []
    date_pattern = re.compile(
        r"(cipro|dna_decoder|decoder|falsifier|leakage|mash)", re.IGNORECASE
    )
    today = _date.today()
    hits: list[str] = []
    for f in downloads_dir.iterdir():
        if not f.is_file():
            continue
        if not date_pattern.search(f.name):
            continue
        if f.suffix.lower() not in {".md", ".json", ".py", ".tab", ".txt"}:
            continue
        # mtime within lookback_days
        try:
            mtime = _date.fromtimestamp(f.stat().st_mtime)
        except Exception:
            continue
        if (today - mtime).days <= lookback_days:
            hits.append(f.name)
    return sorted(hits)


def check_spec_divergence(repo: Path) -> list[dict[str, str]]:
    """For each KNOWN_DIVERGENCE_TARGETS entry, check whether the marker exists in local file."""
    out: list[dict[str, str]] = []
    for relpath, marker in KNOWN_DIVERGENCE_TARGETS:
        f = repo / relpath
        if not f.exists():
            out.append({"file": relpath, "marker": marker, "status": "FILE_MISSING"})
            continue
        try:
            text = f.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            out.append({"file": relpath, "marker": marker, "status": f"READ_ERROR: {e!r}"})
            continue
        if marker in text:
            out.append({"file": relpath, "marker": marker, "status": "PRESENT"})
        else:
            out.append({"file": relpath, "marker": marker, "status": "MISSING_LIKELY_STALE"})
    return out


def check_pytest_collected(repo: Path) -> int | None:
    """Returns the test count from `pytest --collect-only -q`, or None on failure."""
    rc, out, _ = _run(
        ["uv", "run", "pytest", "--collect-only", "-q"], cwd=repo
    )
    if rc != 0 and rc != 5:  # 5 = no tests collected (acceptable info, not failure)
        return None
    # Last non-empty line typically reads "NNN tests collected in X.XXs"
    last_lines = [l for l in out.splitlines() if l.strip()]
    if not last_lines:
        return None
    last = last_lines[-1]
    m = re.search(r"(\d+)\s+tests?\s+collected", last)
    if m:
        return int(m.group(1))
    return None


def render_markdown(report: DriftReport) -> str:
    lines: list[str] = []
    lines.append(f"# Cross-machine sync report — {report.run_date}")
    lines.append("")
    status = "IN SYNC" if report.in_sync else "DRIFT DETECTED"
    lines.append(f"**Status:** {status}")
    lines.append(f"**Repo:** `{report.repo_path}`")
    lines.append("")
    lines.append("## Commit gap")
    lines.append("")
    if report.commit_gap_ahead == 0 and report.commit_gap_behind == 0:
        lines.append("Local and origin in lockstep.")
    else:
        lines.append(f"- Local AHEAD of origin: **{report.commit_gap_ahead}** commits")
        lines.append(f"- Local BEHIND origin: **{report.commit_gap_behind}** commits")
    lines.append("")
    lines.append("## Working tree")
    lines.append("")
    if not report.working_tree_dirty:
        lines.append("Clean.")
    else:
        lines.append(f"- Modified: {len(report.modified_files)} file(s)")
        for f in report.modified_files[:20]:
            lines.append(f"  - {f}")
        if len(report.modified_files) > 20:
            lines.append(f"  - ... and {len(report.modified_files) - 20} more")
        lines.append(f"- Untracked: {len(report.untracked_files)} file(s)")
        for f in report.untracked_files[:20]:
            lines.append(f"  - {f}")
        if len(report.untracked_files) > 20:
            lines.append(f"  - ... and {len(report.untracked_files) - 20} more")
    lines.append("")
    lines.append("## Downloads/ recent project artifacts")
    lines.append("")
    if not report.downloads_recent:
        lines.append("None.")
    else:
        lines.append("Files in Downloads/ that may need to land in the repo:")
        for f in report.downloads_recent:
            lines.append(f"- {f}")
    lines.append("")
    lines.append("## Spec-divergence spot-check")
    lines.append("")
    lines.append("| File | Marker | Status |")
    lines.append("|---|---|---|")
    for d in report.spec_divergences:
        lines.append(f"| `{d['file']}` | `{d['marker']}` | {d['status']} |")
    lines.append("")
    lines.append("## Test-count snapshot")
    lines.append("")
    if report.pytest_collected is None:
        lines.append("Could not collect (pytest may not be installed or failed).")
    else:
        lines.append(f"`pytest --collect-only` reports **{report.pytest_collected}** tests.")
    lines.append("")
    if report.notes:
        lines.append("## Notes")
        lines.append("")
        for n in report.notes:
            lines.append(f"- {n}")
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Cross-machine sync diagnostic.")
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--downloads", type=Path,
                        default=Path.home() / "Downloads",
                        help="Downloads directory to scan for un-landed project artifacts")
    parser.add_argument("--skip-pytest", action="store_true",
                        help="Skip pytest --collect-only step (faster)")
    parser.add_argument("--output-prefix", type=Path,
                        default=Path(f"reports/cross_machine_sync_{_date.today().isoformat()}"))
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    if not (repo / ".git").exists():
        print(f"[sync-check] not a git repo: {repo}", file=sys.stderr)
        return 2

    print(f"[sync-check] inspecting {repo}")

    report = DriftReport(
        run_date=_date.today().isoformat(),
        repo_path=str(repo),
    )

    try:
        ahead, behind, sample = check_commit_gap(repo)
        report.commit_gap_ahead = ahead
        report.commit_gap_behind = behind
        if sample:
            report.notes.append(f"Local-ahead commits (sample): {sample[0]}")
        if ahead > 0 or behind > 0:
            report.in_sync = False

        dirty, modified, untracked = check_working_tree(repo)
        report.working_tree_dirty = dirty
        report.modified_files = modified
        report.untracked_files = untracked
        # Heuristic: many untracked files in `Downloads/` style suggests an active transfer in flight
        if dirty and len(modified) > 0:
            report.in_sync = False

        report.downloads_recent = check_downloads_recent(args.downloads)
        if report.downloads_recent:
            report.notes.append(
                f"Downloads/ has {len(report.downloads_recent)} recent project artifact(s) "
                f"— may need to land in repo."
            )
            report.in_sync = False

        report.spec_divergences = check_spec_divergence(repo)
        if any(d["status"].startswith("MISSING") for d in report.spec_divergences):
            report.in_sync = False

        if not args.skip_pytest:
            print("[sync-check] running pytest --collect-only ...")
            report.pytest_collected = check_pytest_collected(repo)

    except Exception as e:
        print(f"[sync-check] tool error: {e!r}", file=sys.stderr)
        return 2

    # Write outputs
    out_md = args.output_prefix.with_suffix(".md")
    out_json = args.output_prefix.with_suffix(".json")
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(render_markdown(report), encoding="utf-8")
    out_json.write_text(json.dumps(asdict(report), indent=2), encoding="utf-8")

    # Console summary
    status = "IN SYNC" if report.in_sync else "DRIFT DETECTED"
    print(f"[sync-check] {status}")
    print(f"  commit-gap ahead={report.commit_gap_ahead} behind={report.commit_gap_behind}")
    print(f"  working tree {'DIRTY' if report.working_tree_dirty else 'clean'} "
          f"(modified={len(report.modified_files)}, untracked={len(report.untracked_files)})")
    print(f"  Downloads/ recent artifacts: {len(report.downloads_recent)}")
    miss = [d for d in report.spec_divergences if d["status"].startswith("MISSING")]
    print(f"  spec-divergence missing markers: {len(miss)}/{len(report.spec_divergences)}")
    if report.pytest_collected is not None:
        print(f"  tests collected: {report.pytest_collected}")
    print(f"  wrote {out_md}")
    print(f"  wrote {out_json}")

    return 0 if report.in_sync else 1


if __name__ == "__main__":
    sys.exit(main())
