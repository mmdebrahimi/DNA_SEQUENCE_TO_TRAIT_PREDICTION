"""Tests for the live-vs-sg8 delta diagnostic and for the in-place git init.

The classifier tests are pure and reference neither real tree. The git tests assert
the state Step 9 of the reconciliation plan establishes, and skip (with a reason)
where that init has not happened -- they are assertions about this host, not about
every checkout.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.live_vs_sg8_delta import (
    CONTENT_DIFFERS,
    IDENTICAL,
    LIVE_ONLY,
    SG8_ONLY,
    classify,
    index_tree,
    looks_like_text,
)

REPO = Path(__file__).resolve().parent.parent
SG8 = Path(r"C:/Users/Laptop Portable/sg8")

# Commands that must never appear in this step's scripts. Each would destroy live
# artifacts that exist in no other copy at their current state.
FORBIDDEN = ("reset --hard", "clean -fd", "checkout .")
STEP9_SCRIPTS = ("scripts/live_vs_sg8_delta.py",)


# --------------------------------------------------------------------------- #
# the pure classifier
# --------------------------------------------------------------------------- #

def test_classify_covers_all_four_states():
    live = {"both_same.py": "aaa", "both_diff.py": "bbb", "only_live.py": "ccc"}
    sg8 = {"both_same.py": "aaa", "both_diff.py": "zzz", "only_sg8.py": "ddd"}
    delta = classify(live, sg8)
    assert delta.as_dict() == {
        LIVE_ONLY: ["only_live.py"],
        SG8_ONLY: ["only_sg8.py"],
        CONTENT_DIFFERS: ["both_diff.py"],
        IDENTICAL: ["both_same.py"],
    }


def test_classify_counts_partition_the_union():
    live = {f"l{i}": str(i) for i in range(5)}
    sg8 = {f"l{i}": str(i) for i in range(3)} | {"s0": "x"}
    delta = classify(live, sg8)
    counts = delta.counts()
    assert sum(counts.values()) == len(set(live) | set(sg8))


def test_classify_on_two_empty_trees():
    delta = classify({}, {})
    assert delta.counts() == {LIVE_ONLY: 0, SG8_ONLY: 0, CONTENT_DIFFERS: 0, IDENTICAL: 0}


def test_classify_is_not_order_sensitive():
    live = {"a": "1", "b": "2"}
    sg8 = {"b": "2", "a": "9"}
    assert classify(live, sg8).as_dict()[CONTENT_DIFFERS] == ["a"]


# --------------------------------------------------------------------------- #
# the EOL-normalised second index
# --------------------------------------------------------------------------- #

def test_text_detection_is_nul_based():
    assert looks_like_text(b"def f():\n    return 1\n")
    assert not looks_like_text(b"\x89PNG\r\n\x1a\n\x00\x00")


def test_index_tree_normalises_crlf_for_text_only(tmp_path):
    crlf = tmp_path / "a.txt"
    crlf.write_bytes(b"one\r\ntwo\r\n")
    lf = tmp_path / "b.txt"
    lf.write_bytes(b"one\ntwo\n")
    raw, norm = index_tree(tmp_path)
    assert raw["a.txt"] != raw["b.txt"], "raw hashes must see the line endings"
    assert norm["a.txt"] == norm["b.txt"], "normalised hashes must not"


def test_binary_is_never_eol_normalised(tmp_path):
    """Collapsing CRLF inside a binary could manufacture false identity."""
    blob = b"\x00\x01\r\n\x02"
    p = tmp_path / "x.bin"
    p.write_bytes(blob)
    raw, norm = index_tree(tmp_path)
    assert raw["x.bin"] == norm["x.bin"] == hashlib.sha256(blob).hexdigest()


def test_skip_dirs_are_honoured(tmp_path):
    (tmp_path / "keep").mkdir()
    (tmp_path / "keep" / "f.py").write_text("x", encoding="utf-8")
    for skipped in ("__pycache__", "data", ".git"):
        (tmp_path / skipped).mkdir()
        (tmp_path / skipped / "f.py").write_text("x", encoding="utf-8")
    raw, _ = index_tree(tmp_path)
    assert set(raw) == {"keep/f.py"}


def test_absent_root_yields_empty_indexes(tmp_path):
    raw, norm = index_tree(tmp_path / "nope")
    assert raw == {} and norm == {}


# --------------------------------------------------------------------------- #
# the diagnostic writes only its own artifact
# --------------------------------------------------------------------------- #

def test_diagnostic_writes_nothing_into_either_tree(tmp_path):
    """'Read-only' as a MEASURED property -- on synthetic trees, so it is airtight.

    This REPLACES a whole-tree byte-snapshot of the real repo, which was flaky by
    construction: it compared ~3,400 live files around a 30 s operation, so any
    concurrent writer failed it with "the live tree was modified". That is a FALSE
    ACCUSATION about the exact safety property the test exists to defend, and a test
    that cries wolf about its own subject trains the reader to discount it. It was
    observed firing for real on this host: a line-ending normaliser rewrote
    `scripts/live_vs_sg8_delta.py` and `wiki/live_vs_sg8_delta_2026-10-19.json`
    mid-run, and the diagnostic had touched neither.

    Synthetic trees have no third-party writer, so byte-identity is fair here -- and
    this runs on hosts with no sg8 clone, closing a skip.
    """
    from scripts.live_vs_sg8_delta import main

    live, sg8 = tmp_path / "live", tmp_path / "sg8"
    for root in (live, sg8):
        (root / "pkg").mkdir(parents=True)
        (root / "pkg" / "mod.py").write_text("x = 1\n", encoding="utf-8")
    (live / "pkg" / "only_live.py").write_text("y = 2\n", encoding="utf-8")
    (sg8 / "pkg" / "mod.py").write_text("x = 2\n", encoding="utf-8")

    before = (index_tree(live), index_tree(sg8))
    out = tmp_path / "delta.json"          # deliberately OUTSIDE both trees
    argv = ["--live", str(live), "--sg8", str(sg8), "--out", str(out), "--quiet"]
    assert main(argv) == 0
    assert main(argv) == 0, "a second run must also be read-only"
    assert (index_tree(live), index_tree(sg8)) == before, "a tree was modified"

    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["_schema"]
    assert payload["counts"][LIVE_ONLY] == 1
    assert payload["counts"][SG8_ONLY] == 0
    assert payload["counts"][CONTENT_DIFFERS] == 1


@pytest.mark.skipif(not SG8.is_dir(), reason="sg8 clone not present on this host")
def test_the_real_run_creates_and_deletes_no_path_in_either_tree(tmp_path):
    """The real-tree half, scoped to invariants a concurrent writer cannot forge.

    Byte-identity of the LIVE tree is deliberately not asserted -- it is actively
    edited on this host. What is asserted is that the diagnostic creates no file and
    removes none from either tree, and that the sg8 clone (quiescent; nothing has it
    open) stays byte-identical. Deletion is the destructive direction the FORBIDDEN
    command list guards, and creation is how a stray write would show up.
    """
    from scripts.live_vs_sg8_delta import main

    before_live, _ = index_tree(REPO)
    before_sg8, _ = index_tree(SG8)
    assert main(["--out", str(tmp_path / "delta.json"), "--quiet"]) == 0
    after_live, _ = index_tree(REPO)
    after_sg8, _ = index_tree(SG8)

    assert set(after_live) == set(before_live), (
        f"paths created/removed under live: {set(after_live) ^ set(before_live)}"
    )
    assert after_sg8 == before_sg8, "the sg8 clone was modified"


@pytest.mark.skipif(not SG8.is_dir(), reason="sg8 clone not present on this host")
def test_diagnostic_is_deterministic(tmp_path):
    from scripts.live_vs_sg8_delta import main

    a, b = tmp_path / "a.json", tmp_path / "b.json"
    main(["--out", str(a), "--quiet"])
    main(["--out", str(b), "--quiet"])
    da = json.loads(a.read_text(encoding="utf-8"))
    db = json.loads(b.read_text(encoding="utf-8"))
    assert da["delta"] == db["delta"]


def test_missing_sg8_exits_2_rather_than_reporting_an_empty_delta(tmp_path):
    """An absent clone must not read as 'everything is live-only'."""
    from scripts.live_vs_sg8_delta import main

    assert main(["--sg8", str(tmp_path / "absent"), "--out",
                 str(tmp_path / "o.json"), "--quiet"]) == 2


# --------------------------------------------------------------------------- #
# the forbidden-command guard
# --------------------------------------------------------------------------- #

def test_no_destructive_git_command_appears_in_the_step_scripts():
    for rel in STEP9_SCRIPTS:
        text = (REPO / rel).read_text(encoding="utf-8")
        for bad in FORBIDDEN:
            assert bad not in text, f"{rel} contains the forbidden {bad!r}"


# --------------------------------------------------------------------------- #
# the in-place git init
# --------------------------------------------------------------------------- #

def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)


def _is_repo() -> bool:
    return (REPO / ".git").exists()


needs_init = pytest.mark.skipif(
    not _is_repo(), reason="the project has not been git-init'd in place on this host"
)


@needs_init
def test_toplevel_now_resolves_to_the_project_not_the_home_directory():
    """The anchoring fix.

    Before the in-place init, `git rev-parse --show-toplevel` returned
    `C:/Users/Laptop Portable`, whose remote is an unrelated repository -- which is
    how a dna_decode plan entry nearly got written into another project's index.
    """
    proc = _git("rev-parse", "--show-toplevel")
    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()).resolve() == REPO.resolve()


@needs_init
def test_the_four_soraya_runs_negations_are_re_included():
    """Asserted by comment in .gitignore for months; never checked by git until now.

    `--no-index` IS LOAD-BEARING. `git check-ignore` consults the index by default
    and reports any TRACKED path as not-ignored, so without this flag the assertion
    could pass because a file happens to be tracked in origin/main rather than
    because the negation works -- which is precisely a vacuous test. Measured while
    writing this: `soraya_runs/self-init-potential.json` reads un-ignored WITH the
    index (it is tracked) and ignored WITHOUT it (it carries no negation).
    """
    # Exactly the four `!soraya_runs/...` negations the live .gitignore carries.
    # `constraint-expiry.json` is the load-bearing one: sg8's .gitignore has only
    # three and is missing precisely that negation, which is the measured reason
    # this project was git-init'd in place rather than committed under sg8's rules.
    for name in (
        "soraya_runs/standing-positions.json",
        "soraya_runs/negative-results.json",
        "soraya_runs/arms-ledger.json",
        "soraya_runs/constraint-expiry.json",
    ):
        assert (REPO / name).exists(), f"{name} should exist to be worth negating"
        proc = _git("check-ignore", "--no-index", "-q", name)
        assert proc.returncode != 0, (
            f"{name} is IGNORED -- the trailing-/* re-inclusion is not working, and "
            "committing under these rules would silently drop it"
        )


@needs_init
def test_an_unnegated_soraya_runs_sidecar_is_still_ignored():
    """The negation list is an allowlist, not a blanket un-ignore of the directory.

    `continuity.json` and `self-init-potential.json` sit beside the four negated
    sidecars and are NOT negated, so the rules must still exclude them. Without this
    the re-inclusion test above would pass just as well if the `soraya_runs/*` rule
    had simply stopped applying.
    """
    for name in ("soraya_runs/continuity.json", "soraya_runs/self-init-potential.json"):
        if not (REPO / name).exists():
            continue
        proc = _git("check-ignore", "--no-index", "-q", name)
        assert proc.returncode == 0, f"{name} is not negated, so it should be ignored"


@needs_init
def test_origin_tracks_paths_the_current_rules_would_exclude():
    """A recorded fact, not a defect: ignore rules do not apply to tracked files.

    origin/main tracks dated `soraya_runs/<run-id>/` artifacts and
    `self-init-potential.json`, all of which the live `.gitignore` now excludes. So
    those paths stay in history and in `git status` regardless of the ignore rules,
    and anyone reading `check-ignore` output without `--no-index` will get a
    misleading answer for exactly these paths.
    """
    tracked = subprocess.run(
        ["git", "ls-files", "soraya_runs/"], cwd=REPO, capture_output=True, text=True
    ).stdout.split()
    dated = [p for p in tracked if p.startswith("soraya_runs/20")]
    assert dated, "expected origin/main to track dated run directories"
    # Each is ignored at the RULE level yet tracked -- that is the whole point.
    sample = dated[0]
    assert _git("check-ignore", "--no-index", "-q", sample).returncode == 0
    assert _git("check-ignore", "-q", sample).returncode != 0


@needs_init
def test_the_negation_block_is_non_vacuous():
    """The other half of the pair: things that MUST stay ignored.

    A re-inclusion test alone would pass just as well if nothing were ignored at all.
    """
    for name in ("data/processed/x.parquet", "data/raw/y.h5"):
        proc = _git("check-ignore", "-q", name)
        assert proc.returncode == 0, f"{name} should be ignored but is not"


@needs_init
def test_a_dated_run_directory_is_ignored_while_its_sidecars_are_not():
    proc = _git("check-ignore", "-q", "soraya_runs/2026-10-19-1200-some-run/result.md")
    assert proc.returncode == 0, "a dated run directory should be ignored"


@needs_init
def test_the_remote_points_at_the_project_repo_not_the_home_repo():
    proc = _git("remote", "get-url", "origin")
    if proc.returncode != 0:
        pytest.skip("no origin configured yet")
    url = proc.stdout.strip()
    assert "DNA_SEQUENCE_TO_TRAIT_PREDICTION" in url, url
    assert "Guided-Capstone" not in url, (
        "origin points at the home directory's unrelated repository"
    )
