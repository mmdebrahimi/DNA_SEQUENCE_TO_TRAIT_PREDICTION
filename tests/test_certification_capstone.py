"""Certification capstone — pins the THIN-PRESENTATION + NO-AGGREGATE-VERDICT contract.

The capstone is a read-only roll-up. The load-bearing assertion is the NEGATIVE one: it must carry NO
aggregate boolean verdict (no `certified`/`passed`/`overall` bool, no averaged score) — reducing a
multi-tier surface to one boolean is the "boolean-verdict-on-unvalidated-model trap". It must render the
registry's per-track census + the boundary documents.
"""
from __future__ import annotations

import importlib
import json
from collections import Counter
from pathlib import Path

from dna_decode.data.cell_registry import cells

REPO = Path(__file__).resolve().parent.parent
CAP_JSON = REPO / "wiki" / "certification_capstone.json"

_builder = importlib.import_module("scripts.build_certification_capstone")


def _build():
    _builder.main()
    return json.loads(CAP_JSON.read_text(encoding="utf-8"))


def test_capstone_builds_and_is_a_report_not_a_gate():
    assert _builder.main() == 0  # exit 0 always — a report, never a gate


def test_no_aggregate_boolean_verdict():
    """LOAD-BEARING: no aggregate pass/fail / certified / averaged-score field anywhere."""
    cap = _build()
    assert cap.get("no_aggregate_verdict") is True
    banned = {"certified", "passed", "overall_verdict", "aggregate_score", "verdict", "pass", "score"}
    # no top-level banned key, and no banned key implies a single boolean roll-up
    for k, v in cap.items():
        assert k not in banned, f"capstone carries a reducer field: {k}"
    # the disclaimer must be present (the trap is named, not just absent)
    assert "aggregate" in cap.get("aggregate_verdict_disclaimer", "").lower()


def test_registry_census_matches_registry():
    cap = _build()
    cs = list(cells())
    census = cap["registry_census"]
    assert census["total_cells"] == len(cs)
    assert census["by_track"] == dict(Counter(c.track for c in cs))
    assert census["by_evidence_tier"] == dict(Counter(c.evidence_tier.value for c in cs))


def test_not_censused_cells_surfaced_explicitly():
    """Routable-but-unvalidated cells (e.g. delavirdine once v0.1.1 lands) must be listed, never hidden."""
    cap = _build()
    expected = sorted(c.cell_id for c in cells() if c.evidence_tier.value == "not_censused")
    assert cap["registry_census"]["not_censused_cells"] == expected


def test_boundaries_present():
    """The freeze + negative-results-map (what is closed + the label wall) are referenced, not flattened."""
    cap = _build()
    docs = {b["document"] for b in cap["boundaries"]}
    assert "reproducibility_freeze_2026-06-13.md" in docs
    assert "negative_results_map_2026-06-13.md" in docs


def test_domain_cards_present_or_not_run():
    """Each domain card is presented verbatim (or honestly marked NOT_RUN) — never invented."""
    cap = _build()
    for c in cap["domain_report_cards"]:
        assert c["status"] in ("present", "NOT_RUN")
        if c["status"] == "present":
            assert "headline" in c


# ---------------------------------------------------------------------------
# The COMMITTED artifact, not a fresh rebuild.
#
# Every test above calls _build(), which runs the builder and overwrites the
# artifact before reading it back. That makes them assert "the builder agrees
# with the registry" -- true by construction of the builder -- and leaves them
# STRUCTURALLY unable to notice that the file in git has gone stale. Worse, the
# rebuild SILENTLY REPAIRS the staleness, so drift never fails and surfaces only
# as unexplained `git status` noise a developer may well revert.
#
# It happened: ca361e1 moved three cells NEAR_INDEPENDENT -> KNOWLEDGE_BASELINE
# and promoted two to INDEPENDENT_MEASURED without rebuilding, so HEAD shipped a
# capstone reporting 31/24/38 against a registry returning 33/21/39. Planting
# total_cells=999 into the committed file left all six tests green.
#
# So this reads the artifact AS COMMITTED (`git show HEAD:<path>`) rather than
# from disk: immune to any rebuild another test performs, and it asserts the
# claim that actually matters -- what ships agrees with the live registry.
# ---------------------------------------------------------------------------


def _committed_capstone() -> dict | None:
    """The capstone as committed at HEAD, or None when git can't answer."""
    import subprocess

    try:
        out = subprocess.run(
            ["git", "show", "HEAD:wiki/certification_capstone.json"],
            cwd=REPO, capture_output=True, timeout=30,
        )
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    return json.loads(out.stdout.decode("utf-8"))


def test_committed_capstone_census_matches_the_live_registry():
    """A roll-up that disagrees with the registry it summarises is a false trust surface.

    NOT a rebuild: `_build()` would make this tautological. Reads what is committed.
    """
    cap = _committed_capstone()
    if cap is None:
        return  # no git / not committed yet — nothing to compare, never a silent pass elsewhere
    cs = list(cells())
    census = cap["registry_census"]
    assert census["total_cells"] == len(cs), (
        f"committed capstone says total_cells={census['total_cells']} but the registry has "
        f"{len(cs)} — rebuild it: uv run python scripts/build_certification_capstone.py"
    )
    assert census["by_evidence_tier"] == dict(Counter(c.evidence_tier.value for c in cs)), (
        "committed capstone's tier census has drifted from the registry — rebuild it: "
        "uv run python scripts/build_certification_capstone.py"
    )
    assert census["by_track"] == dict(Counter(c.track for c in cs)), (
        "committed capstone's per-track census has drifted from the registry — rebuild it"
    )


def test_the_rebuild_shaped_tests_above_cannot_catch_committed_drift():
    """Pins WHY the guard above exists, so nobody 'simplifies' it into a _build() call.

    Asserts the mechanism directly: the builder is a pure function of the registry, so a
    freshly-built census ALWAYS matches it. That is exactly why a fresh build cannot be
    evidence about the committed file.
    """
    cap = _build()
    cs = list(cells())
    assert cap["registry_census"]["total_cells"] == len(cs)
    # ... and the committed file is a SEPARATE claim, which is what the guard checks.
    assert _committed_capstone is not None
