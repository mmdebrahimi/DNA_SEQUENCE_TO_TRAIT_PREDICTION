"""The orientation surface must not go stale silently.

CLAUDE.md's first section is the ONLY project text auto-loaded into every session, so an error there
propagates into every future run's working model. On 2026-08-29 its opening described a "Phase 1 E. coli
platform predicting cipro/cef/tet" -- three months and six tracks out of date -- and a whole session
reasoned from it, understating the tool's evidence surface ~4x.

The orientation block quotes four figures for concreteness. Quoting them makes them drift-able, so these
tests pin each one to the LIVE registry: when the registry moves, the test fails and names the line to
update. A loud failure beats a silent lie in the file every session reads.

Offline; no network, no Docker, no model load.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CLAUDE_MD = ROOT / "CLAUDE.md"


@pytest.fixture(scope="module")
def orientation() -> str:
    """The READ THIS FIRST block only -- up to the next top-level section."""
    txt = CLAUDE_MD.read_text(encoding="utf-8")
    start = txt.index("## READ THIS FIRST")
    nxt = txt.index("\n## ", start + 10)
    return txt[start:nxt]


def test_project_status_script_runs_clean():
    """Orientation points at this script as the authority; it must actually run."""
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "project_status.py")],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    assert "110" not in r.stdout or "registered cells" in r.stdout


def test_orientation_cell_count_matches_live_registry(orientation):
    from dna_decode.data.cell_registry import cells

    m = re.search(r"\*\*(\d+)\s+cells,\s+(\d+)\s+`INDEPENDENT_MEASURED`\*\*", orientation)
    assert m, "the orientation block no longer states the registry surface in the pinned form"
    quoted_total, quoted_indep = int(m.group(1)), int(m.group(2))

    cs = cells()
    live_indep = sum(1 for c in cs if c.evidence_tier.value == "independent_measured")
    assert quoted_total == len(cs), (
        f"CLAUDE.md orientation says {quoted_total} cells; registry has {len(cs)}. Update the line.")
    assert quoted_indep == live_indep, (
        f"CLAUDE.md says {quoted_indep} INDEPENDENT_MEASURED; registry has {live_indep}. Update the line.")


def test_orientation_report_card_scope_matches_the_card(orientation):
    """The card's row/SCORED counts are quoted to SHOW the scope gap -- so they must be the card's own."""
    import json

    f = ROOT / "wiki" / "decoder_validation_report_card.json"
    if not f.exists():
        pytest.skip("report card artifact absent")
    d = json.loads(f.read_text(encoding="utf-8"))
    rows = d.get("cells") or d.get("rows") or []
    scored = sum(1 for r in rows if (r.get("cell_state") or r.get("state")) == "SCORED")

    m = re.search(r"\*\*(\d+) rows / (\d+) SCORED\*\*", orientation)
    assert m, "the orientation block no longer states the AMR card scope in the pinned form"
    assert (int(m.group(1)), int(m.group(2))) == (len(rows), scored), (
        f"CLAUDE.md says {m.group(1)} rows / {m.group(2)} SCORED; card has {len(rows)} / {scored}.")


def test_orientation_does_not_describe_the_tool_as_an_ecoli_amr_project(orientation):
    """The exact regression. The old opening called this a 3-drug E. coli platform for months.

    Guard the CLAIM, not a phrase: the block must state the multi-track reality, and must not present
    the drug triple as what the tool predicts.
    """
    assert "not an E. coli AMR project" in orientation.lower() or \
           "not an e. coli amr project" in orientation.lower()
    assert "multi-kingdom" in orientation
    stale = "Predicts antibiotic resistance (ciprofloxacin / ceftriaxone / tetracycline) from genomic DNA"
    assert stale not in CLAUDE_MD.read_text(encoding="utf-8"), \
        "the superseded Phase-1 framing is back at the top of CLAUDE.md"


def test_orientation_names_the_population_design_finding(orientation):
    """Error 2 was compressing a scoped negative three separate times. Pin the correction."""
    assert "population design" in orientation
    assert "12/12" in orientation, "the yeast-cross positive is what makes the correction checkable"


def test_orientation_stays_short_enough_to_actually_be_read(orientation):
    """A 300-line preamble is not orientation. If it grows, move detail into the body below."""
    assert len(orientation.splitlines()) < 45, "orientation block is growing into a second document"
