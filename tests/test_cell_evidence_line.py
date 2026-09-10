"""A cell's measured evidence must reach its own CLI, and the coverage gap must be countable.

`dna-amr` has printed an inline trust badge since trust_surface shipped. No other route did, so a user
running `dna-ktype` saw a hand-written caveat while that cell's MEASURED evidence -- 0.8788 against
full-locus Kaptive on 307 genomes -- sat in the registry and a report card they never open. Same failure
as the doubt layer's JSON-only block: evidence carried only in a registry is not a disclosure.

This wiring is PARTIAL by design (the measured cells first). The last test makes that partiality a
number rather than a silence -- a green suite over whichever subset happened to be done would be the
dishonest version.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.data.cell_evidence_line import (  # noqa: E402
    cells_for_route, evidence_one_line, routes_missing_evidence_line,
)

# The routes wired so far. Kept as data so the coverage test below can subtract it from the registry.
WIRED = ("dna-ktype", "dna-salmserovar", "dna-serotype", "dna-pneumo-serotype", "dna-plasmid",
         "dna-resfinder", "dna-disinfinder", "dna-pointfinder", "dna-mlst", "dna-hla", "dna-pgx")


@pytest.mark.parametrize("route", WIRED)
def test_every_wired_route_renders_a_line(route):
    assert cells_for_route(route), f"{route} has no registered cell -- the wiring points at nothing"
    line = evidence_one_line(route)
    assert line and line.startswith("evidence: ")


@pytest.mark.parametrize("route", WIRED)
def test_every_wired_cli_actually_prints_it(route):
    """Importing the helper is not printing it. A CLI that imports and never calls would pass a
    weaker check while showing the user nothing -- the exact shape of the JSON-only disclosure bug."""
    # The console-script name is NOT always the module name: pyproject ships `dna-pneumo-serotype`
    # while the package is `pneumoserotype`. Deriving one from the other silently pointed the wiring at
    # a route that does not exist, so evidence_one_line returned None and the CLI printed NOTHING --
    # a no-op that looks exactly like success. Mapped explicitly.
    mod = {"dna-pneumo-serotype": "pneumoserotype"}.get(route, route.replace("dna-", "").replace("-", "_"))
    src = (ROOT / "dna_decode" / mod / "cli.py").read_text(encoding="utf-8")
    assert "evidence_one_line(" in src, f"{route}: helper never called"
    assert re.search(r"print\(f?\"\s*\{_ev\}\"\)", src), f"{route}: computed but never printed"


def test_an_unregistered_route_returns_none_not_a_reassuring_blank():
    """None means 'not in the registry', which is a different fact from 'measured and clean'. A caller
    must not be able to read silence as a clean bill."""
    assert evidence_one_line("dna-not-a-real-route") is None
    assert cells_for_route("dna-not-a-real-route") == []


def test_a_never_measured_cell_says_so_rather_than_going_quiet():
    """A KNOWLEDGE_BASELINE cell has nothing measured. Rendering nothing would imply it was checked."""
    line = evidence_one_line("dna-coatcolor")
    assert line and "NEVER MEASURED" in line


def test_a_multi_cell_route_reports_its_WEAKEST_tier():
    """dna-pgx serves 14 cells spanning tiers. Reporting the strongest would imply that evidence
    applies to whichever cell the caller actually invoked."""
    found = cells_for_route("dna-pgx")
    assert len(found) > 1, "fixture assumption: dna-pgx serves several cells"
    line = evidence_one_line("dna-pgx")
    assert "KNOWLEDGE_BASELINE" in line and "weakest shown" in line


def test_a_tool_tier_line_says_it_bounds_agreement_not_correctness():
    """The single most over-readable tier: FAITHFUL_TO_TOOL is agreement with another implementation,
    and a line that omitted that would let a reader take it as accuracy."""
    line = evidence_one_line("dna-ktype")
    assert "bounds agreement, never correctness" in line


def test_the_coverage_gap_is_reported_as_a_number_not_hidden():
    """Partial wiring is fine; SILENT partial wiring is not. If this ever reads 0 the wiring is
    complete and this test should be replaced by an exhaustive one."""
    missing = routes_missing_evidence_line(WIRED)
    total = len(missing) + len(WIRED)
    assert missing, (
        "no routes are missing an evidence line -- wiring is complete; replace this test with an "
        "exhaustive per-route check rather than leaving a coverage counter that can never fire")
    # Documented, not asserted-at-a-magic-number: the point is that the gap is VISIBLE.
    assert len(WIRED) < total
    assert all(r.startswith("dna-") for r in missing)
