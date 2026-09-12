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
WIRED = (
    # typing / finder cells (first wave)
    "dna-ktype", "dna-salmserovar", "dna-serotype", "dna-pneumo-serotype", "dna-plasmid",
    "dna-resfinder", "dna-disinfinder", "dna-pointfinder", "dna-mlst", "dna-hla", "dna-pgx",
    # the 14 mammal coat-colour routes -- ONE call site in mammal_color_cli, route name DERIVED
    # (f"dna-{organism}color") rather than hand-listed, so a new organism cannot ship unwired.
    "dna-rabbitcolor", "dna-mousecolor", "dna-cattlecolor", "dna-pigcolor", "dna-sheepcolor",
    "dna-goatcolor", "dna-alpacacolor", "dna-guineapigcolor", "dna-foxcolor", "dna-donkeycolor",
    "dna-buffalocolor", "dna-camelcolor", "dna-minkcolor", "dna-roedeercolor",
    # remaining single-module routes
    "dna-catcolor", "dna-coatcolor", "dna-horsecolor", "dna-pigeoncolor", "dna-plumage",
    "dna-pigment", "dna-morphology", "dna-clinvar", "dna-pathotype", "dna-phage",
    "dna-essentiality", "dna-fba", "dna-metabolic", "dna-motility", "dna-kleb", "dna-flowering",
)

# DELIBERATELY NOT WIRED, and this is the interesting half. Each of these already prints a MORE
# SPECIFIC evidence disclosure than the generic route-level line, so adding it would be redundant at
# best and SELF-CONTRADICTING at worst. Measured, not assumed:
#
#   dna-amr             prints a per-drug/organism `validation:` line from trust_block -- e.g.
#                       "INDEPENDENT_MEASURED -- acc 0.95 (N=16007)" for ciprofloxacin -- plus four
#                       disclosure layers. The generic line would append "NOT_CENSUSED (59 cells,
#                       weakest shown) -- shipped but never censused" DIRECTLY BENEATH it. That is
#                       the `lenacapavir` defect exactly: a human-facing disclosure contradicting
#                       itself two lines apart.
#   dna-decode-forward  prints method-specific measured evidence (BLOSUM62 vs a learned method, with
#                       the measured comparison between them).
#   dna-decode-inverse  prints its own "evidence: ESM (learned) beats a no-oracle null on ..." line.
#                       Both are `dna-decode` SUBCOMMANDS, not console scripts, so they also have no
#                       standalone CLI module of the wired shape.
EXCLUDED_ALREADY_DISCLOSE_BETTER = {
    "dna-amr": "per-drug trust_block line + 4 disclosure layers; generic line would contradict it",
    "dna-decode-forward": "prints method-specific measured evidence; a dna-decode subcommand",
    "dna-decode-inverse": "prints its own evidence line; a dna-decode subcommand",
}


@pytest.mark.parametrize("route", WIRED)
def test_every_wired_route_renders_a_line(route):
    assert cells_for_route(route), f"{route} has no registered cell -- the wiring points at nothing"
    line = evidence_one_line(route)
    assert line and line.startswith("evidence: ")


def _module_path_for(route: str) -> Path:
    """Resolve a console-script route to the file that implements it, DERIVED from pyproject.

    Never guessed from the route name. `dna-pneumo-serotype` ships from the package `pneumoserotype`,
    and deriving one from the other once pointed the wiring at a route that does not exist --
    `evidence_one_line` returned None and the CLI printed NOTHING, a no-op indistinguishable from
    success. pyproject is the only authority on which module a script actually runs.
    """
    import tomllib
    scripts = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["scripts"]
    ep = scripts.get(route)
    assert ep, f"{route} is not a console script in pyproject"
    return ROOT / (ep.split(":")[0].replace(".", "/") + ".py")


@pytest.mark.parametrize("route", WIRED)
def test_every_wired_cli_actually_prints_it(route):
    """Importing the helper is not printing it. A CLI that imports and never calls would pass a
    weaker check while showing the user nothing -- the exact shape of the JSON-only disclosure bug."""
    src = _module_path_for(route).read_text(encoding="utf-8")
    assert "evidence_one_line(" in src, f"{route}: helper never called"
    assert re.search(r"print\(f?\"\s*\{_ev\}\"\)", src), f"{route}: computed but never printed"


@pytest.mark.parametrize("route", sorted(EXCLUDED_ALREADY_DISCLOSE_BETTER))
def test_an_excluded_route_really_does_disclose_on_its_own(route):
    """NON-VACUITY for the exclusions. 'It already discloses' is the entire justification for leaving
    these three unwired, so if one ever stops disclosing, the exclusion becomes a silent gap and this
    must fail rather than sit quietly in a comment."""
    if route == "dna-amr":
        src = (ROOT / "dna_decode" / "amr" / "cli.py").read_text(encoding="utf-8")
        assert "trust_block(" in src and "one_line(" in src
        return
    name = route.replace("dna-decode-", "")
    src = (ROOT / "dna_decode" / "forward" / ("cli.py" if name == "forward" else "inverse_cli.py")
           ).read_text(encoding="utf-8")
    assert re.search(r"print\(.*(evidence|measured to beat|BLOSUM62)", src), f"{route}: no own disclosure"


def test_the_generic_line_would_contradict_the_amr_route():
    """The measured reason dna-amr is excluded, pinned so the exclusion cannot decay into folklore.
    Its per-drug line says INDEPENDENT_MEASURED with a real accuracy; the route-level line reports the
    WEAKEST of 59 cells, which is NOT_CENSUSED. Printing both is a self-contradicting disclosure."""
    from dna_decode.data.trust_surface import one_line as trust_one_line, trust_block
    per_drug = trust_one_line(trust_block("ciprofloxacin", "Escherichia_coli_Shigella")) or ""
    route_level = evidence_one_line("dna-amr") or ""
    assert "INDEPENDENT_MEASURED" in per_drug
    assert "NOT_CENSUSED" in route_level
    assert "weakest shown" in route_level


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


def test_the_only_unwired_routes_are_the_ones_deliberately_excluded():
    """The wiring is now EXHAUSTIVE, so this asserts an exact set rather than 'a gap exists'.

    The old form only checked that SOME routes were missing, which a new unwired route would have
    satisfied silently -- the coverage counter could never distinguish 'known exclusion' from 'someone
    shipped a route and forgot'. Now a new route fails here by name until it is either wired or given
    an explicit exclusion reason.
    """
    missing = set(routes_missing_evidence_line(WIRED))
    assert missing == set(EXCLUDED_ALREADY_DISCLOSE_BETTER), (
        f"unwired routes changed.\n"
        f"  newly unwired (wire them, or add an exclusion reason): "
        f"{sorted(missing - set(EXCLUDED_ALREADY_DISCLOSE_BETTER))}\n"
        f"  no longer unwired (drop them from EXCLUDED): "
        f"{sorted(set(EXCLUDED_ALREADY_DISCLOSE_BETTER) - missing)}")
    assert every_reason_is_recorded()


def every_reason_is_recorded() -> bool:
    """An exclusion without a stated reason is indistinguishable from an oversight."""
    return all(isinstance(v, str) and len(v) > 20
               for v in EXCLUDED_ALREADY_DISCLOSE_BETTER.values())


def test_the_wired_set_matches_the_registry_and_holds_no_phantoms():
    """Every WIRED entry must be a real registry route. A typo would silently shrink the coverage
    assertion above (the typo'd route stays in `missing`, the real one looks wired)."""
    for route in WIRED:
        assert cells_for_route(route), f"{route} is in WIRED but has no registered cell"
