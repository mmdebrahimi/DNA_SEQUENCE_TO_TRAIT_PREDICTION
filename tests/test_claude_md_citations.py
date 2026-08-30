"""CLAUDE.md's pointers must resolve, and the checker must not invent breakage.

CLAUDE.md is auto-loaded into every session, so a bullet promising `wiki/foo.md` that is not on disk is
worse than no pointer at all: it reads as authoritative provenance and sends a future session looking for
something that does not exist. This caught exactly that -- `plans/Genome_Map_Virulence_Overlay_Plan/`
completed and moved to `executed_plans/`, and the citation did not follow.

The parser tests matter as much as the guard. The first run of this checker reported SIX broken citations
of which FIVE were its own bugs; a checker whose false positives look exactly like the defect it hunts is
worse than not having one. Each of those bugs is pinned below.

Offline; pure path resolution, no network.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from claude_md_weight import analyse, bullets, expand_braces, resolves  # noqa: E402


def test_every_citation_in_claude_md_resolves():
    """The guard. A pointer to a file that does not exist is a false provenance claim."""
    text = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    broken = [(r["head"], m) for r in map(analyse, bullets(text)) for m in r["missing"]]
    assert not broken, "CLAUDE.md cites paths that do not exist:\n" + "\n".join(
        f"  {m}  <- {h[:60]}" for h, m in broken)


def test_the_guard_is_not_vacuous():
    """If the citation regex ever stops matching, the test above passes by finding nothing at all."""
    text = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    total = sum(r["n_cited"] for r in map(analyse, bullets(text)))
    assert total > 40, f"only {total} citations found -- the regex is probably broken, not the file"


def test_brace_expansion_handles_every_group_not_just_the_first():
    """BUG 1: one-group expansion left `.{md,json}` literal in the path, so three real citations were
    reported broken."""
    got = expand_braces("wiki/x_2026-07-1{6,7}.{md,json}")
    assert set(got) == {"wiki/x_2026-07-16.md", "wiki/x_2026-07-16.json",
                        "wiki/x_2026-07-17.md", "wiki/x_2026-07-17.json"}


def test_a_path_with_no_braces_is_returned_unchanged():
    assert expand_braces("wiki/plain.md") == ["wiki/plain.md"]


def test_template_placeholder_resolves_by_prefix():
    """BUG 2: `wiki/foo_<date>.json` is a TEMPLATE, not a filename. Treating it literally reported two
    satisfied citations as broken."""
    assert resolves("wiki/fba_within_gene_ranking_<date>.json"), \
        "a dated template should resolve when any file shares its literal prefix"
    assert not resolves("wiki/definitely_not_a_real_prefix_<date>.json")


def test_a_citation_group_is_satisfied_when_ANY_expansion_exists():
    """BUG 3: flattening the groups made a satisfied brace citation report its non-existent siblings.
    `wiki/x_{a,b}.md` means the work is findable, not that both files exist."""
    bullet = ("- **X** see `wiki/fba_within_gene_ranking_{2026-08-29,2099-01-01}.json` for the thing.")
    r = analyse(bullet)
    assert r["n_cited"] == 1, "one citation as written, not one per expansion"
    assert r["missing"] == [], "the 2026 file exists, so the group is satisfied"


def test_a_bullet_with_no_external_store_is_never_a_compression_candidate():
    """The safety property: an uncited bullet is the only copy of what it knows, however long."""
    long_uncited = "- **Y** " + ("word " * 400)
    assert analyse(long_uncited)["compressible"] is False


def test_a_long_bullet_with_a_real_store_is_a_candidate():
    long_cited = ("- **Z** " + ("word " * 400)
                  + " see `wiki/fba_within_gene_ranking_2026-08-29.md`")
    assert analyse(long_cited)["compressible"] is True
