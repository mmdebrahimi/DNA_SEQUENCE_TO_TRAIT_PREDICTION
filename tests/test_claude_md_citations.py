"""Every repo path CLAUDE.md cites must exist — or say plainly that it does not.

WHY THIS EXISTS
`CLAUDE.md` is loaded into EVERY session, so a stale line there misdirects every future run — not once,
repeatedly. That is not hypothetical: it happened TWICE in the 2026-08-25 session alone.

  * The candidate table said the TB independent gold-set was still to be fetched. It had been DELIVERED
    two months earlier via a cohort 48x larger, so an `--advance` run nearly spent hours re-fetching a
    number that already existed.
  * CLAUDE.md said the ProSST "real forward pass is deferred to a Kaggle run". It had RUN the same day it
    shipped, locally on CPU, with an exact column reproduction and a powered lift
    (`wiki/prosst_lift_2026-07-18.md`) -- and that stale line was repeated as a recommendation before
    anyone opened the artifact.

Neither was catchable by a test at the time, because nothing checked CLAUDE.md against the filesystem.
This closes the mechanically-checkable half: a cited path either resolves, or the text next to it admits
it does not.

WHAT THIS CANNOT DO (stated so the guard is not over-trusted)
It checks that cited FILES exist. It cannot check that a cited file still says what CLAUDE.md claims it
says -- the ProSST case would NOT have been caught here, because `wiki/prosst_lift_2026-07-18.md` existed
the whole time; what was stale was the CLAIM ABOUT it. Reading the artifact before repeating a claim
remains model discipline. `test_report_card_doc_sync.py` is the one place a specific NUMBER is pinned.

SCOPE NOTE (a real false-positive class, learned while building this)
Only paths containing a directory separator are checked. CLAUDE.md refers to many files by bare basename
in prose (`cache.py`, `cli.py`, `amr_rules.py`), and a naive regex flags 72 "missing" files of which
almost all are that. Requiring a `/` takes it from 72 noise hits to 3 real ones.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

CLAUDE_MD = ROOT / "CLAUDE.md"

# A cited path may be absent ONLY when the surrounding text says so. Each entry records the reason, and
# the test VERIFIES the disclaimer is still present -- an allowlist that stops explaining itself is just
# a mute button.
_DISCLOSED_ABSENT = {
    "reports/cipro_v0_scope_limit_decision_2026-05-23.md": "not yet on origin",
    ".claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json": "no longer exists on this host",
}

_PATH = re.compile(r"`([A-Za-z0-9_./-]+/[A-Za-z0-9_./-]+\.(?:py|md|json|yaml|yml|tsv|fna|csv|lock|toml|bed))`")


def _cited() -> set[str]:
    return set(_PATH.findall(CLAUDE_MD.read_text(encoding="utf-8", errors="replace")))


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_every_repo_path_cited_in_claude_md_exists():
    cited = _cited()
    missing = sorted(p for p in cited if not (ROOT / p).exists())
    undisclosed = [p for p in missing if p not in _DISCLOSED_ABSENT]
    assert not undisclosed, (
        f"CLAUDE.md cites repo paths that do not exist: {undisclosed}. CLAUDE.md is loaded into every "
        f"session, so a stale citation misdirects every future run. Fix the path, or state plainly that "
        f"the file is absent and add it to _DISCLOSED_ABSENT with the reason.")


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_each_disclosed_absent_path_still_carries_its_disclaimer():
    """An allowlisted absence must keep SAYING it is absent, right where it is cited."""
    text = CLAUDE_MD.read_text(encoding="utf-8", errors="replace")
    for path, phrase in _DISCLOSED_ABSENT.items():
        if (ROOT / path).exists():
            continue                      # it came back -- nothing to disclose, and the guard above passes
        i = text.find(path)
        assert i != -1, f"{path} is allowlisted but no longer cited -- drop it from _DISCLOSED_ABSENT"
        window = text[i:i + 400]
        assert phrase in window, (
            f"CLAUDE.md cites the ABSENT {path} without saying so nearby (expected {phrase!r}). "
            f"An absent file must be described as absent, not silently cited.")


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_the_citation_scan_is_not_vacuous():
    """If the regex ever stopped matching, every assertion above would pass while checking nothing."""
    cited = _cited()
    assert len(cited) >= 100, f"only {len(cited)} paths matched -- the citation regex likely broke"
    # and it must be finding real, resolvable files, not just strings
    assert sum(1 for p in cited if (ROOT / p).exists()) >= 100


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_the_renamed_guard_file_is_not_still_cited_under_its_old_name():
    """REGRESSION: `tests/test_decode_router_commands.py` was renamed to `test_advertised_commands.py`
    on 2026-08-23 when its scope outgrew the router, and one CLAUDE.md citation kept the old name."""
    text = CLAUDE_MD.read_text(encoding="utf-8", errors="replace")
    assert "test_decode_router_commands" not in text
    assert (ROOT / "tests" / "test_advertised_commands.py").exists()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))


# --------------------------------------------------------------------------------------------------
# STALE-DEFERRAL regressions (2026-08-25).
#
# The guard above checks a cited FILE exists. It cannot check that CLAUDE.md's CLAIM about the file is
# still true -- and that is the half which has now bitten four times: the TB gold-set candidate row, the
# ProSST "deferred to a Kaggle run" line, and the three pinned below.
#
# WHY THERE IS NO GENERAL SCREEN FOR THIS (measured, not assumed)
# The obvious rule -- "a deferral marker within N chars of an existing repo path is a contradiction" --
# was built and run: it returned 5 hits, ALL false positives. The marker and the path sit in the same
# prose region but refer to different things (`tests/test_genome_map_browser.py` near "Still deferred:
# pathway/KEGG"), and three of the five were inside the correction text explaining a fix. A guard with a
# ~100% false-positive rate gets disabled, so it was NOT shipped. Natural-language claims about artifacts
# resist mechanical checking; reading the artifact before repeating a claim stays model discipline.
#
# What IS worth pinning is each specific correction, so a stale claim cannot silently return.
# --------------------------------------------------------------------------------------------------

def _md() -> str:
    return CLAUDE_MD.read_text(encoding="utf-8", errors="replace")


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_the_genome_map_browser_is_not_described_as_deferred():
    """CLAUDE.md CONTRADICTED ITSELF: one line said "a visual browser deferred" while the next bullet said
    "GRAPHICAL BROWSER SHIPPED 2026-07-11". Both files exist."""
    text = _md()
    assert (ROOT / "dna_decode" / "genome_map" / "browser.py").exists()
    assert (ROOT / "scripts" / "genome_map_browser.py").exists()
    assert "a visual browser deferred" not in text
    assert "GRAPHICAL BROWSER SHIPPED" in text


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_the_bvbrc_census_is_not_described_as_a_deferred_3_drug_run():
    """It RAN 2026-05-18, as a FOUR-drug census -- the old line was stale on both count and status."""
    assert (ROOT / "wiki" / "bvbrc_strict_mic_4drug_census_2026-05-18.md").exists()
    text = _md()
    assert "strict-MIC 3-drug feasibility census" not in text
    assert "RAN 2026-05-18 as a FOUR-drug census" in text


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_the_tb_pending_data_runs_header_records_that_both_blockers_resolved():
    """The header still read "PENDING ... BLOCKED-gated by design" while its OWN sub-bullets already
    recorded the resolution: the parquet adapter sidestepped the regeno fetch, the AMR-Portal cohort
    delivered the independent number, and the callability probe measured the correction."""
    text = _md()
    assert (ROOT / "wiki" / "tb_independent_amr_portal_scores.json").exists()
    assert (ROOT / "wiki" / "tb_callability_probe_2026-07-10.json").exists()
    assert "PENDING DATA RUNS (BLOCKED-gated by design, not incomplete code)" not in text
    assert "BOTH BLOCKERS SINCE RESOLVED" in text


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_the_prosst_forward_pass_is_not_described_as_deferred():
    """The 2026-08-25 case: the line said "deferred to a Kaggle run"; it had run locally on CPU the same
    day it shipped, and that stale claim was repeated as a recommendation."""
    assert (ROOT / "wiki" / "prosst_lift_2026-07-18.md").exists()
    text = _md()
    assert "the real forward pass is deferred to a Kaggle run" not in text
    assert "RAN the same day, LOCALLY on CPU" in text
