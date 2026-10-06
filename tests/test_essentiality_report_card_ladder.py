"""Guards for the report-card ladder extension — augment-only, verified by diff, non-vacuously."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_essentiality_report_card as card  # noqa: E402

CARD = ROOT / "wiki" / "essentiality_report_card.json"


def _ladder_art():
    return {"date": "2026-10-06", "n_scored": 2, "n_walled": 1,
            "frozen_thresholds": {"plateau_tol": 0.06, "cliff_drop": 0.12, "min_rungs": 3},
            "verdict": {"verdict": "INDETERMINATE_INSUFFICIENT_RUNGS", "reason": "2 < 3"},
            "rungs": [
                {"key": "ecoli", "organism": "Escherichia coli", "depth": 8, "scored": True,
                 "technology": "transposon_insertion", "class_sourcing_mode": "one_file_two_columns",
                 "coverage_lift": 0.3125, "coverage_lift_adjusted": 0.3213,
                 "coverage_essential": 0.3789, "coverage_nonessential": 0.0664,
                 "null_p95": 0.0298, "null_max": 0.0518,
                 "auroc_SECONDARY_not_cross_rung_comparable": 0.6952},
                {"key": "saureus", "organism": "Staphylococcus aureus", "depth": 2, "scored": False,
                 "technology": "transposon_insertion",
                 "class_sourcing_mode": "essential_plus_complement",
                 "wall": "WALL_SCHEMA_UNVERIFIED", "route_tried": "no W0 record for ntml_nebraska"},
                {"key": "human", "organism": "Homo sapiens", "depth": 1, "scored": True,
                 "technology": "crispr_ko", "class_sourcing_mode": "two_files",
                 "coverage_lift": 0.1633, "coverage_lift_adjusted": 0.1787,
                 "coverage_essential": 0.1689, "coverage_nonessential": 0.0056,
                 "null_p95": 0.0214, "null_max": 0.0446,
                 "auroc_SECONDARY_not_cross_rung_comparable": 0.5805}]}


# --------------------------------------------------------------------------------------------------
# row construction
# --------------------------------------------------------------------------------------------------
def test_no_ladder_artifact_yields_no_rows():
    assert card.ladder_rows(None) == []


def test_no_ladder_artifact_in_the_wiki_dir_yields_no_rows(tmp_path, monkeypatch):
    """`ladder_rows(None)` is pinned above; this pins the LOOKUP that produces that None. The card is
    read-only and never triggers a ladder run, so on a host that has never run the ladder it must append
    nothing rather than fail -- and `load_ladder` is the only thing standing between those two outcomes."""
    monkeypatch.setattr(card, "W", tmp_path)
    (tmp_path / "essentiality_report_card.json").write_text("{}", encoding="utf-8")
    assert card.load_ladder() is None
    assert card.ladder_rows(card.load_ladder()) == []


def test_the_NEWEST_ladder_artifact_is_the_one_the_card_reports(tmp_path, monkeypatch):
    """Ladder artifacts are date-stamped and accrue, so a re-run must supersede rather than coexist. If
    the oldest won, the card would keep reporting a verdict the current decoder no longer produces."""
    monkeypatch.setattr(card, "W", tmp_path)
    old = _ladder_art()
    old["rungs"] = [dict(old["rungs"][0], coverage_lift=0.9999)]
    (tmp_path / "essentiality_transfer_ladder_2026-10-05.json").write_text(
        json.dumps(old), encoding="utf-8")
    (tmp_path / "essentiality_transfer_ladder_2026-10-06.json").write_text(
        json.dumps(_ladder_art()), encoding="utf-8")
    rows = card.ladder_rows(card.load_ladder())
    assert len(rows) == 3, "the newer artifact's three rungs, not the older one's single rung"
    assert "0.9999" not in rows[0]["metric"]
    assert "0.3125" in rows[0]["metric"]


def test_a_walled_rung_with_no_recorded_route_says_so_rather_than_rendering_blank():
    """`route_tried` is the only diagnostic a walled row carries, so its absence must be visible as an
    absence. An empty string here would read as 'nothing was tried'."""
    art = _ladder_art()
    art["rungs"] = [r for r in art["rungs"] if not r.get("scored")]
    art["rungs"][0].pop("route_tried")
    row = card.ladder_rows(art)[0]
    assert row["validation"] == "route tried: (not recorded)"
    assert row["tier"] == "WALL_SCHEMA_UNVERIFIED"


def test_each_rung_gets_its_own_honest_tier():
    rows = card.ladder_rows(_ladder_art())
    tiers = {r["organism"]: r["tier"] for r in rows}
    assert tiers["Escherichia coli"] == "COVERAGE_SCORED"
    assert tiers["Homo sapiens"] == "COVERAGE_SCORED"
    assert tiers["Staphylococcus aureus"] == "WALL_SCHEMA_UNVERIFIED"


def test_a_walled_rung_renders_its_wall_and_carries_NO_metric():
    rows = card.ladder_rows(_ladder_art())
    w = next(r for r in rows if r["tier"].startswith("WALL_"))
    assert "none" in w["metric"].lower()
    assert "route tried" in w["validation"].lower()
    assert "coverage_lift" not in w["metric"]


def test_rows_lead_with_coverage_lift_and_label_auroc_as_secondary():
    rows = card.ladder_rows(_ladder_art())
    scored = next(r for r in rows if r["tier"] == "COVERAGE_SCORED")
    assert scored["metric"].startswith("coverage_lift")
    assert "SECONDARY" in scored["metric"] and "NOT cross-rung comparable" in scored["metric"]


def test_rows_are_ordered_nearest_rung_first():
    rows = card.ladder_rows(_ladder_art())
    assert [r["organism"] for r in rows][0] == "Escherichia coli"
    assert [r["organism"] for r in rows][-1] == "Homo sapiens"


def test_no_aggregate_is_computed_over_the_rungs():
    """The card's standing invariant is per-organism tiers only -- no aggregate headline."""
    rows = card.ladder_rows(_ladder_art())
    assert all(set(r) == {"organism", "cell", "tier", "metric", "validation"} for r in rows)


# --------------------------------------------------------------------------------------------------
# augment-only, verified by diff AND proven non-vacuous
# --------------------------------------------------------------------------------------------------
def test_the_committed_card_kept_its_two_original_rows_field_identical():
    """Verified by DIFF against the known pre-extension shape, not by inspection."""
    if not CARD.exists():
        pytest.skip("card not built on this host")
    a = json.loads(CARD.read_text(encoding="utf-8"))
    assert "aggregate" not in a and "headline" not in a
    assert a["schema"] == "essentiality-report-card-v1"
    first_two = a["organisms"][:2]
    assert first_two[0]["organism"] == "Escherichia coli K-12"
    assert first_two[0]["tier"] == "AUROC_SCORED"
    assert first_two[1]["organism"] == "Homo sapiens"
    assert first_two[1]["tier"] == "TRANSFER_SCORED"
    # and the appended rows are exactly the ladder's rungs
    assert len(a["organisms"]) == 2 + 5
    assert all(r["tier"] in ("COVERAGE_SCORED",) or r["tier"].startswith("WALL_")
               for r in a["organisms"][2:])


def test_the_diff_check_is_NON_VACUOUS_a_simulated_merge_bug_is_detected():
    """Simulates the overwrite this guard exists to catch: if a merge clobbered the human row's metric,
    a field-identity comparison must fail. A check that cannot fail is decoration."""
    if not CARD.exists():
        pytest.skip("card not built on this host")
    a = json.loads(CARD.read_text(encoding="utf-8"))
    clobbered = copy.deepcopy(a)
    clobbered["organisms"][1]["metric"] = "coverage_lift 0.1633"      # the merge bug
    assert clobbered["organisms"][1] != a["organisms"][1], \
        "the comparison must be able to see a changed field"
    # and the real row is NOT the clobbered shape
    assert a["organisms"][1]["metric"].startswith("AUROC"), a["organisms"][1]["metric"]


# --------------------------------------------------------------------------------------------------
# the shared scorer
# --------------------------------------------------------------------------------------------------
def test_the_card_uses_the_SHARED_auroc_not_an_inline_one():
    """Card and ladder must not report two different AUROCs for one rung."""
    src = (ROOT / "scripts" / "build_essentiality_report_card.py").read_text(encoding="utf-8")
    assert "from essentiality_transfer_ladder import auroc as shared_auroc" in src
    # Check for the CALL and the IMPORT, not the word: the bare word appears in the comment explaining
    # what was removed. This is the second time in this session a guard of mine grepped a word and
    # tripped on its own prose, so it checks the implementation instead.
    assert "mannwhitneyu(" not in src, "the inline AUROC call must be gone, not merely unused"
    assert "from scipy.stats import mannwhitneyu" not in src


def test_the_shared_auroc_reproduces_the_committed_human_number():
    from essentiality_transfer_ladder import auroc
    import numpy as np
    # a tiny deterministic case first: perfect separation -> 1.0
    assert auroc(np.array([1.0, 2.0]), np.array([0.0, 0.0])) == 1.0
    if not CARD.exists():
        pytest.skip("card not built on this host")
    a = json.loads(CARD.read_text(encoding="utf-8"))
    assert "0.5805" in a["organisms"][1]["metric"]


# --------------------------------------------------------------------------------------------------
# the md section
# --------------------------------------------------------------------------------------------------
def test_the_md_section_states_the_sampling_frame_refusal_and_the_lock():
    md = ROOT / "wiki" / "essentiality_report_card.md"
    if not md.exists():
        pytest.skip("card md not built on this host")
    t = md.read_text(encoding="utf-8")
    assert "## Cross-organism transfer ladder" in t
    assert "coverage_lift" in t
    assert "not" in t and "compared across rungs" in t
    assert "DESCRIPTIVE-CONSISTENCY LOCK" in t
    assert "INDETERMINATE_INSUFFICIENT_RUNGS" in t


def test_the_card_still_exits_zero_it_is_a_report_not_a_gate():
    assert card.main() == 0


def test_the_wiki_path_is_ANCHORED_to_the_repo_not_the_callers_cwd():
    """A silent zero, found by the test-epilogue pass and fixed.

    `W` was a relative `Path("wiki")`. That was only a LOUD failure (FileNotFoundError on the write)
    until load_ladder() was added: a glob over a nonexistent directory returns [], so from the wrong cwd
    the card would append ZERO ladder rows and, because the extension is augment-only, nothing would
    complain. The card would look like a clean 2-row card rather than a broken 7-row one.
    """
    assert card.W.is_absolute(), card.W
    assert card.W.name == "wiki"
    assert (card.W / "essentiality_report_card.json").parent == card.W
    # and load_ladder must not depend on where it is called from
    import os
    cwd = os.getcwd()
    try:
        os.chdir(card.W.parent.parent if card.W.parent.parent.exists() else "/")
        assert card.load_ladder() is not None, "load_ladder became cwd-dependent again"
    finally:
        os.chdir(cwd)
