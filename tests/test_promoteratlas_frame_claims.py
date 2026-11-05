"""Guards that every surface restating the frame-conditioning claim agrees with the measurement.

The pattern these enforce is the one that has bitten this repo repeatedly: a number is measured
correctly, then a prose surface restates it and drifts, and nothing notices because no test looks at
prose. So these check the RELATIONSHIPS -- that the memo's headline figures are the artifact's, that
the bounded verdict is never flattened into "position excluded", that the arm's metric survived the
additive write, and that the retired candidate carries the marker its own id function derives.
"""

from __future__ import annotations

import json
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
MEMO = REPO / "wiki/promoteratlas_frame_result_2026-11-01.md"
RESULT = REPO / "wiki/promoteratlas_frame_conditioning_2026-11-01.json"
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
EXPOSURE_MEMO = REPO / "wiki/promoteratlas_exposure_result_2026-11-01.md"
LEDGER = REPO / "project_state/glm-learned-representation-2026-10-07.md"

pytestmark = pytest.mark.skipif(
    not (MEMO.is_file() and RESULT.is_file()), reason="frame result not recorded yet"
)


def _flat(text: str) -> str:
    """Collapse whitespace so a phrase assertion does not depend on where a line happens to wrap.

    Earned: `as the authors ran it` wraps between "the" and "authors" in the memo, so a plain
    substring check failed on a sentence that was present and correct. Prose assertions must be
    wrap-insensitive or they test the formatter, not the claim.
    """
    return " ".join(text.split())


@pytest.fixture(scope="module")
def memo() -> str:
    return _flat(MEMO.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def res() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


# ------------------------------------------------------- the memo's numbers are the artifact's


def test_the_memos_headline_figures_are_the_artifacts(memo, res):
    c1 = res["c1_does_the_encoder_know_corpus_membership"]
    c2 = res["c2_matched_strength_series"]
    c3 = res["c3_stratified"]
    for value in (
        c1["auroc_embedding_median"], c1["auroc_gc_only_median"],
        c2["gc_median"], c2["gc_plus_frame_median"], c2["frame_alone_median"],
        c2["margin_over_gc"], c2["margin_over_gc_plus_frame"], c2["frame_absorption"],
        c3["margin_in_window"], c3["margin_out_window"],
        c3["random_stratum_min_margin_median"],
    ):
        assert "%.4f" % value in memo, "%.4f is in the artifact but not in the memo" % value


def test_the_memo_states_the_verdict_the_grader_computed(memo, res):
    assert res["verdict"] in memo


def test_the_memo_reports_the_control_premise_figures(memo, res):
    mc = res["c2_matched_strength_series"]["matched_control_premise"]
    for value in (mc["dinucleotide_alone_median"], mc["kmer4_alone_median"],
                  mc["frame_alone_median"]):
        assert "%.4f" % value in memo


def test_the_memo_reports_the_seed_spread_rather_than_only_the_median(memo, res):
    st = res["c3_stratified"]["stability_margin_out_window"]
    assert "%.4f" % st["min"] in memo or ("%.3f" % st["min"]) in memo
    assert "%.4f" % st["max"] in memo or ("%.3f" % st["max"]) in memo
    assert str(st["n_seeds_negative"]) in memo


# ----------------------------------------------------------------------- the bound is carried


def test_the_memo_never_flattens_the_verdict_into_position_excluded(memo):
    """The pattern that turned EXPOSURE_EXCLUDED_LINEAR into 'exposure excluded' downstream."""
    low = memo.lower()
    for phrase in ("position excluded", "position explains it"):
        # the phrase may appear ONLY inside an explicit prohibition
        for idx in range(len(low)):
            i = low.find(phrase, idx)
            if i < 0:
                break
            window = low[max(0, i - 220):i]
            assert ("neither" in window or "never" in window or "not" in window), (
                "%r appears without a prohibition in front of it" % phrase
            )
            idx = i + 1


def test_the_memo_states_the_run_did_NOT_resolve_the_distinction(memo, res):
    assert res["verdict_stability"]["resolved"] is False
    low = memo.lower()
    assert "does not resolve" in low or "does not\nresolve" in low or "not resolve" in low
    assert "coin flip" in low


def test_the_memo_scopes_the_replica_to_the_rule_as_specified(memo):
    assert "as specified" in memo.lower()
    assert "as the authors ran it" in memo
    assert "200 bp" in memo and "5 bp" in memo


def test_the_memo_says_the_metric_was_not_retro_edited(memo):
    assert "keeps its measured value" in memo
    assert "0.5076" in memo


# ------------------------------------------------------- the machine-readable surface agrees


def test_the_arm_artifact_carries_the_frame_block_and_it_matches_the_result():
    from scripts.record_frame_result import build_block

    arm = json.loads(ARM.read_text(encoding="utf-8"))
    res = json.loads(RESULT.read_text(encoding="utf-8"))
    bar = json.loads((REPO / "wiki/promoteratlas_frame_acceptance_bar.json")
                     .read_text(encoding="utf-8"))
    have = arm["corpus_selection_exposure_resolved"]["frame_conditioning"]
    assert json.dumps(have, sort_keys=True) == json.dumps(build_block(res, bar), sort_keys=True), (
        "the arm artifact's frame block has drifted from the result artifact"
    )


def test_the_arms_block_still_carries_its_measured_metric():
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    assert arm["arms"]["native150"]["median"] == 0.5076
    assert arm["gc_anchor"]["median"] == 0.3112


def test_the_frame_block_does_not_claim_the_ceiling_moved():
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    blk = arm["corpus_selection_exposure_resolved"]["frame_conditioning"]
    assert "NEW ARM under R6" in blk["not_a_new_arm"]
    assert blk["resolved_by_this_run"] is False


# ------------------------------------------------------------- the parent memo was amended


def test_the_exposure_memos_original_limit_wording_is_PRESERVED():
    """Amended, never deleted -- a reader must still see what the limit said before."""
    text = EXPOSURE_MEMO.read_text(encoding="utf-8")
    assert "it remains unexcluded NONLINEARLY." in text, "the original wording was removed"
    assert "AMENDED 2026-11-05 with a POSITIONAL clause" in text
    assert "the original wording above preserved" in text


def test_the_exposure_amendment_does_not_retract_the_linear_finding():
    text = EXPOSURE_MEMO.read_text(encoding="utf-8")
    assert "different interventions" in text


# -------------------------------------------------------------------- the ledger was journaled


def test_candidate_7_is_retired_with_the_marker_its_own_id_function_derives():
    import sys

    skills = pathlib.Path.home() / ".claude/skills/soraya/scripts"
    if not (skills / "candidate_replay.py").is_file():
        pytest.skip("candidate_replay not available on this host")
    sys.path.insert(0, str(skills))
    import candidate_replay as cr

    text = LEDGER.read_text(encoding="utf-8")
    row = next(ln for ln in text.splitlines() if ln.startswith("| 7 | CDS-frame"))
    cells = [c.strip() for c in row.strip().strip("|").split("|")]
    action, status, journal = cells[1], cells[3], cells[4]
    assert "done" in status.lower(), "candidate 7 is still listed pending"
    assert cr.done_marker(action) in journal, (
        "the journal marker does not match the id derived from this row's own text"
    )


def test_one_action_log_row_was_appended_with_a_valid_class():
    """The 7-value enum is the contract for the sanctioned cheap path."""
    valid = {"propose", "research", "write-plan", "edit-local-code", "run-tests", "ask-user",
             "stop"}
    text = LEDGER.read_text(encoding="utf-8")
    row = next(ln for ln in text.splitlines() if ln.startswith("| 22 | 2026-11-05 |"))
    cells = [c.strip() for c in row.strip().strip("|").split("|")]
    assert cells[2] in valid, "action class %r is not in the enum" % cells[2]
    assert "[done:" in cells[4], "the log row carries no journal marker"


def test_the_five_pre_existing_done_rows_were_left_alone():
    """Re-doing the 2026-10-19 session's retirement work is foreclosed."""
    text = LEDGER.read_text(encoding="utf-8")
    for marker in ("[done:52c5f50e]", "[done:f8023500]", "[done:d79c6971]",
                   "[done:405ce40b]", "[done:f960fd53]"):
        assert marker in text, "%s went missing" % marker
