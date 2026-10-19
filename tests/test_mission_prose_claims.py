"""Prose-claim guards for `project_state/trait-to-edit/big-idea.md`.

Three kinds of guard:

  1. RELATIONSHIP guards, not substring guards. Asserting `VALID_DIFFERENTIAL_RESPONSE`
     is present would pass even if `NOT CALIBRATED` sat in the next sentence, so the
     superseded status must be ABSENT and the new one PRESENT. Same shape as the
     existing essentiality over-claim guard, which pins a relation between adjacent
     clauses rather than the existence of a string.
  2. A CROSS-REFERENCE guard: the north-star step-4 paragraph names an arm id
     verbatim, so renaming the arm breaks this test instead of silently orphaning
     the reference.
  3. PROVENANCE for every number the correction added -- each must appear
     WHOLE-NUMBER in a named artifact, or be declared `derived` with its arithmetic.
     Whole-number matching is the point: a bare-substring check lets `0.185` match
     inside `0.1852` and a wrong number wears the costume of a cited one.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")

REPO = Path(__file__).resolve().parent.parent
MISSION = REPO / "project_state" / "trait-to-edit" / "big-idea.md"
STACK_ARM = "orthogonal-modality-stack-promcalc-plus-N-channels"
EUK_ARM = "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN"

# Numbers the correction introduced -> the artifact that carries them.
NUMBER_PROVENANCE = {
    # modality framing
    "0.4254": "wiki/glm_tile_hybrid_arm_2026-10-14.json",
    "-0.1232": "wiki/glm_tile_hybrid_arm_2026-10-14.json",
    "0.3968": "wiki/glm_tile_hybrid_arm_2026-10-14.json",
    "0.3234": "wiki/glm_tile_hybrid_arm_2026-10-14.json",
    "0.3167": "soraya_runs/arms-ledger.json",
    "0.3112": "soraya_runs/arms-ledger.json",
    "0.185": "soraya_runs/arms-ledger.json",
    "0.7928": "wiki/glm_genomewide_oracle_2026-10-07.json",
    "0.056": "wiki/forward_modality_hybrid_2026-07-17.json",
    "90.5": "wiki/forward_modality_hybrid_2026-07-17.md",
    # cabbage source concentration
    "0.3534": "wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json",
    "0.9637": "wiki/cabbage_scored_cell_campylobacter_ciprofloxacin_2026-10-10.json",
    "0.9775": "wiki/cabbage_scored_cell_salmonella_ceftriaxone_2026-10-11.json",
    "0.9597": "wiki/cabbage_scored_cell_salmonella_ciprofloxacin_2026-10-11.json",
}

# Numbers that are a RENDERING of an artifact value rather than a literal in it.
# Each carries the arithmetic, which a test below recomputes.
DERIVED_NUMBERS = {
    "0.0441": ("wiki/f3_deconfound_calibration_1001g_2026-10-10.json",
               "round(positive_control_result.median_drop, 4)"),
    "58.4": ("wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json",
             "largest_country_share 0.5844 rendered as a percentage"),
}


def _text() -> str:
    return MISSION.read_text(encoding="utf-8")


def _whole_number_present(value: str, haystack: str) -> bool:
    """`value` as a complete number, not as a substring of a longer one."""
    return re.search(rf"(?<![\d.]){re.escape(value)}(?![\d])", haystack) is not None


SUPERSESSION_MARKERS = (
    "superseded", "supersedes", "read as", "no longer", "was wrong",
    "is not true", "retract",
)


def _live_claim_occurrences(phrase: str, text: str, window: int = 400) -> list[str]:
    """Occurrences of `phrase` that are NOT inside a supersession context.

    A plain absence check cannot distinguish a stale claim from a RETRACTION of it,
    and the retraction is worth keeping -- a reader who remembers the old status
    needs to see it named as superseded. So the guard is relational: every
    occurrence must sit near a marker that denies it.
    """
    out = []
    for match in re.finditer(re.escape(phrase), text, flags=re.IGNORECASE):
        start = max(0, match.start() - window)
        context = text[start:match.end() + window].lower()
        if not any(marker in context for marker in SUPERSESSION_MARKERS):
            out.append(text[start:match.end() + 120])
    return out


def _squash(s: str) -> str:
    """Collapse whitespace so a quote that wraps across lines still matches."""
    return " ".join(s.split())


# --------------------------------------------------------------------------- #
# 1. F3 -- relationship guard
# --------------------------------------------------------------------------- #

def test_the_superseded_f3_status_survives_only_as_a_retraction():
    text = _text()
    stale = _live_claim_occurrences("NOT CALIBRATED", text)
    assert not stale, f"NOT CALIBRATED appears as a LIVE claim: {stale}"
    stale = _live_claim_occurrences("positive_control_not_fold_compatible", text)
    assert not stale, (
        f"that blocked: verdict is an outcome the run did NOT return: {stale}"
    )
    assert "VALID_DIFFERENTIAL_RESPONSE" in text


def test_the_supersession_guard_is_not_vacuous():
    """It must actually fire on the pre-correction wording."""
    pre_fix = (
        "**F3's substrate is resolved.** The de-confounding instrument is built and "
        "**NOT CALIBRATED** (`blocked:positive_control_not_fold_compatible`), but "
        "1001-Genomes is 330 MB, MD5-verified."
    )
    assert _live_claim_occurrences("NOT CALIBRATED", pre_fix), (
        "the guard would not have caught the status it exists to catch"
    )
    assert _live_claim_occurrences("deliberately absent", "Step 4 is deliberately absent.")


def test_f3_correction_carries_the_artifacts_own_limit_verbatim():
    """A calibration result without its scope invites exactly the over-read it forbids."""
    import json

    text = _text()
    art = json.loads(
        (REPO / "wiki/f3_deconfound_calibration_1001g_2026-10-10.json")
        .read_text(encoding="utf-8")
    )
    limit = next(l for l in art["honest_limits"] if "bounds the INSTRUMENT" in l)
    # The quoted fragment must be the artifact's own wording, not a paraphrase.
    # Whitespace is squashed because the quote wraps across lines in the document.
    core = "It does not validate any downstream organism-level prediction."
    assert core in _squash(limit)
    assert core in _squash(text)
    assert "ONE positive and ONE negative panel" in text


def test_f3_correction_names_both_arms_of_the_differential():
    """A differential claim needs BOTH numbers; one alone is not a differential."""
    text = _text()
    assert _whole_number_present("0.0441", text)
    assert "0.0000" in text
    assert "bi-parental" in text


def test_f3_points_at_the_live_status_rather_than_asserting_permanence():
    assert "DECONFOUNDER_CALIBRATION" in _text()


# --------------------------------------------------------------------------- #
# 2. north-star step 4 -- cross-reference guard
# --------------------------------------------------------------------------- #

def test_step_four_is_carried_by_a_named_arm_not_declared_absent():
    text = _text()
    stale = _live_claim_occurrences("deliberately absent", text)
    assert not stale, f"'deliberately absent' still reads as a live claim: {stale}"
    assert EUK_ARM in text, "the carrier arm must be named verbatim"
    assert "CARRIED by a declared, unrun arm" in text


def test_the_named_carrier_arm_actually_exists_in_the_ledger():
    """So the cross-reference cannot dangle."""
    import json

    ledger = json.loads(
        (REPO / "soraya_runs" / "arms-ledger.json").read_text(encoding="utf-8")
    )
    names = {n for arms in ledger["criteria"].values() for n in arms}
    assert EUK_ARM in names
    assert STACK_ARM in names


def test_step_four_rewording_does_not_reopen_the_embedding_negative():
    text = _text()
    assert "closed negative on de-confounded EMBEDDINGS" in text
    assert "CONSTRUCTED eukaryotic panel" in text


def test_sigma_budget_is_unchanged_by_the_rewording():
    """Answering substrate-is-not-the-subject must not cost a duplicate gap."""
    import big_idea

    bi = big_idea.parse_big_idea(_text())
    assert sum(g.strategy_budget for g in bi.gaps) == 12
    assert len(bi.gaps) == 4
    assert "Σ stays 12" in _text() or "stays 12" in _text()


# --------------------------------------------------------------------------- #
# 3. the modality correction
# --------------------------------------------------------------------------- #

def test_modality_correction_keeps_the_generator_bottleneck_claim():
    text = _text()
    assert "GENERATOR, not the scorer" in text
    assert "implied LEVER was wrong" in text


def test_modality_correction_states_the_shuffled_control_hurts():
    """Without the control the 0.4254 win says nothing about channel content."""
    text = _text()
    assert _whole_number_present("-0.1232", text)
    assert "channel CONTENT" in text


def test_modality_correction_does_not_overturn_the_encoder_gate_negative():
    text = _text()
    assert "does NOT overturn the G-A encoder-gate negative" in text
    assert "Do not quote +0.0122 as a win." in text, (
        "the pre-existing encoder-gate caveat must survive the edit"
    )


def test_the_non_comparable_arm_is_flagged_as_non_comparable():
    """0.185 is against a different ceiling; printing it beside 0.4254 invites a false rank."""
    text = _text()
    assert _whole_number_present("0.185", text)
    assert _whole_number_present("0.7928", text)
    assert "not\ncomparable" in text or "not comparable" in text


# --------------------------------------------------------------------------- #
# 4. the cabbage diversity caveat
# --------------------------------------------------------------------------- #

def test_diversity_caveat_reports_three_of_four_concentrated():
    text = _text()
    assert "3 of 4" in text
    assert "source_concentrated" in text and "source_diverse" in text


def test_the_usa_claim_is_not_overstated_as_all_four():
    """The plan's own draft said 'all four are 96-98% USA'. It is not true.

    Three cells are 96-98% USA; the named source_diverse cell is 58.4%. Asserting
    both is what stops the inaccurate blanket version creeping back.
    """
    text = _text()
    assert "96-98% USA" in text
    assert _whole_number_present("58.4", text)
    assert "all four are 96-98%" not in text


def test_diversity_caveat_explains_why_the_remaining_predicate_is_open():
    text = _text()
    assert "cabbage_bioproject_diversity.json" in text
    assert "single-country majority" in text


# --------------------------------------------------------------------------- #
# 5. number provenance
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("value,artifact", sorted(NUMBER_PROVENANCE.items()))
def test_every_added_number_appears_whole_in_its_cited_artifact(value, artifact):
    path = REPO / artifact
    assert path.exists(), f"{artifact} does not exist"
    assert _whole_number_present(value, _text()), f"{value} is not in the mission text"
    haystack = path.read_text(encoding="utf-8", errors="replace")
    assert _whole_number_present(value, haystack), (
        f"{value} is cited to {artifact} but does not appear there as a whole number"
    )


def test_the_derived_numbers_recompute_from_their_artifacts():
    """A `derived` provenance kind must carry arithmetic that actually checks out."""
    import json

    f3 = json.loads(
        (REPO / DERIVED_NUMBERS["0.0441"][0]).read_text(encoding="utf-8")
    )
    assert round(f3["positive_control_result"]["median_drop"], 4) == 0.0441

    cell = json.loads(
        (REPO / DERIVED_NUMBERS["58.4"][0]).read_text(encoding="utf-8")
    )
    share = cell["source_concentration_of_scored_set"]["largest_country_share"]
    assert round(share * 100, 1) == 58.4


def test_whole_number_matching_is_actually_strict():
    """Non-vacuity for the matcher itself."""
    assert _whole_number_present("0.185", "metric 0.185 here")
    assert not _whole_number_present("0.185", "metric 0.1852 here")
    assert not _whole_number_present("0.185", "metric 10.185 here")
    assert _whole_number_present("-0.1232", "delta -0.1232 measured")
