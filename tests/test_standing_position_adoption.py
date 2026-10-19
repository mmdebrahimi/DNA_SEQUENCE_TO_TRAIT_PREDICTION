"""Guards for the ADOPTED standing position.

The dangerous failure here is not a wrong status -- it is a FALSE PROVENANCE: closing
a position about GENERATION by pointing at a DISCRIMINATIVE result. So the strongest
test asserts `tested_exactly` does NOT present a wiki measurement as the resolver.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")

import standing_positions as sp  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
REGISTER = REPO / "soraya_runs" / "standing-positions.json"
MEMO = REPO / "wiki" / "standing_position_adoption_2026-10-19.md"
PID = "learned-model-is-the-product"
OTHER = "substrate-is-not-the-subject"

# The user's three restatements, pinned verbatim. Adoption must not rewrite them.
QUOTES = [
    "just build it... I don't care if it is not getting good results",
    "I don't want you to compare it with the lookup table, that was just built as a bridge",
    "I had qwen do some research on the topic, and a talk on genomics by a professor "
    "that seems very related. please analyze and tell me what you think",
]


def _register() -> dict:
    return json.loads(REGISTER.read_text(encoding="utf-8"))


def _position() -> dict:
    return _register()["positions"][PID]


def test_the_position_is_adopted_with_the_module_constant():
    pos = _position()
    assert pos["status"] == sp.ADOPTED == "adopted-user-was-right"
    assert pos["resolved_at"]


def test_the_register_still_round_trips_through_its_owning_dataclass():
    """An out-of-schema field makes `--show` crash: Position(**row) rejects extras.

    This happened while writing this step -- two explanatory fields were added and
    the owning CLI stopped rendering. The rationale lives in a memo instead.
    """
    for row in _register()["positions"].values():
        sp.Position(**row)  # raises TypeError on any stray key
    assert sp.PositionRegister(cwd=str(REPO)).render()


def test_the_loud_three_times_open_flag_no_longer_fires():
    rendered = sp.PositionRegister(cwd=str(REPO)).render()
    assert "[ADOPTED — user was right] learned-model-is-the-product" in rendered
    assert f"[OPEN] {PID}" not in rendered


def test_tested_exactly_states_that_nothing_was_measured():
    tested = _position()["tested_exactly"]
    assert "NO experiment tested this claim" in tested
    assert "ADOPTION IS NOT EVIDENCE" in tested
    assert "authority" in tested.lower()


def test_tested_exactly_does_not_present_a_wiki_measurement_as_the_resolver():
    """The proxy trap, guarded at string level.

    A `wiki/...` artifact named in `tested_exactly` would read as 'this is what was
    measured'. The one relevant artifact is DISCRIMINATIVE while the position is
    about GENERATION, so it is referenced only as further reading, and the single
    mention must be the memo pointer -- never a result artifact.
    """
    tested = _position()["tested_exactly"]
    cited = [tok for tok in tested.replace(";", " ").split()
             if tok.startswith("wiki/")]
    assert cited == ["wiki/standing_position_adoption_2026-10-19.md."], cited
    for forbidden in ("glm_headroom_allrt", "glm_alphabet_headroom", "0.8141", "0.8142"):
        assert forbidden not in tested, (
            f"{forbidden} must not appear in tested_exactly -- it measures a "
            "DISCRIMINATIVE quantity, not the generative claim"
        )


def test_proxy_warning_is_silent_on_the_adoption_text():
    pos = _position()
    assert sp.proxy_warning(pos["claim"], pos["tested_exactly"]) == ""


def test_proxy_warning_would_have_fired_on_an_empty_field():
    """Non-vacuity for the guard above."""
    assert sp.proxy_warning("x", "") != ""
    assert sp.proxy_warning("x", "used a k-mer baseline for the learned model") != ""


def test_resolved_by_names_a_file_that_exists_and_contains_the_claim():
    pos = _position()
    target = REPO / pos["resolved_by"]
    assert target.exists(), pos["resolved_by"]
    text = target.read_text(encoding="utf-8")
    # The exact premise the adoption rests on must be present in that file.
    assert "The product is the generative engine" in text
    assert "deterministic catalog was a BRIDGE" in text


def test_the_users_own_words_are_byte_unchanged():
    assert _position()["quotes"] == QUOTES


def test_times_raised_is_unchanged_by_adoption():
    assert _position()["times_raised"] == 3


def test_the_other_position_stays_open():
    other = _register()["positions"][OTHER]
    assert other["status"] == sp.OPEN
    assert other["times_raised"] == 2
    assert other["resolved_by"] == ""
    assert other["tested_exactly"] == ""


def test_the_memo_records_the_evidence_that_was_deliberately_not_cited():
    text = MEMO.read_text(encoding="utf-8")
    assert "ADOPTION IS NOT EVIDENCE" in text
    assert "DISCRIMINATIVE" in text and "GENERATION" in text
    assert "orientation-error 0" in text
    assert "proxy_warning()" in text


def test_the_memo_pairs_each_number_with_its_own_artifact():
    """0.8141 and 0.8142 live in DIFFERENT artifacts and are one digit apart."""
    text = MEMO.read_text(encoding="utf-8")
    assert "median 0.8141" in text
    assert "glm_headroom_allrt_2026-09-30.md" in text
    assert "glm_alphabet_headroom_2026-09-30.md" in text
    for rel, value in (
        ("wiki/glm_headroom_allrt_2026-09-30.md", "0.8141"),
        ("wiki/glm_alphabet_headroom_2026-09-30.md", "0.8142"),
    ):
        assert value in (REPO / rel).read_text(encoding="utf-8", errors="replace"), (
            f"{value} is attributed to {rel} but does not appear there"
        )


def test_the_memo_states_the_modality_consistency():
    text = MEMO.read_text(encoding="utf-8")
    assert "contains a learned channel" in text.lower()
    assert "-0.1232" in text


def test_the_memo_records_what_stays_open():
    text = MEMO.read_text(encoding="utf-8")
    assert OTHER in text
    assert "does not resolve it" in text
