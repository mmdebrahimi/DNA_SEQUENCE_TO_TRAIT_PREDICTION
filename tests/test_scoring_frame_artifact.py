"""Guards for the Step 5 resolution artifact and the ledger annotation.

Three of these exist because of specific recorded incidents.

An invalid hand-edited JSON artifact was committed and cited once, so the artifact is re-parsed.
A curve trend quoted without its null band would be the retracted statistic's mistake in a new
costume, so no trend may appear without its band in the same block. And the ledger is read by
nine scripts that treat `len(criteria[<crit>])` as the arm count, so the annotation must leave
that block byte-identical.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ART_JSON = Path("wiki/glm_scoring_frame_resolution.json")
ART_MD = Path("wiki/glm_scoring_frame_resolution.md")
LEDGER = Path("soraya_runs/arms-ledger.json")

pytestmark = pytest.mark.skipif(
    not ART_JSON.exists(), reason="resolution artifact not written yet")


@pytest.fixture(scope="module")
def art():
    return json.loads(ART_JSON.read_text(encoding="utf-8"))


# ----------------------------------------------------------------------------------------
# shape and provenance
# ----------------------------------------------------------------------------------------

def test_the_artifact_parses_as_json():
    """A hand-edited JSON artifact was committed invalid and cited once."""
    json.loads(ART_JSON.read_text(encoding="utf-8"))


def test_the_markdown_companion_exists():
    assert ART_MD.exists()
    assert ART_MD.read_text(encoding="utf-8").strip()


def test_required_keys_are_present(art):
    for k in ("schema", "date", "question", "headline", "supersedes",
              "never_cite_for_provenance", "honest_limits", "cites",
              "finding_1_reproducibility_rises_across_the_range",
              "finding_2_the_decisive_gate_did_not_grade",
              "finding_3_the_margin_grows_on_the_active_frame",
              "finding_4_the_human_label_is_not_similarly_concentrated",
              "the_published_eukaryotic_verdict_is_UNTOUCHED_and_why",
              "which_of_the_four_eukaryotic_differences_this_resolves"):
        assert k in art, f"missing required key {k}"


def test_every_cited_wiki_path_exists(art):
    missing = [p for p in art["cites"] if not Path(p).exists()]
    assert not missing, f"cited but absent: {missing}"


def test_the_supersedes_block_names_the_retracted_quantity(art):
    s = art["supersedes"]
    assert "0.8194" in s["retracted_statistic"]
    assert "closed form" in s["why_vacuous"]
    assert "reproducibility_curve" in s["replacement"]
    assert "NOT deleted" in " ".join(s.keys()) or s["the_retracted_value_is_NOT_deleted"]


# ----------------------------------------------------------------------------------------
# no trend without its null band
# ----------------------------------------------------------------------------------------

def test_no_curve_trend_appears_without_its_null_band_in_the_same_block(art):
    """Quoting a trend bare is the retracted statistic's mistake in a new costume."""
    f1 = art["finding_1_reproducibility_rises_across_the_range"]
    assert "trend" in f1
    assert "attenuation_null" in f1, "a trend with no null band in the same block"
    for key in ("mean", "p95", "max", "min"):
        assert key in f1["attenuation_null"], f"null band lacks {key}"
    assert "beyond_null_max" in f1
    mirror = f1["mirror_replication"]
    assert "trend" in mirror and "null_max" in mirror, "the mirror trend has no band either"


def test_the_trend_is_actually_beyond_the_null_max_as_claimed(art):
    f1 = art["finding_1_reproducibility_rises_across_the_range"]
    assert f1["trend"] > f1["attenuation_null"]["max"]
    assert f1["beyond_null_max"] is True
    assert f1["mirror_replication"]["trend"] > f1["mirror_replication"]["null_max"]


def test_the_reads_as_field_insists_on_reading_the_trend_against_the_null(art):
    f1 = art["finding_1_reproducibility_rises_across_the_range"]
    assert "null" in f1["reads_as"].lower()


# ----------------------------------------------------------------------------------------
# the ungraded gate must not be reported as a result
# ----------------------------------------------------------------------------------------

def test_the_ungraded_gate_numbers_are_flagged_as_not_a_result(art):
    f2 = art["finding_2_the_decisive_gate_did_not_grade"]
    assert f2["verdict"].startswith("INDETERMINATE")
    assert "ungraded_numbers_DO_NOT_QUOTE_AS_A_RESULT" in f2
    block = f2["ungraded_numbers_DO_NOT_QUOTE_AS_A_RESULT"]
    assert "why_it_is_withheld" in block
    assert "what_the_bar_would_have_returned" in block


def test_the_null_diagnosis_owns_the_bar_specification_error(art):
    d = art["finding_2_the_decisive_gate_did_not_grade"]["null_diagnosis"]
    assert "INHERITED" in d["my_error_in_the_bar"]
    assert d["changing_the_gate_is_a_user_authority_call"] is True
    assert "what_is_NOT_claimed" in d, "the diagnosis must bound what it does not show"


def test_the_headline_says_the_gate_did_not_grade(art):
    assert "INDETERMINATE" in art["headline"]


# ----------------------------------------------------------------------------------------
# no bare ratio without the frame caveat
# ----------------------------------------------------------------------------------------

def test_the_artifact_never_states_the_ratio_without_a_frame_caveat():
    """A bare 20-25x anywhere in the memo would undo the point of the memo."""
    md = ART_MD.read_text(encoding="utf-8")
    for para in md.split("\n\n"):
        if "20-25" in para or "20–25" in para:
            low = para.lower()
            assert any(w in low for w in ("frame", "confound", "caveat", "unresolved",
                                          "partially")), (
                f"a paragraph states the ratio with no frame caveat:\n{para[:200]}"
            )


def test_the_published_verdict_is_recorded_as_untouched(art):
    v = art["the_published_eukaryotic_verdict_is_UNTOUCHED_and_why"]
    assert v["so_no_reinterpretation_is_needed"] is True
    assert "CI" in v["verdict_detail_is_a_CI_COMPARISON"]
    assert art["frozen_surface_touched"] is False


# ----------------------------------------------------------------------------------------
# the ledger annotation is augment-only
# ----------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ledger():
    return json.loads(LEDGER.read_text(encoding="utf-8"))


def test_the_ledger_still_holds_exactly_sixteen_arms(ledger):
    crit = next(k for k in ledger["criteria"] if k.startswith("ceiling:"))
    assert len(ledger["criteria"][crit]) == 16, (
        "len(criteria[crit]) IS the arm count for scripts/declare_successor_arms.py; a "
        "criterion-level annotation would read as a 17th arm"
    )


def test_every_criterion_value_is_still_arm_shaped(ledger):
    """The collision this annotation was moved out of the criterion to avoid."""
    crit = next(k for k in ledger["criteria"] if k.startswith("ceiling:"))
    for name, arm in ledger["criteria"][crit].items():
        assert isinstance(arm, dict), f"{name} is not a dict; consumers call arm.get('status')"
        assert "status" in arm, f"{name} has no status; it would read as a malformed arm"


def test_the_ledger_carries_NO_disclosure_and_the_artifact_says_so():
    """The annotation was REVERTED, and this pins the reverted state plus its explanation.

    These tests previously asserted a top-level `frame_disclosures` key. It was implemented and
    then reverted: there is no augment-only location in this file. A criterion-level key reads as
    a 17th arm (`len(criteria[crit])` IS the arm count,
    `scripts/declare_successor_arms.py:172-173`) and a top-level key breaks a pinned key set
    (`tests/test_record_successor_arm_results.py:367` asserts
    `set(data) == {"schema_version", "criteria"}`). Pinning the ABSENCE keeps the reason on the
    record instead of letting a future session re-attempt it and re-discover the same wall.
    """
    ledger = json.loads(LEDGER.read_text(encoding="utf-8"))
    assert set(ledger.keys()) == {"schema_version", "criteria"}, (
        "a new top-level key appeared; test_record_successor_arm_results.py:367 pins this set "
        "exactly, so anything added here fails there instead"
    )
    assert "frame_disclosures" not in ledger


def test_the_artifact_does_not_CLAIM_a_ledger_annotation_it_did_not_make():
    """The failure this exists for actually happened.

    After the revert, the memo still read "`soraya_runs/arms-ledger.json` gained one top-level
    `frame_disclosures` key" -- a false statement about a frozen-by-convention file, on the
    artifact that is this work's provenance. The claim and the file disagreed for the length of
    one execution run.
    """
    md = ART_MD.read_text(encoding="utf-8")
    assert "gained one **top-level**" not in md
    section = md[md.index("## Ledger"):]
    assert "byte-identical" in section
    assert "no augment-only location" in section.lower()
    # and the consequence must be stated, not just the mechanism
    assert "0.5013" in section and "0.9823" in section


def test_no_cell_registry_contract_cites_a_never_cite_path(art):
    """Pins today's clean state and fires if a future contract adds the citation."""
    from dna_decode.data.cell_registry import cells

    banned = set(art["never_cite_for_provenance"]["paths"])
    offenders = []
    for c in cells():
        blob = json.dumps(c, default=str) if not hasattr(c, "__dict__") else json.dumps(
            vars(c), default=str)
        for b in banned:
            if b in blob:
                offenders.append((getattr(c, "cell_id", "?"), b))
    assert not offenders, (
        f"a contract cites a never-cite path: {offenders}. The contract-number audit verifies "
        f"by whole-number presence, so citing that memo would certify the retracted 0.8194"
    )
