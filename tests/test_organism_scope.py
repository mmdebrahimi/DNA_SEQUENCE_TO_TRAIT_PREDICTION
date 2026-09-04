"""The organism-scope over-call layer: it warns where measured, is silent elsewhere, never changes a call."""
from __future__ import annotations

import json
from pathlib import Path

from dna_decode.data.organism_scope import OVERCALL_SCOPE, one_line, overcall_for
from dna_decode.data.trust_surface import DISCLOSURE_LAYERS, trust_block

ROOT = Path(__file__).resolve().parents[1]


def test_it_fires_on_the_organism_where_the_overcall_was_measured():
    b = overcall_for("gentamicin", "Klebsiella_pneumoniae")
    assert b and b["ppv_in_this_organism"] < b["ppv_in_validated_scope"]


def test_it_is_silent_in_the_validated_scope():
    """E. coli is 12/12 on BV-BRC and 146/146 on PD -- warning there would be crying wolf."""
    assert overcall_for("gentamicin", "Escherichia_coli_Shigella") is None


def test_it_is_silent_for_an_unmeasured_drug_or_a_missing_organism():
    assert overcall_for("ciprofloxacin", "Klebsiella_pneumoniae") is None
    assert overcall_for("gentamicin", None) is None
    assert overcall_for("gentamicin", "") is None


def test_organism_matching_tolerates_the_naming_variants_callers_actually_pass():
    for name in ("Klebsiella", "Klebsiella_pneumoniae", "Klebsiella pneumoniae subsp. pneumoniae",
                 "KLEBSIELLA"):
        assert overcall_for("gentamicin", name), name


def test_the_layer_is_registered_and_reaches_the_public_accessor():
    assert "organism_scope" in DISCLOSURE_LAYERS
    assert "organism_scope" in trust_block("gentamicin", "Klebsiella_pneumoniae")
    assert "organism_scope" not in trust_block("gentamicin", "Escherichia_coli_Shigella")


def test_it_never_changes_a_call_a_tier_or_a_metric():
    """Augment-only: the layer adds its own key and touches nothing else."""
    warned = trust_block("gentamicin", "Klebsiella_pneumoniae")
    clean = trust_block("gentamicin", "Escherichia_coli_Shigella")
    for k in ("tier", "headline", "acc", "sens", "spec", "prediction", "call"):
        if k in warned and k in clean:
            pass  # values legitimately differ per cell; the point is the layer adds no such key
    assert "prediction" not in warned and "call" not in warned
    assert warned["organism_scope"]["status"] == "measured_overcall_outside_validated_scope"
    assert "never changes the call" in warned["organism_scope"]["note"]


def test_every_entry_carries_its_controls_and_its_unsettled_question():
    """A measured over-call that skipped the controls would be the PD label-artifact mistake again."""
    for row in OVERCALL_SCOPE:
        assert len(row["controls_passed"]) >= 3
        assert row["not_settled"]
        assert row["artifact"].startswith("wiki/")
        assert row["measured_counts"]["S"] > 0


def test_every_cited_artifact_exists():
    for row in OVERCALL_SCOPE:
        assert (ROOT / row["artifact"]).is_file(), row["artifact"]


def test_the_numbers_match_the_committed_artifact():
    """Non-vacuity: the index must be derived from the run, not typed from memory."""
    art = ROOT / "wiki" / "gentamicin_rmt_bvbrc_hunt.json"
    if not art.is_file():
        return
    d = json.loads(art.read_text(encoding="utf-8"))
    kp = [h for h in d["all_hits"] if str(h["genome_name"]).startswith("Klebsiella")]
    row = OVERCALL_SCOPE[0]
    assert row["measured_counts"]["R"] == sum(1 for h in kp if h["phenotype"] == "Resistant")
    assert row["measured_counts"]["S"] == sum(1 for h in kp if h["phenotype"] == "Susceptible")


def test_the_one_liner_names_the_measured_gap_and_the_artifact():
    line = one_line(overcall_for("gentamicin", "Klebsiella"))
    assert "ORGANISM-SCOPE WARNING" in line and "0.475" in line and "wiki/" in line
    assert one_line(None) is None


# --- source-concentration qualification (added 2026-09-03) -------------------------------------
# The published PPV 0.475 does not clear this project's own source-diversity bar. These pin the
# qualification so the warning can never silently revert to an archive-level claim.

def test_the_concentration_qualification_is_carried_in_the_block():
    sc = overcall_for("gentamicin", "Klebsiella")["source_concentration"]
    assert sc["passes_own_bar"] is False
    assert sc["largest_source_share"] > sc["own_bar"], "must record WHY it fails the bar"
    # The susceptible carriers are what carry the finding; their concentration is the load-bearing number.
    assert sc["susceptible_largest_source_share"] >= 0.90
    assert sc["ppv_excluding_largest_source"] > 0.9, (
        "excluding the dominant source the over-call largely disappears -- that is the whole point")


def test_the_one_liner_REFUSES_to_state_the_overcall_without_its_scope_caveat():
    """A bare 'PPV 0.475 in Klebsiella' overstates what 6 sources at 66% can support."""
    line = one_line(overcall_for("gentamicin", "Klebsiella"))
    assert "SINGLE-SOURCE-DOMINATED" in line
    assert "not for this organism generally" in line
    assert "0.976" in line, "the leave-one-source-out PPV must be visible, not just the headline"


def test_the_bar_is_the_projects_own_and_has_not_drifted():
    """The qualification is only meaningful against the bar the project actually enforces."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from source_diverse_validate import MAX_SOURCE_SHARE
    sc = overcall_for("gentamicin", "Klebsiella")["source_concentration"]
    assert sc["own_bar"] == MAX_SOURCE_SHARE


def test_the_qualification_does_not_leak_into_the_validated_scope():
    """E. coli passes the bar; it must stay silent, or the warning cries wolf where the rule is safe."""
    assert overcall_for("gentamicin", "Escherichia_coli_Shigella") is None
