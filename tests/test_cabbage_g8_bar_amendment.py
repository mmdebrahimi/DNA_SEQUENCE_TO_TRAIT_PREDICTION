"""Guards for the 2026-10-19 amendment that admitted sourmash as a G8 distance source.

The amendment widens a FROZEN pre-commitment, which is exactly the move that can quietly
destroy a bar's value. These guards pin the four things the amendment promised not to do:
keep the forbidden-proxy clause verbatim, keep the bar imported by reference rather than
restated, cite real evidence, and not let the amendment inherit the original freeze's
"before any number existed" properties.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
BAR = REPO / "wiki" / "cabbage_g8_acceptance_bar.json"

# The clause that makes the bar non-circular: the decoder's prediction is a function of the
# same determinants, so a determinant-derived partition would manufacture agreement.
FORBIDDEN_PROXY = "The portal's pre-computed AMRFinderPlus tree is FORBIDDEN as a lineage proxy"


def _bar() -> dict:
    return json.loads(BAR.read_text(encoding="utf-8"))


def test_the_bar_parses_and_is_still_the_v1_schema():
    d = _bar()
    assert d["schema"] == "cabbage-g8-acceptance-bar-v1"


def test_the_forbidden_proxy_clause_SURVIVES_the_amendment():
    """Widening the distance METHOD must not drop the anti-circularity clause.

    Checked on the live `distance_source`, not merely somewhere in the file, because the
    clause only does work where the method is declared.
    """
    d = _bar()
    assert FORBIDDEN_PROXY in d["method_commitments"]["distance_source"]


def test_the_bar_VALUE_is_still_imported_by_reference_not_restated():
    """A literal threshold in this file could drift from the gate that owns the verdict."""
    d = _bar()
    src = d["the_bar_is_not_restated_here"]["source"]
    assert src == "dna_decode.eval.rejection_gates.MIN_EFFECTIVE_LINEAGES"
    # and the referenced symbol must actually exist, at the recorded value
    from dna_decode.eval.rejection_gates import MIN_EFFECTIVE_LINEAGES

    assert d["the_bar_is_not_restated_here"]["value_at_freeze"] == MIN_EFFECTIVE_LINEAGES


def test_the_amendment_cites_an_artifact_that_EXISTS_and_actually_passed():
    """An `agreement_evidence` block naming a missing or failing artifact is worse than none."""
    d = _bar()
    ev = d["agreement_evidence"]
    art = REPO / ev["artifact"]
    assert art.is_file(), f"cited evidence artifact missing: {ev['artifact']}"
    cc = json.loads(art.read_text(encoding="utf-8"))
    assert cc["status"] == "OK" and cc["exit_code"] == 0, "cited artifact did not pass"
    assert all(c["met"] is True for c in cc["clauses"]), "cited artifact has an unmet clause"
    # the recorded pair count must match the artifact, not be restated from prose
    assert ev["n_pairs"] == f'{cc["n_pairs_observed"]}/{cc["n_pairs_expected"]}'
    assert ev["cohorts_blocked"] == len(cc["blocked"]) == 0


def test_the_amendment_does_NOT_inherit_the_original_freeze_properties():
    """The two `frozen_before_*` booleans describe the 2026-10-17 freeze ONLY.

    By amendment time the crosscheck had run and 580 assemblies were fetched, so an
    amendment claiming to pre-date any number would be false. The scope note is what
    stops a reader inheriting them.
    """
    d = _bar()
    assert d["amended"] == "2026-10-19"
    of = d["original_freeze"]
    assert of["frozen"] == "2026-10-17", "original freeze date must be preserved"
    assert of["frozen_before_any_lineage_number_existed"] is True
    assert of["frozen_before_any_assembly_was_fetched"] is True
    note = of["scope_note"]
    assert "ORIGINAL" in note and "does NOT inherit" in note


def test_clause_1_is_not_presented_as_the_evidence():
    """Clause 1 is self-labelled nearly vacuous; the amendment must say so.

    Without this the amendment reads as '20/20 agreement', which is the one summary the
    bar's own text warns against.
    """
    ev = _bar()["agreement_evidence"]
    guidance = ev["read_clause_2_not_clause_1"]
    assert "VACUOUS" in guidance.upper()
    assert "clause 2" in guidance.lower()


def test_the_honest_limits_travel_WITH_the_amendment():
    """The limits live in the artifact; the amendment must carry them, not strip them."""
    ev = _bar()["agreement_evidence"]
    limits = ev["what_this_does_NOT_license"]
    assert isinstance(limits, list) and len(limits) >= 3
    blob = " ".join(limits)
    # the three load-bearing ones
    assert "FracMinHash" in blob, "the not-the-same-estimator limit must be carried"
    assert "CABBAGE" in blob, "the extrapolation-across-cohorts limit must be carried"
    assert "SINGLE TRIP" in blob.upper() or "single trip" in blob.lower()


def test_the_five_prospective_lock_files_are_NOT_touched_by_this_bar():
    """The G8 bar is deliberately NOT one of the sha256-pinned files.

    If it ever became one, amending it would retire every prospective number — so this
    asserts the separation rather than assuming it.
    """
    manifests = sorted((REPO / "wiki").glob("prospective_lock_manifest_*.json"))
    assert manifests, "no lock manifest found"
    for m in manifests:
        pinned = json.loads(m.read_text(encoding="utf-8")).get("surface_sha256", {})
        assert not any("cabbage_g8" in p for p in pinned), (
            f"{m.name} pins the G8 bar; amending it would retire the lock"
        )
