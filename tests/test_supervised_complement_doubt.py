"""Guards for the supervised blind-spot complement wired into the L2 doubt layer (2026-09-30).

WHY THIS EXISTS. `dna_decode/data/hiv_supervised_complement.py` shipped 2026-07-12 with leave-one-STUDY-out
blind-spot AUROC 0.81 NNRTI / 0.89 PI / 0.89 INSTI -- and was reachable from NOWHERE: imported only by its
own builder, zero mentions in `cell_registry`, no CLI route, not wired into the doubt layer whose entire job
is flagging catalog incompleteness. Third "shipped but unreachable" instance in this repo. These tests pin
the wiring AND the three ways I got it wrong while doing it.
"""

from __future__ import annotations

import pytest

from dna_decode.data import hiv_supervised_complement as C
from dna_decode.data.hiv_amr import call_hiv_observed
from dna_decode.eval.doubt import (STRONG, _complement_token, doubt_one_line,
                                   supervised_complement_signal, target_site_doubt)


def test_the_signal_fires_STRONG_on_a_real_blind_spot_isolate():
    """K103S is the case that motivated the wiring: a REAL Stanford isolate carrying only this substitution,
    which the deployed catalog calls SUSCEPTIBLE (it carries 103N, not 103S) while the measured wet-lab fold
    is log10 +0.833, about 6.8x -- genuinely resistant. Before this wiring nothing escalated it."""
    sig = supervised_complement_signal(["K103S"], "efavirenz", call="S")
    assert sig.tier == STRONG, sig.reason
    assert sig.evidence["assessed"] is True
    assert sig.evidence["blind_spot_risk"] >= C.DEFAULT_THRESHOLD


def test_three_states_are_never_collapsed():
    """not-measured (class outside SUPPORTED_CLASSES) / not-assessable (no observed subs) / assessed.
    Reporting "no doubt" for either of the first two would be a false clean bill of health."""
    unsupported = supervised_complement_signal(["M184V"], "lamivudine")      # NRTI
    assert unsupported.evidence["assessed"] is False
    assert unsupported.evidence["applicable"] is False
    assert "NOT measured" in unsupported.reason

    not_assessable = supervised_complement_signal([], "efavirenz")
    assert not_assessable.evidence["assessed"] is False
    assert not_assessable.evidence["applicable"] is True
    assert "not assessable" in not_assessable.reason.lower()

    assessed = supervised_complement_signal(["K103N"], "efavirenz", call="S")
    assert assessed.evidence["assessed"] is True


def test_it_NEVER_carries_a_call_value_and_the_guard_is_what_caught_me():
    """MY OWN BUG, pinned. The first version put `call_supplied: "R"` in the evidence for context, and
    `assert_no_call` refused it -- correctly: it checks VALUES, not just keys, and echoing the call into its
    own disclosure makes the disclosure call-shaped. The CLI was broken until this was fixed."""
    blk = target_site_doubt("efavirenz", {"RT": {"K103N"}}, call="R").as_dict()   # must not raise
    flat = repr(blk)
    assert "call_was_supplied" in flat and "call_supplied'" not in flat
    sig = [s for s in blk["signals"] if s["kind"] == "supervised_complement"][0]
    assert sig["evidence"]["call_was_supplied"] is True
    assert sig["evidence"]["blind_spot_question_applies"] is False
    for v in sig["evidence"].values():                 # no bare call token anywhere in the evidence
        assert v not in ("R", "S", "I", "RESISTANT", "SUSCEPTIBLE")


def test_the_blind_spot_question_does_not_arise_on_an_R_call():
    """The complement scores the risk that a SUSCEPTIBLE call is wrong. On an R call there is no blind spot,
    so the signal must report NONE with that reason rather than a spurious strong doubt."""
    r = supervised_complement_signal(["K103N"], "efavirenz", call="R")
    s = supervised_complement_signal(["K103N"], "efavirenz", call="S")
    assert r.tier != STRONG and "does not arise" in r.reason
    assert s.tier == STRONG                            # same genotype, unknown/S call -> fires
    assert r.evidence["blind_spot_risk"] == s.evidence["blind_spot_risk"]


def test_low_risk_is_reported_as_NOT_reassurance():
    """179F is the one CONFIRMED gap and the complement scores it near baseline (measured 0.180) while the
    purity screen flags it decisively. A quiet complement must therefore never read as a clean bill."""
    sig = supervised_complement_signal(["V179F"], "efavirenz", call="S")
    assert sig.tier != STRONG
    assert "NOT reassurance" in sig.reason
    assert sig.evidence["low_risk_is_not_reassurance"] is True


def test_the_token_forms_still_agree_so_normalisation_stays_redundant():
    """I documented `_complement_token` as LOAD-BEARING and that was WRONG -- measured, the complement
    normalises internally and both forms score identically. This pins that: if it ever stops, the
    normalisation becomes load-bearing and we learn it from a FAILURE, not from a silently quiet signal."""
    assert _complement_token("K103N") == "103N"
    assert _complement_token("103N") == "103N"
    raw = C.blind_spot_risk({"K103N"}, drug_class="NNRTI")
    norm = C.blind_spot_risk({"103N"}, drug_class="NNRTI")
    base = C.blind_spot_risk(set(), drug_class="NNRTI")
    assert raw == pytest.approx(norm, abs=1e-9), "the complement stopped normalising -- see the docstring"
    assert norm - base > 0.1, "K103N must move the score well off baseline, else the probe is vacuous"


def test_the_signal_is_AUGMENT_ONLY_and_can_only_raise_the_tier():
    """L2 qualifies, never alters. The two deterministic signals must be byte-identical with the third
    present, and adding a signal can raise `max_tier` but never lower it."""
    import json
    for obs in ({"RT": {"K103N"}}, {"RT": {"V179F"}}, {"RT": {"V179I"}}):
        blk = target_site_doubt("efavirenz", obs, call="S").as_dict()
        sigs = blk["signals"]
        assert [s["kind"] for s in sigs[:2]] == ["position_novelty", "target_site_completeness"]
        assert sigs[2]["kind"] == "supervised_complement"
        others = [s for s in sigs if s["kind"] != "supervised_complement"]
        max_without = next((t for t in (STRONG, "weak") if any(s["tier"] == t for s in others)), "none")
        if max_without == STRONG:
            assert blk["max_tier"] == STRONG, "adding a signal LOWERED the reported tier"
        json.dumps(blk)                                 # must stay serialisable


def test_the_same_genotype_cannot_render_differently_for_two_CALLERS():
    """MY OWN BUG, caught by an EXISTING test (`test_the_only_honest_silence_is_assessed_and_quiet`).

    The first wiring treated an unsupplied `call` as "assume susceptible" and FIRED. The CLI passes
    `call.prediction` and so stayed correctly quiet on K103N -- a catalogued major DRM the catalog itself
    calls RESISTANT -- while a library caller passing no call got a STRONG blind-spot doubt on the very same
    genotype. `target_site_doubt` now resolves the call from the same deployed catalog, so the two agree.
    """
    for gene_subs in ({"RT": {"K103N"}}, {"RT": {"K103S"}}, {"RT": {"V179F"}}, {"RT": {"Y181C"}}):
        derived = target_site_doubt("efavirenz", gene_subs)
        explicit = target_site_doubt("efavirenz", gene_subs,
                                     call=call_hiv_observed("efavirenz", gene_subs).prediction)
        assert derived.as_dict() == explicit.as_dict(), (
            f"{gene_subs} renders differently with and without an explicit call -- the caller-divergence "
            "bug is back")

    # And the direction that matters: on a genotype the catalog calls R, the honest answer is SILENCE.
    r_case = target_site_doubt("efavirenz", {"RT": {"K103N"}})
    assert call_hiv_observed("efavirenz", {"RT": {"K103N"}}).prediction == "R"
    assert doubt_one_line(r_case.as_dict()) is None, "fired a blind-spot doubt on an already-R call"


def test_an_unresolvable_call_is_NOT_ASSESSABLE_rather_than_a_doubt():
    """The fallback must not become the old bug in a new place. With no call AND no way to derive one, the
    signal reports the third state -- not a fire, and not a clean bill either."""
    sig = supervised_complement_signal(["K103S"], "efavirenz", call=None)
    assert sig.tier != STRONG
    assert sig.evidence["assessed"] is False and sig.evidence["applicable"] is True
    assert "not assessable" in sig.reason.lower()
    assert "NOT an absence of doubt" in sig.reason


def test_the_human_facing_line_can_actually_show_this_signal():
    """A block carried only in JSON is not a disclosure -- the standard this repo already records. On the
    K103S case the complement is the signal that fires, so it must be what `doubt_one_line` renders."""
    blk = target_site_doubt("efavirenz", {"RT": {"K103S"}}, call="S").as_dict()
    line = doubt_one_line(blk)
    assert line and "complement" in line.lower()
    assert "IN-DISTRIBUTION" in line or "in-distribution" in line
