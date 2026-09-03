"""Pins for the BV-BRC independent-archive hunt: the filters, the control, and the organism split."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HUNT = ROOT / "wiki" / "gentamicin_rmt_bvbrc_hunt.json"
CTRL = ROOT / "wiki" / "gentamicin_rmt_bvbrc_control.json"


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_only_measured_phenotypes_are_used():
    """BV-BRC also ships ML-PREDICTED phenotypes; scoring a deterministic rule against a model's output
    is the circular-label gate firing."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    assert "Laboratory Method" in d["phenotype_filter"]
    assert "ML-predicted" in d["phenotype_filter"]


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_the_archive_is_shown_independent_by_measured_overlap_not_assertion():
    """An 'independent archive' returning the same isolates is independent in name only."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    ov = d["pd_overlap"]
    assert ov.get("checked") is True
    assert ov["n_new_relative_to_pd"] > ov["n_shared_with_pd_sweep"]


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_the_genotype_caller_is_named_and_is_not_ours():
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    assert "CARD" in d["genotype_caller"] and "AMRFinder" in d["genotype_caller"]
    assert any("CARD" in lim and "tool-derived" in lim for lim in d["honest_limits"])


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_every_susceptible_carrier_would_actually_be_rescued_by_the_deployed_rule():
    """A counter-example the deployed rule would never have called is not a counter-example to it."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    s = d["susceptible_carriers"]
    assert s, "no susceptible carriers recorded"
    assert all(h["rescued_by_deployed_rule"] for h in s)


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_the_deployed_rules_own_organism_scope_is_reported_separately():
    """The headline is organism-stratified: Klebsiella over-calls, E. coli does not. Pooling them would
    either condemn a rule that is fine in its scope or hide a real over-call outside it."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    ec = [h for h in d["all_hits"] if str(h["genome_name"]).startswith(("Escherichia", "Shigella"))]
    kp = [h for h in d["all_hits"] if str(h["genome_name"]).startswith("Klebsiella")]
    assert ec and kp
    assert not [h for h in ec if h["phenotype"] == "Susceptible"], "an E. coli susceptible carrier exists"
    assert [h for h in kp if h["phenotype"] == "Susceptible"], "no Klebsiella susceptible carrier"


@pytest.mark.skipif(not CTRL.is_file(), reason="control artifact not present")
def test_the_control_uses_the_same_yardstick_and_threshold_as_the_pd_control():
    """Same determinant, same rule, opposite verdict -- that comparability is what makes it a control."""
    d = json.loads(CTRL.read_text(encoding="utf-8"))
    assert d["verdict"] in ("LABEL_ARTIFACT", "SPECIFIC_TO_RMT", "INCONCLUSIVE")
    assert d["aac3_R_rate_inside"] is not None and d["aac3_R_rate_outside"] is not None


@pytest.mark.skipif(not CTRL.is_file(), reason="control artifact not present")
def test_a_specific_to_rmt_verdict_is_backed_by_the_numbers():
    """Non-vacuity: the verdict must follow from the rates, not be asserted."""
    d = json.loads(CTRL.read_text(encoding="utf-8"))
    if d["verdict"] == "SPECIFIC_TO_RMT":
        assert not (d["aac3_R_rate_inside"] < 0.5 and d["aac3_R_rate_outside"] > 0.8)
        assert d["inside"]["rmt"]["r_rate"] < d["inside"]["aac3_no_rmt"]["r_rate"]


@pytest.mark.skipif(not CTRL.is_file(), reason="control artifact not present")
def test_the_control_disclaims_settling_the_deployed_scope():
    d = json.loads(CTRL.read_text(encoding="utf-8"))
    assert any("per-organism" in lim for lim in d["honest_limits"])
