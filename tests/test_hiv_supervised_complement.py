"""Tests for the HIV supervised blind-spot complement scorer — NNRTI + PI + INSTI (uses committed model JSONs)."""
from __future__ import annotations

import pytest

from dna_decode.data import hiv_supervised_complement as C


def test_registry_covers_three_classes():
    assert set(C.SUPPORTED_CLASSES) == {"NNRTI", "PI", "INSTI"}


@pytest.mark.parametrize("cls,gene", [("NNRTI", "RT"), ("PI", "PR"), ("INSTI", "IN")])
def test_each_model_loads_v1(cls, gene):
    info = C.model_info(cls)
    assert info["schema"] == "hiv-supervised-complement-v1" and info["gene"] == gene
    assert info["n_features"] > 100


def test_norm_accepts_multiple_token_forms():
    assert C._norm({"103N"}) == {"103N"}
    assert C._norm({"K103N"}) == {"103N"}          # strips WT letter
    assert C._norm({(103, "N")}) == {"103N"}       # tuple form


# Per-class high-risk genotype. NNRTI + INSTI are LOW genetic barrier (a single major DRM confers resistance);
# PI is HIGH genetic barrier (single mutation insufficient — resistance needs STACKED mutations). The model
# correctly encodes this: single V82F ~0.33, stacked PI DRMs -> resistant.
@pytest.mark.parametrize("cls,resistant", [
    ("NNRTI", {"103N"}), ("INSTI", {"155H"}), ("PI", {"82F", "54V", "47A", "84V"})])
def test_class_resistant_genotype_scores_high_empty_low(cls, resistant):
    r_res = C.blind_spot_risk(resistant, drug_class=cls)
    r_none = C.blind_spot_risk(set(), drug_class=cls)
    assert r_res > 0.7 and r_none < 0.5
    assert C.is_flagged(resistant, drug_class=cls) is True and C.is_flagged(set(), drug_class=cls) is False


def test_pi_high_genetic_barrier_single_mutation_moderate():
    # a single PI major DRM is NOT sufficient for resistance (high barrier) — the model encodes it
    assert C.blind_spot_risk({"82F"}, drug_class="PI") < 0.5
    assert C.blind_spot_risk({"82F", "54V", "47A", "84V"}, drug_class="PI") > 0.9


def test_default_class_is_nnrti():
    assert C.blind_spot_risk({"103N"}) == C.blind_spot_risk({"103N"}, drug_class="NNRTI")


def test_unknown_tokens_contribute_zero():
    assert C.blind_spot_risk({"103N"}, "NNRTI") == pytest.approx(
        C.blind_spot_risk({"103N", "999Z"}, "NNRTI"))


def test_unknown_class_and_missing_model_raise(tmp_path):
    with pytest.raises(C.ComplementUnavailable):
        C.blind_spot_risk({"103N"}, drug_class="BOGUS")
    with pytest.raises(C.ComplementUnavailable):
        C.blind_spot_risk({"103N"}, model_path=str(tmp_path / "nope.json"))
