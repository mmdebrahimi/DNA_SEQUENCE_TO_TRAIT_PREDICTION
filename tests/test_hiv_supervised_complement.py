"""Tests for the HIV NNRTI supervised blind-spot complement scorer (offline; uses the committed model JSON)."""
from __future__ import annotations

import math
import pytest

from dna_decode.data import hiv_supervised_complement as C


def test_model_loads_and_is_v1():
    info = C.model_info()
    assert info["schema"] == "hiv-nnrti-supervised-complement-v1"
    assert info["n_features"] > 100 and info["drug_trained"] == "EFV"


def test_norm_accepts_multiple_token_forms():
    assert C._norm({"103N"}) == {"103N"}
    assert C._norm({"K103N"}) == {"103N"}          # strips WT letter
    assert C._norm({(103, "N")}) == {"103N"}       # tuple form


def test_major_drm_scores_high_polymorphism_low():
    # K103N is the archetypal NNRTI DRM -> high risk; an empty genotype -> low (just the intercept)
    r_drm = C.blind_spot_risk({"103N"})
    r_none = C.blind_spot_risk(set())
    assert 0.0 <= r_none <= 1.0 and 0.0 <= r_drm <= 1.0
    assert r_drm > 0.8 and r_none < 0.5
    assert C.is_flagged({"103N"}) is True and C.is_flagged(set()) is False


def test_accessory_stack_raises_risk_over_single():
    # accessory mutations should ADD risk (the weighted-combination point) — two < ... monotone in known DRMs
    base = C.blind_spot_risk({"179D"})
    stacked = C.blind_spot_risk({"179D", "98G", "227L"})
    assert stacked >= base


def test_unknown_tokens_contribute_zero():
    # a residue not in the feature set must not change the score
    assert C.blind_spot_risk({"103N"}) == pytest.approx(C.blind_spot_risk({"103N", "999Z"}))


def test_missing_model_raises(tmp_path):
    with pytest.raises(C.ComplementUnavailable):
        C.blind_spot_risk({"103N"}, model_path=str(tmp_path / "nope.json"))
