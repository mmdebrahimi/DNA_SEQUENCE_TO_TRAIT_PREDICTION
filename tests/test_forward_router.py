"""Offline tests for the forward edit router (dna_decode/forward/router) — regime classification + routing.
No torch / no network (Regime-B uses a mock ESM table; Regime-A uses an injected determinant key set)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import (  # noqa: E402
    REGIME_A,
    REGIME_B,
    REGIME_C,
    REGIME_UNKNOWN,
    catalogue_call,
    classify_edit,
    predict_edit,
    variant_key,
)


def test_classify_edit_precedence():
    # organismal wins over everything
    r, _ = classify_edit("acsA", drug="glucose", determinant_locus=True, molecular_predictor=True,
                         organismal=True)
    assert r == REGIME_C
    # determinant locus + drug -> A
    assert classify_edit("gyrA", drug="ciprofloxacin", determinant_locus=True)[0] == REGIME_A
    # determinant locus WITHOUT drug -> not A (falls through)
    assert classify_edit("gyrA", determinant_locus=True)[0] == REGIME_UNKNOWN
    # molecular predictor -> B
    assert classify_edit("blaTEM-1", molecular_predictor=True)[0] == REGIME_B
    # nothing -> unknown
    assert classify_edit("yfjX")[0] == REGIME_UNKNOWN


def test_variant_key_and_catalogue_call():
    keys = {variant_key("rpoB", "S450L"), variant_key("gyrA", "S83L")}
    assert catalogue_call("rpoB", "s450l", keys)[0] == "R"       # case-insensitive
    assert catalogue_call("gyrA", "S83L", keys)[0] == "R"
    assert catalogue_call("gyrA", "D87N", keys)[0] == "S"        # not catalogued -> S
    assert catalogue_call("katG", "S315T", keys)[0] == "S"


def test_predict_edit_regime_a_catalogue():
    keys = {variant_key("rpoB", "S450L")}
    r = predict_edit("rpoB", "S450L", regime=REGIME_A, drug="rifampicin", resistance_keys=keys)
    assert r["prediction"] == "R" and r["predictor"] == "determinant_catalogue" and not r["abstain"]
    s = predict_edit("rpoB", "H445Y", regime=REGIME_A, drug="rifampicin", resistance_keys=keys)
    assert s["prediction"] == "S"


def test_predict_edit_regime_b_uses_predictor():
    # mock ESM table so no torch: K at pos 2 in "MKV"
    r = predict_edit("blaTEM-1", "K2R", regime=REGIME_B, protein_seq="MKV", method="esm2",
                     esm_table={2: {"K": -1.0, "R": -0.5}})
    assert r["prediction"] == "preserved" and r["predictor"] == "esm2_DMS_validated"
    assert abs(r["raw_score"] - 0.5) < 1e-9


def test_predict_edit_regime_c_and_unknown_abstain():
    c = predict_edit("acsA", "X10Y", regime=REGIME_C)
    assert c["abstain"] and c["prediction"] == "ABSTAIN" and c["predictor"] == "none"
    u = predict_edit("yfjX", "A5T")     # auto -> unknown
    assert u["abstain"] and u["regime"] == REGIME_UNKNOWN


def test_predict_edit_auto_routes():
    # no regime passed -> classify from flags; a resistance edit never reaches the likelihood predictor
    keys = {variant_key("gyrA", "S83L")}
    r = predict_edit("gyrA", "S83L", drug="ciprofloxacin", determinant_locus=True, resistance_keys=keys)
    assert r["regime"] == REGIME_A and r["prediction"] == "R" and r["predictor"] == "determinant_catalogue"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
