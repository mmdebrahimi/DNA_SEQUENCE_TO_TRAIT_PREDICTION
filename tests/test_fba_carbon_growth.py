"""FBA carbon-source growth: pure name/exchange helpers + real-model medium-swap smoke (slow)."""
from __future__ import annotations

import pytest

from dna_decode.fba.carbon_growth import (
    match_carbon_exchange,
    normalize_carbon_name,
)


def test_normalize_carbon_name_alias():
    assert normalize_carbon_name("D-Maltose") == "maltose"
    assert normalize_carbon_name("a-Ketoglutaric") == "2-oxoglutarate"
    # unaliased names pass through lowercased
    assert normalize_carbon_name("D-Glucose") == "d-glucose"


def test_match_carbon_exchange():
    idx = {"d-glucose": "EX_glc__D_e", "maltose": "EX_malt_e"}
    assert match_carbon_exchange("D-Glucose", idx) == "EX_glc__D_e"       # exact (lowercased)
    assert match_carbon_exchange("D-Maltose", idx) == "EX_malt_e"          # via alias -> maltose
    assert match_carbon_exchange("Unobtainium", idx) is None               # unresolved


@pytest.mark.slow
def test_predict_growth_real_iml1515():
    pytest.importorskip("cobra")
    from dna_decode.fba.carbon_growth import build_exchange_name_index, predict_growth
    from dna_decode.fba.model import load_model

    m = load_model()
    idx = build_exchange_name_index(m)
    # glucose grows (~0.877); a pentose (arabinose) grows but slower than glucose
    glc = predict_growth(m, idx["d-glucose"])
    ara = predict_growth(m, idx["l-arabinose"])
    assert glc > 0.5
    assert ara > 1e-4 and ara < glc                    # pentoses grow slower -> quantitative + ordered
    # a carbon source with no exchange -> 0.0 (no crash)
    assert predict_growth(m, "EX_not_a_real_metabolite_e") == 0.0
