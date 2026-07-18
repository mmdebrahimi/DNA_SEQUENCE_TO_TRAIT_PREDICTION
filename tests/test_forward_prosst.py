"""Offline seam tests for the ProSST structure method (dna_decode/forward/prosst_scorer.py).

The real ProSST forward needs `transformers` + `prosst` (+ torch_geometric) and a GPU — NOT exercised in CI;
the real run is `scripts/prosst_lift.py` on Kaggle (Step 6 of the plan). These pin the seam: the
unavailable-signal, the predict_effect('prosst') branch on a stub table, the tier thresholds, the reused
AlphaFold URL, the hybrid composition, and the prosst_lift pure analysis helpers.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import prosst_tier, prosst_variant_table, quantize_structure  # noqa: E402  (exports)
from dna_decode.forward.prosst_scorer import (  # noqa: E402
    StructureMethodUnavailable, alphafold_pdb_url, prosst_model_name,
)
from dna_decode.forward.variant_effect import predict_effect, rank_average_hybrid  # noqa: E402


# ---- the unavailable-signal (this host has no prosst / torch_geometric) -----------------------------------

def test_quantize_structure_raises_when_prosst_absent():
    with pytest.raises(StructureMethodUnavailable, match="quantizer unavailable"):
        quantize_structure("nonexistent.pdb")


def test_load_prosst_raises_when_transformers_absent(monkeypatch):
    monkeypatch.setitem(sys.modules, "transformers", None)   # force `from transformers import ...` to fail
    from dna_decode.forward import prosst_scorer
    prosst_scorer._BUNDLE.clear()
    with pytest.raises(StructureMethodUnavailable, match="transformers"):
        prosst_scorer._load_prosst(2048)


def test_prosst_variant_table_needs_tokens_or_pdb():
    # reaches the input guard before any model load
    with pytest.raises((ValueError, StructureMethodUnavailable)):
        prosst_variant_table("MK", ["M1L"])


# ---- tiers + reused helpers -------------------------------------------------------------------------------

def test_prosst_tier_thresholds():
    assert prosst_tier(0.0) == "preserved"       # >= -1.0
    assert prosst_tier(-5.0) == "damaging"        # <= -3.0
    assert prosst_tier(-2.0) == "uncertain"       # in between


def test_alphafold_url_and_model_name():
    assert alphafold_pdb_url("P12345").endswith("AF-P12345-F1-model_v4.pdb")
    assert prosst_model_name(2048) == "AI4Protein/ProSST-2048"
    assert prosst_model_name(512) == "AI4Protein/ProSST-512"


# ---- the predict_effect seam (stub table -> no model) -----------------------------------------------------

def test_predict_effect_prosst_on_a_stub_table():
    p = predict_effect("M", "M1W", method="prosst", prosst_table={"M1W": -4.0})
    assert p.predicted_effect == "damaging" and p.regime == "B_molecular" and p.method == "prosst"
    q = predict_effect("M", "M1L", method="prosst", prosst_table={"M1L": 0.2})
    assert q.predicted_effect == "preserved"


def test_predict_effect_prosst_refuses_missing_table_or_variant():
    with pytest.raises(ValueError, match="prosst_table"):
        predict_effect("M", "M1W", method="prosst")                       # no table
    with pytest.raises(ValueError, match="prosst_table"):
        predict_effect("M", "M1W", method="prosst", prosst_table={"M1L": 0.1})   # variant absent


def test_unknown_method_message_names_prosst():
    with pytest.raises(NotImplementedError, match="prosst"):
        predict_effect("M", "M1L", method="magic")


def test_prosst_table_composes_into_the_hybrid():
    esm = {"M1L": -1.0, "M1W": -3.0, "M1V": 0.0}
    prosst = {"M1L": -0.5, "M1W": -4.0, "M1V": 0.1}
    combined = rank_average_hybrid([esm, prosst])
    assert set(combined) == {"M1L", "M1W", "M1V"}
    assert combined["M1V"] > combined["M1L"] > combined["M1W"]     # concordant order preserved


# ---- prosst_lift.py pure analysis helpers -----------------------------------------------------------------

def test_prosst_lift_pure_helpers():
    from scripts.prosst_lift import paired_in_subset, by_category, reproduction_median
    recs = [
        {"dms": "A_Activity", "status": "OK", "hybrid_minus_esm": 0.05, "prosst_vs_pg_repro": 0.9},
        {"dms": "B_Stability", "status": "OK", "hybrid_minus_esm": 0.07, "prosst_vs_pg_repro": 0.8},
        {"dms": "C_Stability", "status": "OK", "hybrid_minus_esm": -0.01, "prosst_vs_pg_repro": 0.85},
        {"dms": "D_bad", "status": "ERROR"},
    ]
    cats = {"A_Activity": "Activity", "B_Stability": "Stability", "C_Stability": "Stability"}
    p = paired_in_subset(recs)
    assert p["n"] == 3 and p["win"] == "2/3"
    bc = by_category(recs, cats)
    assert bc["Stability"]["n"] == 2 and bc["Activity"]["n"] == 1
    assert reproduction_median(recs) == 0.85


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
