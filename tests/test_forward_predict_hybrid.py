"""Offline tests for the flagship convenience predictor predict_variant_hybrid + the ESM pos->variant adapter."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import (  # noqa: E402
    esm_pos_table_to_variant_table, predict_variant_hybrid,
)


def test_esm_pos_table_to_variant_table_converts_and_deltas():
    # seq "MA": pos1 M, pos2 A. logp tables -> deltas logP(alt)-logP(wt)
    pos_table = {1: {"M": 0.0, "L": -1.0, "W": -3.0}, 2: {"A": 0.0, "G": -0.5}}
    v = esm_pos_table_to_variant_table(pos_table, "MA")
    assert v["M1L"] == -1.0 and v["M1W"] == -3.0 and v["A2G"] == -0.5
    assert "M1M" not in v                                  # self-sub excluded
    assert all(k[0] in "MA" for k in v)


def test_predict_variant_hybrid_ranks_in_saturation_context():
    # 3 variants; ESM + ProSST both rank W1P worst, W1L best -> hybrid agrees; percentile reflects it
    esm = {"W1L": 2.0, "W1V": 0.0, "W1P": -4.0}
    prosst = {"W1L": 1.5, "W1V": 0.1, "W1P": -3.0}
    worst = predict_variant_hybrid("W", "W1P", esm_table=esm, prosst_table=prosst)
    best = predict_variant_hybrid("W", "W1L", esm_table=esm, prosst_table=prosst)
    assert worst.predicted_effect == "damaging" and best.predicted_effect == "preserved"
    assert worst.method == "hybrid_esm2_prosst" and worst.regime == "B_molecular"
    assert worst.raw_score < best.raw_score
    assert any("percentile" in n for n in worst.notes)
    assert any("RANKS within THIS protein" in n for n in worst.notes)


def test_predict_variant_hybrid_accepts_esm_pos_table():
    # esm given as a {pos:{aa:logp}} table -> auto-converted
    esm_pos = {1: {"W": 0.0, "L": 1.0, "P": -4.0, "V": -0.5}}
    prosst = {"W1L": 1.0, "W1P": -3.0, "W1V": 0.0}
    p = predict_variant_hybrid("W", "W1L", esm_table=esm_pos, prosst_table=prosst)
    assert 0.0 <= p.raw_score <= 1.0 and p.predicted_effect in ("preserved", "uncertain", "damaging")


def test_predict_variant_hybrid_wt_mismatch_fails_loudly():
    with pytest.raises(ValueError, match="WT mismatch"):
        predict_variant_hybrid("K", "W1L", esm_table={"W1L": 1.0, "W1P": -1.0},
                               prosst_table={"W1L": 1.0, "W1P": -1.0})


def test_predict_variant_hybrid_missing_variant_fails():
    with pytest.raises(ValueError, match="not present in all"):
        predict_variant_hybrid("W", "W1Y", esm_table={"W1L": 1.0, "W1P": -1.0},
                               prosst_table={"W1L": 1.0, "W1P": -1.0})   # W1Y absent


def test_predict_variant_hybrid_three_way_with_extra_table():
    esm = {"W1L": 2.0, "W1P": -4.0, "W1V": 0.0}
    prosst = {"W1L": 1.5, "W1P": -3.0, "W1V": 0.1}
    gemme = {"W1L": 0.5, "W1P": -2.0, "W1V": -0.1}
    p = predict_variant_hybrid("W", "W1L", esm_table=esm, prosst_table=prosst, extra_tables=[gemme])
    assert p.predicted_effect == "preserved"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
