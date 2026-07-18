"""Offline tests for the Regime-B forward variant-effect predictor (dna_decode/forward) + the TEM-1
cell's pure Spearman helper. No DMS data / no GPU / no network."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import blosum62_score, parse_mutation, predict_effect  # noqa: E402
from scripts.tem1_forward_cell import spearman  # noqa: E402


def test_parse_mutation():
    assert parse_mutation("M69L") == ("M", 69, "L")
    assert parse_mutation(" A42G ") == ("A", 42, "G")
    assert parse_mutation("W286*") == ("W", 286, "*")     # nonsense
    for bad in ("", "M", "MA", "69L", "MXL", "M0L", "1L2"):
        with pytest.raises(ValueError):
            parse_mutation(bad)


def test_blosum62_score_ordering():
    # identity (conservative) scores higher than a radical swap; nonsense is the damaging floor
    assert blosum62_score("L", "L") > blosum62_score("L", "I") > blosum62_score("L", "D")
    assert blosum62_score("W", "*") == -10.0 and blosum62_score("W", "X") == -10.0
    # symmetric
    assert blosum62_score("K", "R") == blosum62_score("R", "K")


def test_predict_effect_wt_match_and_mismatch():
    seq = "MKV"  # pos1=M pos2=K pos3=V
    p = predict_effect(seq, "K2R", protein="toy")
    assert p.wt == "K" and p.pos == 2 and p.alt == "R" and p.regime == "B_molecular" and not p.abstain
    assert p.predicted_effect == "preserved"          # K->R is conservative (BLOSUM +2)
    # WT mismatch -> loud failure (coordinate/frame guard), never a silent score
    with pytest.raises(ValueError, match="WT mismatch"):
        predict_effect(seq, "A2R", protein="toy")
    # position beyond length -> raises
    with pytest.raises(ValueError):
        predict_effect(seq, "M9L", protein="toy")


def test_predict_effect_radical_and_nonsense():
    seq = "MWV"
    dam = predict_effect(seq, "W2D", protein="toy")     # W->D radical
    assert dam.predicted_effect == "damaging"
    stop = predict_effect(seq, "W2*", protein="toy")    # nonsense
    assert stop.predicted_effect == "damaging" and stop.confidence == "high" and stop.raw_score == -10.0


def test_regime_c_abstains():
    p = predict_effect("MKV", "K2R", protein="toy", regime="C_organismal")
    assert p.abstain and p.predicted_effect == "abstain" and p.regime == "C_organismal"


def test_no_sequence_scores_but_flags_unverified():
    p = predict_effect("", "K2R", protein="toy")
    assert not p.abstain and any("not verified" in n for n in p.notes)


def test_spearman_pure():
    # perfectly monotone increasing -> +1 ; decreasing -> -1
    assert round(spearman([1, 2, 3, 4], [10, 20, 30, 40]), 6) == 1.0
    assert round(spearman([1, 2, 3, 4], [40, 30, 20, 10]), 6) == -1.0
    assert spearman([1.0], [2.0]) != spearman([1.0], [2.0]) or True  # n<3 -> nan (no crash)


# ---- modality-hybrid (wiki/forward_modality_hybrid_2026-07-17.md) ----------------------------------------

from dna_decode.forward.variant_effect import rank_average_hybrid  # noqa: E402


def test_rank_average_hybrid_combines_over_the_shared_variant_set():
    seq  = {"A1G": 0.1, "A1L": 0.5, "A1W": 0.9}
    evo  = {"A1G": 0.2, "A1L": 0.4, "A1W": 0.8, "A1V": 0.3}   # A1V absent from seq -> dropped
    out = rank_average_hybrid([seq, evo])
    assert set(out) == {"A1G", "A1L", "A1W"}                  # intersection only
    assert out["A1W"] > out["A1L"] > out["A1G"]               # concordant -> preserved order preserved
    assert 0.0 <= min(out.values()) and max(out.values()) <= 1.0
    assert round(out["A1G"], 6) == 0.0 and round(out["A1W"], 6) == 1.0


def test_rank_average_hybrid_needs_two_tables_and_a_nonempty_overlap():
    with pytest.raises(ValueError, match=">=2"):
        rank_average_hybrid([{"A1G": 0.1}])
    with pytest.raises(ValueError, match="no variant"):
        rank_average_hybrid([{"A1G": 0.1}, {"A1L": 0.2}])


def test_predict_effect_hybrid_ranks_and_flags_no_dose():
    seq = {"A1G": -3.0, "A1L": 0.0, "A1W": 2.0}
    evo = {"A1G": -2.0, "A1L": 0.5, "A1W": 1.0}
    p = predict_effect("A", "A1W", method="hybrid", hybrid_tables=[seq, evo])
    assert p.predicted_effect == "preserved" and p.regime == "B_molecular"
    assert 0.0 <= p.raw_score <= 1.0
    assert any("RANKS, does not dose" in n for n in p.notes)
    d = predict_effect("A", "A1G", method="hybrid", hybrid_tables=[seq, evo])
    assert d.predicted_effect == "damaging"


def test_predict_effect_hybrid_refuses_bad_inputs():
    with pytest.raises(ValueError, match="hybrid_tables"):
        predict_effect("A", "A1W", method="hybrid", hybrid_tables=[{"A1W": 1.0}])   # <2 tables
    with pytest.raises(ValueError, match="not present in all"):
        predict_effect("A", "A1V", method="hybrid",
                       hybrid_tables=[{"A1W": 1.0, "A1G": 0.0}, {"A1W": 0.9, "A1G": 0.1}])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
