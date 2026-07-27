"""Tests for the forward cell's multi-mutant (additive-null) capability."""
import pytest
from dna_decode.forward import predict_multi_effect, predict_effect

# synthetic controlled protein: WT residue = SEQ[pos-1], so mutations verify by RAW position
# (avoids Ambler-style gapped numbering). pos: 1M 2A 3K 4E 5L 6V 7G 8R 9S 10T 11Q 12W 13F 14Y 15H 16N 17D 18C 19P 20I
SEQ = "MAKELVGRSTQWFYHNDCPI" * 2


def test_additive_sum_matches_singles():
    muts = ["A2G", "K3D", "E4V"]
    m = predict_multi_effect(SEQ, muts, protein="synth", method="blosum62")
    singles = [predict_effect(SEQ, x, protein="synth", method="blosum62").raw_score for x in muts]
    assert abs(m.additive_score - sum(singles)) < 1e-9
    assert abs(m.mean_score - sum(singles) / 3) < 1e-9
    assert m.epistasis_model == "additive_null" and m.regime == "B_molecular"
    assert len(m.per_variant) == 3


def test_single_damaging_dominates():
    # G7W: glycine->tryptophan is a large BLOSUM-negative (damaging) substitution
    g7w = predict_effect(SEQ, "G7W", protein="synth", method="blosum62")
    assert g7w.predicted_effect == "damaging"                      # precondition
    m = predict_multi_effect(SEQ, ["A2G", "G7W"], protein="synth", method="blosum62")
    assert m.predicted_effect == "damaging"                        # one damaging edit dominates


def test_conflicting_positions_raise():
    with pytest.raises(ValueError, match="distinct positions"):
        predict_multi_effect(SEQ, ["A2G", "A2K"], protein="synth", method="blosum62")


def test_wt_mismatch_propagates():
    with pytest.raises(ValueError, match="WT mismatch"):
        predict_multi_effect(SEQ, ["M2L"], protein="synth", method="blosum62")  # pos2 is A not M


def test_regime_c_abstains():
    m = predict_multi_effect(SEQ, ["A2G"], protein="x", regime="C_organismal")
    assert m.abstain and m.predicted_effect == "abstain"


def test_empty_raises():
    with pytest.raises(ValueError, match=">=1 mutation"):
        predict_multi_effect(SEQ, [], protein="synth")
