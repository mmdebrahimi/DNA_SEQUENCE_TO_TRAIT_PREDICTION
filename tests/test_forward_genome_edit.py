"""Offline tests for the genome-level forward edit path (dna_decode/forward/genome_edit) + the ESM-method
wiring in predict_effect (mock table — no torch/transformers needed)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import (  # noqa: E402
    cds_point_edit,
    predict_effect,
    predict_genome_edit,
    translate_codon,
)


def test_translate_codon():
    assert translate_codon("ATG") == "M" and translate_codon("atg") == "M"
    assert translate_codon("TAA") == "*" and translate_codon("TGA") == "*" and translate_codon("TAG") == "*"
    assert translate_codon("GAA") == "E"
    for bad in ("AT", "ATGG", "ATX", ""):
        with pytest.raises(ValueError):
            translate_codon(bad)


def test_cds_point_edit_consequences():
    cds = "ATGAAAGTT"  # M K V
    # nt_pos 4 (codon2 start) A->C : AAA -> CAA = K->Q (missense)
    info = cds_point_edit(cds, 4, "A", "C")
    assert info["aa_pos"] == 2 and info["wt_aa"] == "K" and info["alt_aa"] == "Q" and info["within_codon"] == 0
    # nt_pos 6 (codon2 3rd base) A->G : AAA -> AAG = K->K (synonymous)
    syn = cds_point_edit(cds, 6, "A", "G")
    assert syn["wt_aa"] == "K" and syn["alt_aa"] == "K"
    # REF mismatch -> loud failure
    with pytest.raises(ValueError, match="REF mismatch"):
        cds_point_edit(cds, 4, "G", "C")
    with pytest.raises(ValueError):
        cds_point_edit(cds, 99, "A", "C")


def test_predict_genome_edit_missense_silent_nonsense():
    cds = "ATGGAAGTT"          # M E V
    prot = "MEV"
    # missense: nt_pos5 (codon2 2nd base) A->G : GAA -> GGA = E->G
    mis = predict_genome_edit(cds, 5, "A", "G", protein_seq=prot, protein="toy")
    assert mis.consequence == "missense" and mis.aa_mutation == "E2G"
    assert mis.protein_prediction is not None and mis.protein_prediction.wt == "E"
    # silent: nt_pos6 (codon2 3rd base) A->G : GAA -> GAG = E->E
    sil = predict_genome_edit(cds, 6, "A", "G", protein_seq=prot, protein="toy")
    assert sil.consequence == "silent" and sil.aa_mutation is None and sil.protein_prediction is None
    # nonsense: nt_pos4 (codon2 1st base) G->T : GAA -> TAA = E->* (stop)
    non = predict_genome_edit(cds, 4, "G", "T", protein_seq=prot, protein="toy")
    assert non.consequence == "nonsense" and non.alt_aa == "*" and non.aa_mutation == "E2*"
    assert non.protein_prediction.predicted_effect == "damaging"


def test_predict_genome_edit_double_coordinate_check():
    # translated WT AA must also match protein_seq (predict_effect re-verifies) -> a wrong protein_seq raises
    cds = "ATGGAAGTT"          # M E V
    with pytest.raises(ValueError, match="WT mismatch"):
        predict_genome_edit(cds, 5, "A", "G", protein_seq="MQV", protein="toy")  # says Q at pos2, really E


def test_predict_effect_esm2_mock_table():
    seq = "MKV"
    table = {2: {"K": -1.2, "R": -0.4, "A": -7.0}}    # ESM zero-shot log-probs at position 2
    p = predict_effect(seq, "K2R", method="esm2", esm_table=table)
    assert p.method == "esm2" and abs(p.raw_score - 0.8) < 1e-9 and p.predicted_effect == "preserved"
    pd = predict_effect(seq, "K2A", method="esm2", esm_table=table)
    assert abs(pd.raw_score + 5.8) < 1e-9 and pd.predicted_effect == "damaging"
    with pytest.raises(ValueError, match="requires esm_table"):
        predict_effect(seq, "K2R", method="esm2")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
