"""Offline tests for the AMR Portal independent-scoring pure logic (Wilson CI + main.tsv reconstruction +
confusion), plus a faithfulness check that the reconstructed main.tsv drives the FROZEN rule correctly."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.amr_portal_score_independent import (  # noqa: E402
    confusion, genotype_to_main_tsv, wilson_ci,
)


def test_wilson_ci():
    assert wilson_ci(0, 0) is None
    lo, hi = wilson_ci(50, 100)
    assert 0.39 < lo < 0.41 and 0.59 < hi < 0.61      # 50/100 ~ [0.40, 0.60]
    lo2, hi2 = wilson_ci(100, 100)
    assert hi2 > 0.999 and lo2 < 1.0


def test_confusion():
    assert confusion("R", "R") == "TP" and confusion("R", "S") == "FP"
    assert confusion("S", "S") == "TN" and confusion("S", "R") == "FN"
    assert confusion("INDETERMINATE", "R") is None


def test_genotype_to_main_tsv_columns_and_method():
    dets = [
        {"amr_element_symbol": "gyrA_S83L", "gene_symbol": "gyrA", "class": "QUINOLONE",
         "subclass": "QUINOLONE", "element_subtype": "POINT", "evidence_description": "DNA gyrase"},
        {"amr_element_symbol": "blaCTX-M-15", "gene_symbol": "bla", "class": "BETA-LACTAM",
         "subclass": "CEPHALOSPORIN", "element_subtype": "AMR", "evidence_description": "ESBL"},
    ]
    txt = genotype_to_main_tsv(dets)
    header, *rows = txt.strip().split("\n")
    assert header.split("\t") == ["Element symbol", "Method", "Class", "Subclass", "Element name",
                                  "% Identity to reference"]
    assert rows[0].split("\t")[:4] == ["gyrA_S83L", "POINTX", "QUINOLONE", "QUINOLONE"]   # POINT -> POINTX
    assert rows[1].split("\t")[:4] == ["blaCTX-M-15", "EXACTX", "BETA-LACTAM", "CEPHALOSPORIN"]


def test_reconstructed_main_tsv_drives_frozen_rule(tmp_path):
    """End-to-end: a reconstructed main.tsv must produce the SAME frozen call as a hand-written one.

    cipro needs >=2 QRDR POINT mutations for R; ceftriaxone needs >=1 extended-spectrum bla."""
    from dna_decode.eval.amr_rules import call_resistance
    # two QRDR points -> cipro R
    cipro_dets = [
        {"amr_element_symbol": "gyrA_S83L", "gene_symbol": "gyrA", "class": "QUINOLONE",
         "subclass": "QUINOLONE", "element_subtype": "POINT", "evidence_description": ""},
        {"amr_element_symbol": "parC_S80I", "gene_symbol": "parC", "class": "QUINOLONE",
         "subclass": "QUINOLONE", "element_subtype": "POINT", "evidence_description": ""},
    ]
    p = tmp_path / "m.tsv"
    p.write_text(genotype_to_main_tsv(cipro_dets), encoding="utf-8")
    assert call_resistance(p, "ciprofloxacin", organism=None)["prediction"] == "R"
    # one QRDR point only -> below threshold 2 -> S
    p.write_text(genotype_to_main_tsv(cipro_dets[:1]), encoding="utf-8")
    assert call_resistance(p, "ciprofloxacin", organism=None)["prediction"] == "S"
    # extended-spectrum bla -> ceftriaxone R; a plain blaTEM-1 (narrow) -> ceftriaxone S
    p.write_text(genotype_to_main_tsv([{"amr_element_symbol": "blaCTX-M-15", "class": "BETA-LACTAM",
                                        "subclass": "CEPHALOSPORIN", "element_subtype": "AMR",
                                        "gene_symbol": "bla", "evidence_description": ""}]), encoding="utf-8")
    assert call_resistance(p, "ceftriaxone", organism=None)["prediction"] == "R"
    p.write_text(genotype_to_main_tsv([{"amr_element_symbol": "blaTEM-1", "class": "BETA-LACTAM",
                                        "subclass": "BETA-LACTAM", "element_subtype": "AMR",
                                        "gene_symbol": "bla", "evidence_description": ""}]), encoding="utf-8")
    assert call_resistance(p, "ceftriaxone", organism=None)["prediction"] == "S"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
