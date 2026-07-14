"""Offline tests for the AlphaMissense method (dna_decode/forward/am_scorer) + its predict_effect wiring.
No D: / no network (synthetic AM dicts)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import am_table_for_mutants, am_tier, predict_effect  # noqa: E402


def test_am_tier_thresholds():
    assert am_tier(0.10) == "preserved"      # benign (<= 0.34)
    assert am_tier(0.34) == "preserved"
    assert am_tier(0.45) == "uncertain"      # ambiguous band
    assert am_tier(0.564) == "damaging"      # pathogenic (>= 0.564)
    assert am_tier(0.99) == "damaging"


def test_am_table_for_mutants_offset_and_coverage():
    am_by_variant = {"V133A": 0.9, "K10R": 0.2}          # UniProt-numbered
    # DMS numbering with offset +0 -> V133A covered, K10R covered, Q5P not covered
    t = am_table_for_mutants(am_by_variant, 0, ["V133A", "K10R", "Q5P", "bad", "A1B:C2D"])
    assert t == {"V133A": 0.9, "K10R": 0.2}
    # offset: DMS pos + offset = UniProt pos. DMS 'V33A' with offset 100 -> UniProt V133A
    t2 = am_table_for_mutants({"V133A": 0.9}, 100, ["V33A"])
    assert t2 == {"V33A": 0.9}               # keyed by the DMS mutation, value = AM of the UniProt variant


def test_predict_effect_alphamissense_polarity_and_tier():
    am_table = {"V133A": 0.9, "K10R": 0.2, "M50T": 0.45}
    dmg = predict_effect("", "V133A", method="alphamissense", am_table=am_table)
    assert dmg.method == "alphamissense" and dmg.predicted_effect == "damaging"
    assert abs(dmg.raw_score - (1.0 - 0.9)) < 1e-9        # raw_score = 1 - AM (higher = benign)
    ben = predict_effect("", "K10R", method="alphamissense", am_table=am_table)
    assert ben.predicted_effect == "preserved" and abs(ben.raw_score - 0.8) < 1e-9
    unc = predict_effect("", "M50T", method="alphamissense", am_table=am_table)
    assert unc.predicted_effect == "uncertain"


def test_predict_effect_alphamissense_missing_variant_raises():
    with pytest.raises(ValueError, match="not AlphaMissense-covered"):
        predict_effect("", "Q5P", method="alphamissense", am_table={"V133A": 0.9})
    with pytest.raises(ValueError, match="needs am_table"):
        predict_effect("", "V133A", method="alphamissense")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
