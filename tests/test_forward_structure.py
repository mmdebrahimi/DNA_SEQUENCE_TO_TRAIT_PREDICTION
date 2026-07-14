"""Offline tests for the structure-based method seam (dna_decode/forward/structure_scorer) + its
predict_effect wiring. No torch_geometric / torch_scatter / biotite on this host -> the real ESM-IF path
correctly raises StructureMethodUnavailable; the seam is validated with a synthetic esm_if_table."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import (  # noqa: E402
    StructureMethodUnavailable,
    alphafold_pdb_url,
    esm_if_tier,
    esm_if_variant_table,
    predict_effect,
)


def test_esm_if_tier():
    assert esm_if_tier(0.5) == "preserved"       # delta >= -1
    assert esm_if_tier(-1.0) == "preserved"
    assert esm_if_tier(-2.0) == "uncertain"
    assert esm_if_tier(-3.0) == "damaging"       # delta <= -3
    assert esm_if_tier(-9.0) == "damaging"


def test_alphafold_url():
    assert alphafold_pdb_url("P60484") == "https://alphafold.ebi.ac.uk/files/AF-P60484-F1-model_v4.pdb"


def test_predict_effect_esm_if_seam():
    # ESM-IF conditional-LL delta table (higher = structure-compatible = preserved)
    tab = {"V133A": -4.2, "K10R": 0.3, "M50T": -2.0}
    dmg = predict_effect("", "V133A", method="esm_if", esm_if_table=tab)
    assert dmg.method == "esm_if" and dmg.predicted_effect == "damaging" and abs(dmg.raw_score + 4.2) < 1e-9
    ben = predict_effect("", "K10R", method="esm_if", esm_if_table=tab)
    assert ben.predicted_effect == "preserved"
    assert predict_effect("", "M50T", method="esm_if", esm_if_table=tab).predicted_effect == "uncertain"
    with pytest.raises(ValueError, match="needs esm_if_table"):
        predict_effect("", "V133A", method="esm_if")


def test_esm_if_real_path_raises_when_deps_absent():
    """On this host (no torch_scatter) the real ESM-IF loader must fail LOUDLY + clearly, never silently."""
    with pytest.raises(StructureMethodUnavailable, match="structure stack is not installed"):
        esm_if_variant_table(Path("nonexistent.pdb"), "MKV", ["K2R"])


def test_leaderboard_knows_esm_if():
    from scripts.forward_leaderboard import _METHOD
    assert _METHOD["esm_if_structure"] == "esm_if"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
