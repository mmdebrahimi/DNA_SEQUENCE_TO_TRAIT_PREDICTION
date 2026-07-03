"""Offline pins for the ESM catalog-extension test's pure helpers (no torch / no model download)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.esm_catalog_extension_test import (  # noqa: E402
    BENIGN_NEGATIVES, _auc, _load_variants, _mpro_wt, _VAR,
)


def test_variant_regex():
    assert _VAR.match("140A").groups() == ("140", "A")
    assert _VAR.match("54A").groups() == ("54", "A")
    assert _VAR.match("P132H") is None            # WT-carrying form not accepted (position+mut only)
    assert _VAR.match("21I").groups() == ("21", "I")


def test_auc_direction():
    assert _auc([2.0, 3.0], [0.0, 1.0]) == 1.0     # positives all above negatives
    assert _auc([0.0, 1.0], [2.0, 3.0]) == 0.0     # all below
    assert _auc([1.0, 1.0], [1.0]) == 0.5          # ties
    assert _auc([1.0], []) is None                 # empty -> None


def test_load_variants_uses_reference_wt(tmp_path):
    wt = "ACDEFGHIKLMNPQRSTVWY" * 8                # 160-residue synthetic protein
    vj = tmp_path / "v.json"
    vj.write_text(json.dumps({"per_mutation_fold": {
        "10K": ["5", "7"], "20W": ["100"], "999A": ["3"], "bad": ["1"], "15G": ["not_num"]}}))
    got = _load_variants(vj, wt)
    d = {(v["pos"], v["mut"]): v for v in got}
    assert (10, "K") in d and d[(10, "K")]["wt"] == wt[9] and d[(10, "K")]["median_fold"] == 6.0
    assert (20, "W") in d and d[(20, "W")]["median_fold"] == 100.0
    assert (999, "A") not in d                     # out of range dropped
    assert (15, "G") not in d                      # all-non-numeric folds dropped


def test_mpro_reference_translates_to_catalytic_dyad():
    ref = Path(__file__).resolve().parent.parent / "data/sarscov2_ref/SARSCoV2_Mpro_NC045512_cds.fna"
    if not ref.exists():
        import pytest
        pytest.skip("Mpro reference not present")
    wt = _mpro_wt(ref)
    assert len(wt) == 306
    assert wt[40] == "H" and wt[144] == "C" and wt[165] == "E"   # H41/C145 catalytic + E166 nirmatrelvir


def test_benign_negatives_are_sourced_polymorphisms():
    assert ("P", 132, "H") in BENIGN_NEGATIVES and ("K", 90, "R") in BENIGN_NEGATIVES
