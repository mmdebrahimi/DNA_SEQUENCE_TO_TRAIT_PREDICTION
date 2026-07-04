"""Offline pins for the Phase-3 imputation helper (no zip / no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.abo_blood import call_abo_o_status  # noqa: E402
from scripts.impute_determinant_abstain import _impute_map  # noqa: E402


def test_impute_map_majority_per_tag():
    # tag 'GG' -> target mostly DD; tag 'AA' -> target mostly II; tag 'AG' -> mostly DI
    pairs = [("DD", "GG"), ("DD", "GG"), ("DI", "GG"),       # GG -> DD (2 vs 1)
             ("II", "AA"), ("II", "AA"),                     # AA -> II
             ("DI", "AG"), ("DI", "AG"), ("DD", "AG")]       # AG -> DI (2 vs 1)
    m = _impute_map(pairs)
    assert m["GG"] == "DD" and m["AA"] == "II" and m["AG"] == "DI"


def test_imputed_ostatus_robust_to_di_ii_confusion():
    # the O-status rule collapses DI/II both to non-O, so a DI<->II imputation error does NOT flip O-status
    assert call_abo_o_status("DI") == call_abo_o_status("II") == "non-O"
    assert call_abo_o_status("DD") == "O"
    # only a DD<->(DI/II) confusion flips the O-status (the error mode that matters)
    assert call_abo_o_status("DD") != call_abo_o_status("DI")


def test_impute_map_empty():
    assert _impute_map([]) == {}
