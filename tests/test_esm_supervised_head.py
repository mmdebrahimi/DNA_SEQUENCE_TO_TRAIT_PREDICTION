"""Offline pins for the supervised-head Phase-2 helpers (no torch / no model / no ESM)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.esm_supervised_resistance_head import _balacc, _blosum, _catalog_call  # noqa: E402


def test_blosum_symmetry_and_values():
    assert _blosum("A", "A") == 4.0          # BLOSUM62 diagonal
    assert _blosum("W", "W") == 11.0
    assert _blosum("A", "V") == 0.0          # A<->V
    assert _blosum("?", "X") == 0.0          # unknown -> 0 (no crash)


def test_balacc_math():
    # 3 true-R all called R, 3 true-S all called S -> perfect
    assert _balacc([True, True, True, False, False, False], [True, True, True, False, False, False]) == 1.0
    # all called S -> sens 0, spec 1 -> 0.5
    assert _balacc([True, True, False, False], [False, False, False, False]) == 0.5


def test_catalog_call_shape():
    # a known NNRTI DRM (K103N) -> R via the frozen dispatch; a benign RT position -> S
    r = _catalog_call({"gene": "RT", "drug": "efavirenz", "wt": "K", "pos": 103, "mut": "N"})
    s = _catalog_call({"gene": "RT", "drug": "efavirenz", "wt": "A", "pos": 400, "mut": "V"})
    assert r == "R" and s == "S"
    # NRTI M184V -> position-based R
    assert _catalog_call({"gene": "RT", "drug": "3tc" if False else "lamivudine",
                          "wt": "M", "pos": 184, "mut": "V"}) == "R"
