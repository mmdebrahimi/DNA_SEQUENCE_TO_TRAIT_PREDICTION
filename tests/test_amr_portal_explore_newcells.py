"""Pins the exploratory new-cell verdict thresholds (Tier 3/4 cipro-QRDR transfer).

The verdict is the load-bearing honesty surface: it decides which new cells are clean enough to be
promotion CANDIDATES vs need organism-specific curation. Pure thresholds, no parquet/network.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.amr_portal_explore_newcells import _verdict  # noqa: E402


def test_clean_transfer_needs_high_balacc_and_both_axes():
    assert _verdict(0.85, 0.99) == "TRANSFERS_CLEAN"      # Neisseria-like
    assert _verdict(0.94, 0.89) == "TRANSFERS_CLEAN"      # Staph-like


def test_partial_when_one_axis_drags_balacc_mid():
    # Pseudomonas-like: high spec, low sens (efflux under-call) -> balacc 0.78 -> PARTIAL
    assert _verdict(0.577, 0.991) == "TRANSFERS_PARTIAL"


def test_needs_curation_when_recall_collapses():
    # Enterobacter/Serratia-like: spec 1.0 but sens ~0 -> conserved rule does not transfer
    assert _verdict(0.0, 1.0) == "NEEDS_CURATION"
    assert _verdict(0.264, 1.0) == "NEEDS_CURATION"       # balacc 0.632 < 0.70


def test_high_balacc_but_one_axis_below_floor_is_not_clean():
    # balacc 0.85 but sens 0.70 exactly clean; sens 0.69 -> not clean (drops to partial)
    assert _verdict(0.69, 1.0) == "TRANSFERS_PARTIAL"      # balacc 0.845 >=0.70 but sens<0.70 -> not clean
    assert _verdict(0.70, 1.0) == "TRANSFERS_CLEAN"        # both axes >=0.70 + balacc 0.85


def test_unscorable_when_a_class_empty():
    assert _verdict(None, 1.0) == "UNSCORABLE"
    assert _verdict(0.9, None) == "UNSCORABLE"
