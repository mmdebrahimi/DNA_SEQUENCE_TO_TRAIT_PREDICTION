"""Tests for the HIV supervised drug-panel generality run (offline; real-data run skip-guarded)."""
from __future__ import annotations

import importlib.util, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))
pytest.importorskip("sklearn")

spec = importlib.util.spec_from_file_location("hiv_supervised_panel", REPO / "scripts" / "hiv_supervised_panel.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
import dna_decode.data.hiv_amr as H


def test_panel_covers_all_11_rt_drugs():
    assert set(mod.NNRTI) == {"EFV", "NVP", "ETR", "RPV", "DOR"}
    assert set(mod.NRTI) == {"3TC", "ABC", "AZT", "D4T", "DDI", "TDF"}
    assert mod.NRTI["3TC"] == 5.0 and mod.NRTI["TDF"] == 1.5   # validated per-drug lower cutoffs


def test_catalog_calls_use_the_right_catalog():
    pcols = [f"P{p}" for p in range(1, 250)]
    # NNRTI mutant-level: K103N is a major DRM -> called
    r = {c: "-" for c in pcols}; r["P103"] = "N"
    assert mod._nnrti_call(r, pcols) is True
    # NNRTI: K103S is NOT in the major-DRM mutant set -> the blind spot (not called by NNRTI catalog)
    r2 = {c: "-" for c in pcols}; r2["P103"] = "S"
    assert mod._nnrti_call(r2, pcols) is False
    # NRTI position-based: M184V is at major position 184 -> called
    r3 = {c: "-" for c in pcols}; r3["P184"] = "V"
    assert mod._nrti_call(r3, pcols) is True and 184 in H.NRTI_MAJOR_POSITIONS
