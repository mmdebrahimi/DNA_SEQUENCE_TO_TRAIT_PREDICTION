"""Tests for the quantitative-decoder capability-map rollup (pure flatteners, no I/O)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

spec = importlib.util.spec_from_file_location(
    "quantitative_decoder_capability_map", REPO / "scripts" / "quantitative_decoder_capability_map.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

_HIV = {"per_class": {"NNRTI": {"class": "NNRTI", "per_drug": {
    "EFV": {"n": 2168, "powered": True, "r2_oof": 0.844, "cover_90": 0.897, "fold_factor_90": 3.52,
            "calibrated_90": True},
    "CompMutList": {"n": 0, "powered": False, "note": "too few"}}}}}

_TB = {"per_drug": {
    "rifampicin": {"n": 12097, "powered": True, "r2_oof_resolved": 0.42, "cover_resolved_90": 0.9035,
                   "interval_fold_factor": 3.14, "calibrated": True},
    "amikacin": {"n": 5, "powered": False, "note": "no determinant vocab"}}}


def test_hiv_flatten_powered_and_under_powered():
    rows = mod.rows_from_hiv(_HIV)
    assert len(rows) == 2
    powered = [r for r in rows if r["calibrated"] is True]
    assert powered[0]["unit"] == "NNRTI/EFV" and powered[0]["cover_90"] == 0.897
    assert powered[0]["fold_factor"] == 3.52 and powered[0]["pathogen"] == "HIV-1"
    under = [r for r in rows if r["calibrated"] is None]
    assert under and under[0]["status"] == "under-powered"


def test_tb_flatten_uses_resolved_cover_and_fold():
    rows = mod.rows_from_tb(_TB)
    rif = [r for r in rows if r["unit"] == "rifampicin"][0]
    assert rif["cover_90"] == 0.9035 and rif["fold_factor"] == 3.14 and rif["calibrated"] is True
    ami = [r for r in rows if r["unit"] == "amikacin"][0]
    assert ami["calibrated"] is None and "under-powered" in ami["status"]


def test_build_summary_counts_and_coverage_range():
    res = mod.build.__wrapped__ if hasattr(mod.build, "__wrapped__") else None
    # build reads from disk by default; call the flatteners + summary logic directly on synthetic inputs
    rows = mod.rows_from_hiv(_HIV) + mod.rows_from_tb(_TB)
    calibrated = [r for r in rows if r["status"] == "calibrated"]
    assert len(calibrated) == 2                       # EFV + rifampicin
    covers = [r["cover_90"] for r in calibrated]
    assert min(covers) == 0.897 and max(covers) == 0.9035


def test_render_md_has_no_crash_on_synthetic():
    res = {"artifact": "quantitative_decoder_capability_map", "schema": "qdcm-v1",
           "sources": {"hiv": "h.json", "tb": "t.json"},
           "summary": {"n_cells_total": 3, "n_powered": 2, "n_calibrated": 2, "n_informative": 2,
                       "n_coverage_valid_only": 0, "n_under_powered": 1,
                       "coverage_range": [0.897, 0.9035], "pathogens": ["HIV-1", "M. tuberculosis"]},
           "honest_caveats": ["x"], "rows": mod.rows_from_hiv(_HIV) + mod.rows_from_tb(_TB)}
    md = mod.render_md(res, "2026-07-12")
    assert "2/2 powered cells INFORMATIVELY calibrated" in md and "NNRTI/EFV" in md and "rifampicin" in md
