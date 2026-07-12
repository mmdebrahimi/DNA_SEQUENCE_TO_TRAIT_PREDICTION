"""Tests for the HIV blind-spot position-novelty flag (`scripts/hiv_blindspot_position_novelty.py`,
Family D). Uses the deployed HIV catalog (pure Python, CI-safe); real-data test SKIPS when the gitignored
NNRTI dataset is absent.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import hiv_blindspot_position_novelty as bs  # noqa: E402


def test_is_catalogued_known_efv_drm():
    assert bs._is_catalogued("efavirenz", "K103N") is True     # THE canonical EFV DRM


def test_is_catalogued_rejects_non_drm_substitution():
    # a silent/benign substitution far from any NNRTI DRM position must NOT be catalogued as R
    assert bs._is_catalogued("efavirenz", "K11R") is False


def test_prereg_constants_frozen():
    assert bs.SENS_MIN == 0.30 and bs.LIFT_MIN == 1.2


def test_verdict_thresholds_sane():
    # FLAG_RECOVERS_BLINDSPOT iff median(sens)>=SENS_MIN AND median(lift)>=LIFT_MIN — sane bar bounds
    assert 0.0 < bs.SENS_MIN <= 0.5 and bs.LIFT_MIN >= 1.0


def test_module_does_not_modify_frozen_bacterial_surface():
    src = (REPO / "scripts" / "hiv_blindspot_position_novelty.py").read_text(encoding="utf-8")
    assert "amr_rules" not in src and "calibrated_amr_rules" not in src   # reads hiv_amr, not the frozen bacterial surface


_NNRTI = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.txt"


@pytest.mark.skipif(not _NNRTI.exists(), reason="HIV NNRTI dataset absent (gitignored)")
def test_analyse_drug_counts_internally_consistent_real_data():
    rows = bs.load_rows(_NNRTI)
    m = bs.analyse_drug(rows, "efavirenz")
    assert m["n_blindspot_true_R"] >= 1                          # the documented EFV blind spot exists
    assert m["n_catalog_negative"] == m["n_blindspot_true_R"] + m["n_catneg_true_S"]
    assert 0.0 <= (m["flag_sens_on_blindspot"] or 0) <= 1.0


@pytest.mark.skipif(not _NNRTI.exists(), reason="HIV NNRTI dataset absent (gitignored)")
def test_run_verdict_valid():
    res = bs.run(_NNRTI)
    assert res["verdict"] in ("FLAG_RECOVERS_BLINDSPOT", "BLINDSPOT_NOT_POSITION_LOCAL", "NO_POWERED_DRUGS")
    assert "efavirenz" in res["per_drug"]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
