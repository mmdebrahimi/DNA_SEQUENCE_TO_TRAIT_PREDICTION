"""Offline pins for the multi-gene PGx AF-corroboration validator."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import pgx_af_corroboration as m  # noqa: E402


def test_classify_af_bands():
    assert m.classify_af(0.10, 0.05, 0.15) == "IN_BAND"
    assert m.classify_af(0.50, 0.05, 0.15) == "OUT_OF_BAND"
    assert m.classify_af(None, 0.0, 1.0) == "NO_DATA"


def test_all_new_cells_corroborate():
    rep = m.build_report([dict(r) for r in m.PGX_AF])
    assert rep["n_variants"] == 5
    assert rep["n_in_band"] == 5
    assert rep["verdict"] == "AF_CORROBORATED"


def test_covers_the_four_new_genes():
    assert set(m.build_report([dict(r) for r in m.PGX_AF])["genes"]) == {"NUDT15", "UGT1A1", "CYP4F2", "ABCG2"}


def test_getrm_wall_and_tier_honest():
    rep = m.build_report([dict(r) for r in m.PGX_AF])
    assert "EXTERNAL_WALL" in rep["getrm_concordance_status"]
    assert "KNOWLEDGE_BASELINE" in rep["honesty_tier"]
    assert "NOT an independent per-sample concordance" in rep["honesty_tier"]
