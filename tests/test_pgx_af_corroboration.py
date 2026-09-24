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
    """Every tabled variant must land in band.

    The variant count is DERIVED from the table, not pinned to a literal: this test used to assert
    `n_variants == 5` and broke the moment the missing CYP4F2 EAS row was added (2026-09-24), which is
    the hardcoded-count drift pattern. What matters is that nothing is out of band and nothing was
    silently dropped, not that the table has a particular size.
    """
    rep = m.build_report([dict(r) for r in m.PGX_AF])
    assert rep["n_variants"] == len(m.PGX_AF) >= 5
    assert rep["n_in_band"] == rep["n_variants"]
    assert rep["verdict"] == "AF_CORROBORATED"


def test_cyp4f2_is_checked_in_both_asserted_populations():
    """REGRESSION GUARD for the 2026-09-24 defect: this table checked CYP4F2 in EUR only -- the number
    that is correct -- while cyp4f2.py asserted TWO populations, and the UNCHECKED one ('~79% EAS') was
    wrong by 0.58. A one-population-per-variant table cannot catch that by construction, so both must
    stay present. Panel-derived truth: EUR 0.2773, EAS 0.2128.
    """
    pops = {r["pop"] for r in m.PGX_AF if r["gene"] == "CYP4F2"}
    assert {"EUR", "EAS"} <= pops
    eas = [r for r in m.PGX_AF if r["gene"] == "CYP4F2" and r["pop"] == "EAS"][0]
    assert eas["af"] < 0.30, "the EAS *3 frequency is ~0.21; 0.79 was the REFERENCE-allele frequency"


def test_covers_the_four_new_genes():
    assert set(m.build_report([dict(r) for r in m.PGX_AF])["genes"]) == {"NUDT15", "UGT1A1", "CYP4F2", "ABCG2"}


def test_getrm_wall_and_tier_honest():
    rep = m.build_report([dict(r) for r in m.PGX_AF])
    assert "EXTERNAL_WALL" in rep["getrm_concordance_status"]
    assert "KNOWLEDGE_BASELINE" in rep["honesty_tier"]
    assert "NOT an independent per-sample concordance" in rep["honesty_tier"]
