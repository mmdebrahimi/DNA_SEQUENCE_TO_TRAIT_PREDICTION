"""Offline pins for the DPYD AF-corroboration validation (`scripts/dpyd_af_corroboration.py`)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import dpyd_af_corroboration as m  # noqa: E402


def test_classify_af_in_and_out_of_band():
    assert m.classify_af(0.005, 0.002, 0.015) == "IN_BAND"
    assert m.classify_af(0.10, 0.002, 0.015) == "OUT_OF_BAND"
    assert m.classify_af(0.0, 0.0, 0.005) == "IN_BAND"      # boundary inclusive
    assert m.classify_af(None, 0.0, 0.01) == "NO_DATA"


def test_committed_afs_all_corroborate():
    """The 4 actionable DPYD variants' committed EUR AFs all fall in the CPIC/gnomAD-expected band."""
    rep = m.build_report([dict(r) for r in m.DPYD_AF])
    assert rep["n_variants"] == 4
    assert rep["n_in_band"] == 4
    assert rep["verdict"] == "AF_CORROBORATED"


def test_getrm_wall_is_documented_not_faked():
    rep = m.build_report([dict(r) for r in m.DPYD_AF])
    assert "EXTERNAL_WALL" in rep["getrm_concordance_status"]
    assert "KNOWLEDGE_BASELINE" in rep["honesty_tier"]  # never overclaimed as independent concordance


def test_all_four_actionable_haplotypes_present():
    alleles = {r["allele"] for r in m.DPYD_AF}
    assert alleles == {"*2A", "*13", "c.2846A>T", "HapB3"}
