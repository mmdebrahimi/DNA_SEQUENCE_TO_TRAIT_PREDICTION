"""Pins the HIV INSTI v0.1 mutant-catalog derivation: the load-bearing DECONFOUNDING logic (synthetic,
CI-safe) + the real-data MVP criterion (v0.1 improves-or-holds balanced accuracy on all 5 INSTIs at the
uniform fold>=3 cutoff). The real-data test SKIPS when the gitignored INI_DataSet.txt is absent.

Mirrors tests/test_hiv_pi_mutant_catalog.py (the INSTI v0.1 builder is the PI builder's class twin).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.hiv_insti_mutant_catalog import DEFAULT_DATA, derive_resistant_mutants, run  # noqa: E402


# --- deconfounding logic (synthetic; no dataset needed) ---

def test_deconfounding_excludes_co_occurring_accessory():
    recs = []
    recs += [(set(), 1.0)] * 100                       # background, baseline fold
    recs += [({"Q148H"}, 10.0)] * 60                   # the REAL major DRM alone -> large independent effect
    recs += [({"Q148H", "G140S"}, 10.0)] * 25          # accessory G140S ALWAYS rides on Q148H here, adds nothing
    resistant = derive_resistant_mutants(recs)
    assert "Q148H" in resistant                        # independent effect -> kept
    assert "G140S" not in resistant                    # co-occurring, zero independent effect -> deconfounded OUT


def test_independent_effect_is_kept():
    recs = [(set(), 1.0)] * 100 + [({"N155H"}, 8.0)] * 40   # N155H on its own elevates fold
    assert "N155H" in derive_resistant_mutants(recs)


def test_too_few_records_returns_empty():
    assert derive_resistant_mutants([({"Q148H"}, 5.0)] * 10) == set()   # <30 -> empty
    assert derive_resistant_mutants([]) == set()


def test_rare_mutant_below_min_carriers_excluded():
    recs = [(set(), 1.0)] * 100 + [({"Q148R"}, 50.0)] * 3   # only 3 carriers < MIN_CARRIERS(5)
    assert "Q148R" not in derive_resistant_mutants(recs)


# --- real-data MVP criterion (gated on the gitignored INI dataset) ---

@pytest.mark.skipif(not DEFAULT_DATA.exists(), reason="INI_DataSet.txt not present (gitignored)")
def test_insti_v0_1_closes_the_gap_on_real_data():
    result = run(DEFAULT_DATA)
    assert result["n_drugs"] == 5                      # all 5 INSTIs scored
    gains = {d: m["balacc_gain_v0_1_minus_v0"] for d, m in result["per_drug"].items()
             if m.get("balacc_gain_v0_1_minus_v0") is not None}
    assert len(gains) == 5
    # v0.1 deconfounding improves-or-holds every INSTI at the uniform cutoff (the deferral was pessimistic)
    assert result["n_drugs_improved_or_held"] == 5
