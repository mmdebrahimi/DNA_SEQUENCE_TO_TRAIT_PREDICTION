"""Pins the HIV PI v0.1 mutant-catalog derivation: the load-bearing DECONFOUNDING logic (synthetic, CI-safe)
+ the real-data MVP criterion (the v0.1 refinement closes the OLS-vs-catalog gap -> majority of PI drugs
improve-or-hold balanced accuracy at the uniform fold>=3 cutoff). The real-data test SKIPS when the
gitignored PI_DataSet.txt is absent (so CI stays green without it).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.hiv_pi_mutant_catalog import DEFAULT_DATA, derive_resistant_mutants, run  # noqa: E402


# --- deconfounding logic (synthetic; no dataset needed) ---

def test_deconfounding_excludes_co_occurring_accessory():
    recs = []
    recs += [(set(), 1.0)] * 100                       # background, baseline fold
    recs += [({"L90M"}, 10.0)] * 60                    # the REAL major DRM alone -> large independent effect
    recs += [({"L90M", "M46I"}, 10.0)] * 25            # accessory M46I ALWAYS rides on L90M, adds nothing
    resistant = derive_resistant_mutants(recs)
    assert "L90M" in resistant                         # independent effect -> kept
    assert "M46I" not in resistant                     # co-occurring, zero independent effect -> deconfounded OUT


def test_independent_effect_is_kept():
    recs = [(set(), 1.0)] * 100 + [({"I84V"}, 8.0)] * 40   # I84V on its own elevates fold
    assert "I84V" in derive_resistant_mutants(recs)


def test_too_few_records_returns_empty():
    assert derive_resistant_mutants([({"L90M"}, 5.0)] * 10) == set()   # <30 -> empty
    assert derive_resistant_mutants([]) == set()


def test_rare_mutant_below_min_carriers_excluded():
    recs = [(set(), 1.0)] * 100 + [({"L90M"}, 50.0)] * 3   # only 3 carriers < MIN_CARRIERS(5)
    assert "L90M" not in derive_resistant_mutants(recs)


# --- real-data MVP criterion (gated on the gitignored PI dataset) ---

@pytest.mark.skipif(not DEFAULT_DATA.exists(), reason="PI_DataSet.txt not present (gitignored)")
def test_pi_v0_1_closes_the_gap_on_real_data():
    result = run(DEFAULT_DATA)
    assert result["n_drugs"] == 8                      # all 8 PIs scored
    gains = {d: m["balacc_gain_v0_1_minus_v0"] for d, m in result["per_drug"].items()
             if m.get("balacc_gain_v0_1_minus_v0") is not None}
    assert len(gains) == 8
    n_improved = sum(1 for g in gains.values() if g >= 0)
    assert n_improved >= 5, f"only {n_improved}/8 PI drugs improve-or-hold: {gains}"   # majority MVP bar
    assert result["mean_balacc_gain"] >= 0.0           # delta-honest: the refinement does not regress overall


@pytest.mark.skipif(not DEFAULT_DATA.exists(), reason="PI_DataSet.txt not present (gitignored)")
def test_pi_v0_1_catalog_recovers_real_major_drms():
    """The data-derived catalog should surface canonical major PI DRMs (sanity that it learned biology)."""
    result = run(DEFAULT_DATA)
    nfv = set(result["deliverable_catalog_all_data"]["nelfinavir"])
    assert "D30N" in nfv and "L90M" in nfv              # the two signature nelfinavir DRMs
    drv = set(result["deliverable_catalog_all_data"]["darunavir"])
    assert "I50V" in drv                                # the signature darunavir DRM


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        try:
            fn(); print(f"PASS {fn.__name__}")
        except Exception as e:  # pragma: no cover
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns)} tests")
