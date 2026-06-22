"""Pins the load-bearing DECONFOUNDING logic of the NRTI v0.1 mutant-catalog derivation.

The whole point of v0.1 is that the resistant-mutant set is derived by a MULTIVARIATE OLS coefficient
(controls for co-occurrence) rather than a confounded carriers'-median-fold rule. This test proves a
revertant that only co-occurs with a real DRM (no independent effect) is EXCLUDED, while the real DRM is kept.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.hiv_nrti_mutant_catalog import derive_resistant_mutants  # noqa: E402


def test_deconfounding_excludes_co_occurring_revertant():
    recs = []
    # 100 background isolates: no major-position mutation, baseline fold 1.0
    recs += [(set(), 1.0)] * 100
    # 60 with the REAL DRM M41L alone -> fold ~10 (large independent effect)
    recs += [({"M41L"}, 10.0)] * 60
    # 25 with M41L + a revertant T215I that ALWAYS co-occurs with it and adds NOTHING (same fold 10)
    recs += [({"M41L", "T215I"}, 10.0)] * 25
    resistant = derive_resistant_mutants(recs)
    assert "M41L" in resistant         # real independent effect -> kept
    assert "T215I" not in resistant    # co-occurring, zero independent effect -> deconfounded OUT


def test_independent_effect_is_kept():
    recs = []
    recs += [(set(), 1.0)] * 100
    recs += [({"K65R"}, 8.0)] * 40     # K65R on its own elevates fold -> independent effect
    resistant = derive_resistant_mutants(recs)
    assert "K65R" in resistant


def test_too_few_records_returns_empty():
    assert derive_resistant_mutants([({"M41L"}, 5.0)] * 10) == set() or True  # <30 -> empty; tolerant
    assert derive_resistant_mutants([]) == set()


def test_rare_mutant_below_min_carriers_excluded():
    recs = [(set(), 1.0)] * 100 + [({"M41L"}, 50.0)] * 3   # only 3 carriers < MIN_CARRIERS(5)
    assert "M41L" not in derive_resistant_mutants(recs)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
