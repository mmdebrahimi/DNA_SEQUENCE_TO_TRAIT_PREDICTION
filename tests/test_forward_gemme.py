"""Offline tests for the GEMME evolution seam (dna_decode/forward/gemme_scorer.py) + the 3-way hybrid.

GEMME's real run needs the JET2/R/Java toolchain (Windows-hostile) — NOT exercised in CI; the validation
path uses a precomputed column (GEMME is deterministic, so a precomputed column is canonical). These pin the
column adapter, the tier, the unavailable-signal, the predict_effect('gemme') branch, and that the 3-way
`rank_average_hybrid([esm2, prosst, gemme])` composes.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import GemmeUnavailable, gemme_table_from_column, gemme_tier  # noqa: E402  (exports)
from dna_decode.forward.gemme_scorer import run_gemme  # noqa: E402
from dna_decode.forward.variant_effect import predict_effect, rank_average_hybrid  # noqa: E402


def test_gemme_table_from_column_adapts_and_skips_na():
    rows = [
        {"mutant": "M1L", "GEMME": "-0.20"},
        {"mutant": "M1W", "GEMME": "-5.00"},
        {"mutant": "M1V", "GEMME": "NA"},        # skipped
        {"mutant": "M1A", "GEMME": ""},           # skipped
        {"mutant": "M1C"},                          # no GEMME key -> skipped
    ]
    t = gemme_table_from_column(rows)
    assert t == {"M1L": -0.20, "M1W": -5.00}


def test_gemme_table_from_column_raises_on_empty():
    with pytest.raises(ValueError, match="no usable"):
        gemme_table_from_column([{"mutant": "M1L", "GEMME": "NA"}])


def test_gemme_tier_thresholds():
    assert gemme_tier(0.0) == "preserved"        # >= -0.5
    assert gemme_tier(-6.0) == "damaging"         # <= -4.0
    assert gemme_tier(-2.0) == "uncertain"


def test_run_gemme_raises_when_toolchain_absent(monkeypatch):
    # toolchain absent -> GemmeUnavailable (never a silent wrong call). FORCE the absent condition by
    # patching shutil.which -> None so this is deterministic on ANY host (a Docker-present host would
    # otherwise skip the guard and hit a FileNotFoundError reading the nonexistent MSA).
    import dna_decode.forward.gemme_scorer as gs
    monkeypatch.setattr(gs.shutil, "which", lambda _cmd: None)
    with pytest.raises(GemmeUnavailable, match="Docker"):
        run_gemme("some.a3m", "MKAY")


def test_predict_effect_gemme_branch():
    p = predict_effect("M", "M1W", method="gemme", gemme_table={"M1W": -5.0})
    assert p.predicted_effect == "damaging" and p.method == "gemme" and p.regime == "B_molecular"
    with pytest.raises(ValueError, match="gemme_table"):
        predict_effect("M", "M1W", method="gemme")


def test_three_way_hybrid_composes():
    esm = {"A1G": -3.0, "A1L": 0.0, "A1W": 2.0}
    prosst = {"A1G": -4.0, "A1L": 0.1, "A1W": 1.5}
    gemme = {"A1G": -5.0, "A1L": -0.3, "A1W": 0.8}
    combined = rank_average_hybrid([esm, prosst, gemme])
    assert set(combined) == {"A1G", "A1L", "A1W"}
    assert combined["A1W"] > combined["A1L"] > combined["A1G"]      # all concordant


def test_three_way_lift_pure_helpers():
    from scripts.three_way_lift import paired, by_category
    recs = [
        {"dms": "A_Activity", "status": "OK", "three_minus_two": 0.02},
        {"dms": "B_Stability", "status": "OK", "three_minus_two": -0.01},
        {"dms": "C_Stability", "status": "OK", "three_minus_two": -0.02},
        {"dms": "D_bad", "status": "ERROR"},
    ]
    cats = {"A_Activity": "Activity", "B_Stability": "Stability", "C_Stability": "Stability"}
    p = paired(recs, "three_minus_two")
    assert p["n"] == 3 and p["win"] == "1/3"
    bc = by_category(recs, cats, "three_minus_two")
    assert bc["Stability"]["n"] == 2 and bc["Stability"]["win"] == "0/2"
    assert bc["Activity"]["win"] == "1/1"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
