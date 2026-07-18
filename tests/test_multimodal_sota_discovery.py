"""Offline tests for the multimodal-SOTA discovery helpers (scripts/multimodal_sota_discovery.py).

Pins the sign test + the pure-column / hybrid Spearman helpers on synthetic rows (no ProteinGym data)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.multimodal_sota_discovery import _col_spearman, _hybrid_spearman, sign_p  # noqa: E402


def test_sign_p():
    assert round(sign_p(5, 5), 3) == 1.0
    assert sign_p(20, 0) < 0.01
    assert sign_p(0, 0) == 1.0


def _rows(n):
    # a monotone column that perfectly tracks DMS -> |Spearman| ~ 1.0
    return [{"mutant": f"m{i}", "COL": str(i), "A": str(i), "B": str(i * 2), "DMS_score": str(i)}
            for i in range(n)]


def test_col_spearman_perfect_and_underpowered():
    rows = _rows(25)
    dms = {r["mutant"]: float(r["DMS_score"]) for r in rows}
    assert round(_col_spearman(rows, "COL", dms), 3) == 1.0
    assert _col_spearman(_rows(5), "COL", {r["mutant"]: float(r["DMS_score"]) for r in _rows(5)}) is None  # < MIN_N


def test_hybrid_spearman_concordant_columns():
    rows = _rows(25)
    dms = {r["mutant"]: float(r["DMS_score"]) for r in rows}
    rho = _hybrid_spearman(rows, ("A", "B"), dms)   # both concordant with DMS
    assert rho is not None and round(rho, 3) == 1.0


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
