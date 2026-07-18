"""Offline tests for the modality-hybrid sweep's pure analysis helpers (no ProteinGym data / no GPU).

Pins the per-phenotype-category bucketing (`by_category` / `paired_in_subset`) + the sign test + the
rank-average hybrid combine that the sweep constructs, on synthetic per-assay tables.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.forward_modality_hybrid_sweep import (  # noqa: E402
    BASELINE, MIN_ASSAYS_PER_CATEGORY, by_category, hybrid_spearman, paired_in_subset, sign_test_p,
)


def test_sign_test_two_sided():
    assert round(sign_test_p(5, 5), 3) == 1.0        # even split -> not significant
    assert sign_test_p(10, 0) < 0.01                 # a clean sweep is significant
    assert sign_test_p(0, 0) == 1.0                   # no data -> 1.0


def test_paired_in_subset_restricts_to_the_given_assays():
    per_assay = {
        "a": {BASELINE: 0.40, "GEMME": 0.50},        # +0.10 win
        "b": {BASELINE: 0.60, "GEMME": 0.55},        # -0.05 loss
        "c": {BASELINE: 0.50, "GEMME": 0.50},        # tie
        "z": {BASELINE: 0.10, "GEMME": 0.90},        # NOT in the subset -> ignored
    }
    r = paired_in_subset(per_assay, "GEMME", {"a", "b", "c"})
    assert r["n_paired"] == 3
    assert r["win_rate"] == round(1 / 3, 3)          # only 'a' wins; 'c' is a tie (not a win)


def test_by_category_buckets_and_gates_thin_categories():
    # 9 Stability assays (reportable) + 3 Binding (below the floor) -> Binding gated out
    per_assay, cats = {}, {}
    for i in range(9):
        per_assay[f"s{i}"] = {BASELINE: 0.50, "ProSST-2048": 0.60}
        cats[f"s{i}"] = "Stability"
    for i in range(3):
        per_assay[f"b{i}"] = {BASELINE: 0.50, "ProSST-2048": 0.70}
        cats[f"b{i}"] = "Binding"
    out = by_category(per_assay, cats)
    assert out["Stability"]["reportable"] is True and out["Stability"]["n_assays"] == 9
    assert out["Binding"]["reportable"] is False and out["Binding"]["n_assays"] == 3
    assert MIN_ASSAYS_PER_CATEGORY == 8
    # ProSST beats baseline in every Stability assay -> median delta +0.10, win-rate 1.0
    c = out["Stability"]["candidates"]["ProSST-2048"]
    assert c["median_delta"] == 0.10 and c["win_rate"] == 1.0


def test_hybrid_rank_average_orients_and_combines():
    # two concordant modalities over 20 variants (>= MIN_N); hybrid should track them
    rows = [{"mut": f"m{i}", "ESM2_650M": str(i), "GEMME": str(i * 2), "DMS_score": str(i)}
            for i in range(20)]
    res = hybrid_spearman(rows, ("ESM2_650M", "GEMME"))
    assert res is not None
    rho, n = res
    assert n == 20 and round(rho, 3) == 1.0           # perfectly concordant -> |Spearman| 1.0


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
