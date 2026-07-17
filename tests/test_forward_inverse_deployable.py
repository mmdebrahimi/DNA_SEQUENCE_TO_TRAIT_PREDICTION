"""Tests for the degeneracy gate + the deployable (rank/percentile) inverse.

The gate is the load-bearing part: a censored assay does not fail loudly, it FLATTERS the metric. CcdB
scored the best err/span in the whole magnitude sweep (+77.1% vs null) purely because 79.3% of its
variants tie at one ceiling value, so most quantile targets collapse onto it.
"""
from __future__ import annotations

import pytest

from scripts.forward_inverse_deployable import TARGET_PCTS, random_null, rank_inverse
from scripts.forward_inverse_roundtrip import Candidate, assay_degeneracy


# ---- the degeneracy gate ---------------------------------------------------------------------------------

def test_a_censored_assay_is_flagged_degenerate():
    vals = [-2.0] * 80 + [-3.0, -4.0, -5.0, -6.0, -7.0] * 4
    d = assay_degeneracy(vals)
    assert d["degenerate"] is True
    assert d["mode_value"] == -2.0
    assert d["mode_share"] > 0.25
    assert "censored" in d["reason"]


def test_a_well_spread_assay_is_not_flagged():
    vals = [i / 100 for i in range(200)]
    d = assay_degeneracy(vals)
    assert d["degenerate"] is False
    assert d["n_distinct_values"] == 200


def test_too_few_distinct_levels_is_degenerate_even_without_a_dominant_mode():
    """A uniformly-binned assay has no dominant mode yet still cannot express a fine target grid."""
    vals = [float(i % 5) for i in range(500)]        # 5 levels, each 20% -> no mode >25%
    d = assay_degeneracy(vals)
    assert d["mode_share"] <= 0.25                    # the mode rule alone would MISS this
    assert d["degenerate"] is True                    # ...the distinct-levels rule catches it
    assert d["n_distinct_values"] == 5


def test_degeneracy_gate_catches_the_real_ccdb_assay():
    """Real-data pin: the exact assay that flattered the sweep must be excluded."""
    from pathlib import Path

    from scripts.forward_inverse_roundtrip import PG, load_substrate
    from scripts.forward_inverse_sweep import build_candidates

    if not PG.exists():
        pytest.skip("ProteinGym cache (D:) not mounted")
    target, dms, esm, cds = load_substrate("CCDB_ECOLI_Tripathi_2016", None)
    cands, _ = build_candidates(target, dms, cds)
    d = assay_degeneracy([c.measured for c in cands])
    assert d["degenerate"] is True
    assert d["n_distinct_values"] < 20
    assert d["mode_share"] > 0.5


def test_the_committed_sweep_excludes_ccdb_rather_than_scoring_it():
    import json
    from pathlib import Path

    art = Path(__file__).resolve().parent.parent / "wiki" / "forward_inverse_sweep_2026-07-17.json"
    if not art.exists():
        pytest.skip("sweep artifact not present")
    rows = {r["dms_id"]: r for r in json.loads(art.read_text(encoding="utf-8"))["assays"]}
    assert rows["CCDB_ECOLI_Tripathi_2016"]["status"] == "DEGENERATE_CENSORED_ASSAY"


# ---- the rank/percentile inverse -------------------------------------------------------------------------

def _cands(n=200):
    return [Candidate(f"A{i+1}V", i + 1, "A", "V", measured=-i / 10.0) for i in range(n)]


def test_a_perfect_ranker_lands_on_the_requested_percentile():
    c = _cands()
    r = rank_inverse(c, lambda x: -x.pos, TARGET_PCTS, top_k=1)   # perfectly rank-correlated
    assert r["mean_pct_err_top1"] < 0.02


def test_a_useless_ranker_is_no_better_than_the_null():
    c = _cands()
    r = rank_inverse(c, lambda x: 0.0, TARGET_PCTS, top_k=5)
    n = random_null(TARGET_PCTS, top_k=5)
    assert r["mean_pct_err_best_of_k"] >= n["mean_pct_err_best_of_k"] * 0.9


def test_rank_inverse_needs_no_calibrator_and_no_measured_labels_to_propose():
    """The whole point of the deployable variant: the PROPOSAL must depend only on scores. Perturbing the
    measured labels must not change which variant is proposed (it only changes the grading)."""
    a = rank_inverse(_cands(), lambda x: -x.pos, [0.5], top_k=1)
    shifted = [Candidate(c.mutant, c.pos, c.wt, c.alt, c.measured * 3.0 - 7.0) for c in _cands()]
    b = rank_inverse(shifted, lambda x: -x.pos, [0.5], top_k=1)
    assert a["per_target"][0]["proposed"] == b["per_target"][0]["proposed"]


def test_random_null_is_exact_and_deterministic():
    assert random_null(TARGET_PCTS, 5) == random_null(TARGET_PCTS, 5)


def test_random_null_top1_matches_the_closed_form():
    """E|U - p| = p^2 - p + 1/2 for U ~ Uniform(0,1). At p=0.5 that is 0.25."""
    n = random_null([0.5], top_k=1)
    assert n["mean_pct_err_top1"] == pytest.approx(0.25, abs=1e-6)


def test_best_of_k_null_beats_top1_null():
    n1 = random_null(TARGET_PCTS, 1)
    n5 = random_null(TARGET_PCTS, 5)
    assert n5["mean_pct_err_best_of_k"] < n1["mean_pct_err_top1"]
