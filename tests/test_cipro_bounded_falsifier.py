"""Tests for scripts/cipro_bounded_falsifier.py — verdict logic + helpers.

Pins the contract that Codex's runner on the Precision 7780 must satisfy if it
adopts Claude's draft. Specifically: positive-only ranking has noise for
non-positives, and consumers MUST gate on delta > 0.
"""
from __future__ import annotations

import math

import pytest


# ---- _ranked_by contract ----


def test_ranked_by_returns_unique_ranks_descending():
    from scripts.cipro_bounded_falsifier import _ranked_by
    deltas = [0.5, 0.3, 0.1, -0.2]
    ranks = _ranked_by(deltas, descending=True)
    assert ranks == [1, 2, 3, 4]


def test_ranked_by_handles_negative_descending():
    from scripts.cipro_bounded_falsifier import _ranked_by
    deltas = [-0.5, 0.0, 0.5]
    ranks = _ranked_by(deltas, descending=True)
    assert ranks == [3, 2, 1]


def test_ranked_by_non_positives_have_arbitrary_but_contiguous_ranks():
    """Contract: when non-positives are mapped to 0.0 + ranked descending,
    their ranks are arbitrary but contiguous-from-the-bottom. Consumers
    MUST gate on delta > 0 before reading the rank externally."""
    from scripts.cipro_bounded_falsifier import _ranked_by
    pos_replaced = [0.5, 0.0, 0.3, 0.0, -0.0]  # 0.5, 0.3 are real; rest mapped to 0.0
    ranks = _ranked_by(pos_replaced, descending=True)
    # Real positives get distinct top ranks
    assert ranks[0] == 1  # 0.5
    assert ranks[2] == 2  # 0.3
    # The three 0.0 entries get ranks 3, 4, 5 in some order
    zero_ranks = sorted([ranks[1], ranks[3], ranks[4]])
    assert zero_ranks == [3, 4, 5]


# ---- _logit ----


def test_logit_at_half_returns_zero():
    from scripts.cipro_bounded_falsifier import _logit
    assert _logit(0.5) == pytest.approx(0.0, abs=1e-12)


def test_logit_clamps_extremes():
    """logit(0) and logit(1) would be -inf / +inf; clamp via eps=1e-9."""
    from scripts.cipro_bounded_falsifier import _logit
    assert math.isfinite(_logit(0.0))
    assert math.isfinite(_logit(1.0))
    assert _logit(0.99999999) > 10  # saturated -> large positive logit
    assert _logit(0.00000001) < -10  # opposite saturation


# ---- saturation thresholds ----


def test_saturation_thresholds_are_documented():
    from scripts.cipro_bounded_falsifier import (
        SATURATION_PROBA_THRESHOLD,
        SATURATION_MAX_ABS_DELTA_THRESHOLD,
    )
    # Pinned per plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md
    # Changing these values requires an updated plan + commit message rationale.
    assert SATURATION_PROBA_THRESHOLD == 0.95
    assert SATURATION_MAX_ABS_DELTA_THRESHOLD == 0.01


# ---- verdict matrix on synthetic StrainResults ----


def _mk_result(
    bucket: str,
    best_pos: int | None,
    best_abs: int | None,
    saturated: bool = False,
    indeterminate: str | None = None,
    hits: list[dict] | None = None,
):
    from scripts.cipro_bounded_falsifier import StrainResult
    r = StrainResult(
        strain_id="x",
        accession="y",
        label="R",
        bucket=bucket,
        n_cached_genes=5000,
        baseline_proba_R=0.99 if saturated else 0.6,
        baseline_logit=4.0 if saturated else 0.4,
        max_abs_delta_all_genes=0.001 if saturated else 0.1,
        saturation_flag=saturated,
        n_known_loci_hits=len(hits) if hits else 0,
        best_known_locus_rank_pos_delta=best_pos,
        best_known_locus_rank_abs_delta=best_abs,
        indeterminate_reason=indeterminate,
    )
    if hits:
        r.per_known_locus = hits
    return r


def test_compute_verdict_pass_when_all_buckets_pass():
    from scripts.cipro_bounded_falsifier import compute_verdict
    # Bucket A: 4 ERS with QRDR in top-10 positive ranking
    A = [_mk_result("A_ERS", 5, 12, hits=[{"alias": "gyrA", "rank_pos_delta": 5, "rank_abs_delta": 12}]) for _ in range(4)]
    # Bucket B: all 4 in top-10, median rank ratio 1727 / 8 ≈ 216x
    B = [_mk_result("B_ELX", 8, 1727, hits=[{"alias": "gyrA", "rank_pos_delta": 8, "rank_abs_delta": 1727}]) for _ in range(4)]
    # Bucket C: indeterminate-with-saturation
    C = [_mk_result("C_NEGATIVE", None, None, saturated=True, indeterminate="ALL_NEGATIVE_DELTA") for _ in range(4)]
    v, d = compute_verdict(A, B, C)
    assert v == "PASS"
    assert d["bucket_A_pass"]
    assert d["bucket_B_pass"]
    assert d["bucket_C_handled"]


def test_compute_verdict_fail_when_bucket_B_doesnt_improve():
    from scripts.cipro_bounded_falsifier import compute_verdict
    A = [_mk_result("A_ERS", 5, 12, hits=[{"alias": "gyrA", "rank_pos_delta": 5, "rank_abs_delta": 12}]) for _ in range(4)]
    B = [_mk_result("B_ELX", None, 1727) for _ in range(4)]
    C = [_mk_result("C_NEGATIVE", None, None, saturated=True, indeterminate="ALL_NEGATIVE_DELTA") for _ in range(4)]
    v, d = compute_verdict(A, B, C)
    assert v == "FAIL"
    assert d["bucket_A_pass"]
    assert not d["bucket_B_pass"]


def test_compute_verdict_runner_regression_when_bucket_A_fails():
    from scripts.cipro_bounded_falsifier import compute_verdict
    A = [_mk_result("A_ERS", None, 12) for _ in range(4)]
    B = [_mk_result("B_ELX", 8, 1727, hits=[{"alias": "gyrA", "rank_pos_delta": 8, "rank_abs_delta": 1727}]) for _ in range(4)]
    C = [_mk_result("C_NEGATIVE", None, None, saturated=True, indeterminate="ALL_NEGATIVE_DELTA") for _ in range(4)]
    v, d = compute_verdict(A, B, C)
    assert v == "RUNNER_REGRESSION"
    assert not d["bucket_A_pass"]


def test_compute_verdict_revert_when_bucket_C_unhandled():
    """A+B pass but Bucket C is neither recovered nor all-indeterminate-saturated."""
    from scripts.cipro_bounded_falsifier import compute_verdict
    A = [_mk_result("A_ERS", 5, 12, hits=[{"alias": "gyrA", "rank_pos_delta": 5, "rank_abs_delta": 12}]) for _ in range(4)]
    B = [_mk_result("B_ELX", 8, 1727, hits=[{"alias": "gyrA", "rank_pos_delta": 8, "rank_abs_delta": 1727}]) for _ in range(4)]
    # Bucket C: ambiguous — neither recovered (best_pos=None) nor saturated
    C = [_mk_result("C_NEGATIVE", None, None, saturated=False, indeterminate="ALL_NEGATIVE_DELTA") for _ in range(4)]
    v, d = compute_verdict(A, B, C)
    assert v == "REVERT"


# ---- bucket_A_pass / bucket_B_pass / bucket_C_handled ----


def test_bucket_A_pass_requires_at_least_3_qrdr_in_top_10():
    from scripts.cipro_bounded_falsifier import bucket_A_pass
    # 3 of 4 with QRDR in top-10 -> pass
    results = [
        _mk_result("A_ERS", 5, 10, hits=[{"alias": "gyrA", "rank_pos_delta": 5, "rank_abs_delta": 10}]),
        _mk_result("A_ERS", 7, 12, hits=[{"alias": "parC", "rank_pos_delta": 7, "rank_abs_delta": 12}]),
        _mk_result("A_ERS", 9, 14, hits=[{"alias": "gyrB", "rank_pos_delta": 9, "rank_abs_delta": 14}]),
        _mk_result("A_ERS", 200, 300, hits=[{"alias": "gyrA", "rank_pos_delta": 200, "rank_abs_delta": 300}]),
    ]
    assert bucket_A_pass(results) is True
    # Only 2 in top-10 -> fail
    results_fail = results.copy()
    results_fail[2] = _mk_result("A_ERS", 50, 60, hits=[{"alias": "gyrA", "rank_pos_delta": 50, "rank_abs_delta": 60}])
    assert bucket_A_pass(results_fail) is False


def test_bucket_B_pass_requires_2_in_top_10_and_100x_rank_shift():
    from scripts.cipro_bounded_falsifier import bucket_B_pass
    # 2 of 4 in top-10 + median ratio 1727/8 ≈ 216x (passes 100x bar)
    results = [
        _mk_result("B_ELX", 5, 1727),
        _mk_result("B_ELX", 8, 3178),
        _mk_result("B_ELX", 50, 2644),
        _mk_result("B_ELX", 200, 1110),
    ]
    assert bucket_B_pass(results) is True
    # Only 1 in top-10 -> fail
    results_fail = [
        _mk_result("B_ELX", 5, 1727),
        _mk_result("B_ELX", 200, 3178),
        _mk_result("B_ELX", 200, 2644),
        _mk_result("B_ELX", 200, 1110),
    ]
    assert bucket_B_pass(results_fail) is False


def test_bucket_C_handled_when_all_indeterminate_with_saturation():
    from scripts.cipro_bounded_falsifier import bucket_C_handled
    results = [
        _mk_result("C_NEGATIVE", None, None, saturated=True, indeterminate="ALL_NEGATIVE_DELTA")
        for _ in range(4)
    ]
    handled, desc = bucket_C_handled(results)
    assert handled is True


def test_bucket_C_handled_when_some_recovered():
    from scripts.cipro_bounded_falsifier import bucket_C_handled
    results = [
        _mk_result("C_NEGATIVE", 45, 100, hits=[{"alias": "qnr", "rank_pos_delta": 45}]),
        _mk_result("C_NEGATIVE", None, None, indeterminate="ALL_NEGATIVE_DELTA"),
        _mk_result("C_NEGATIVE", None, None, indeterminate="ALL_NEGATIVE_DELTA"),
        _mk_result("C_NEGATIVE", None, None, indeterminate="ALL_NEGATIVE_DELTA"),
    ]
    handled, desc = bucket_C_handled(results)
    assert handled is True  # 1 recovered into top-50


def test_bucket_C_unhandled_when_ambiguous():
    """Bucket C: neither recovered nor saturated -> not handled."""
    from scripts.cipro_bounded_falsifier import bucket_C_handled
    results = [
        _mk_result("C_NEGATIVE", None, None, saturated=False, indeterminate="ALL_NEGATIVE_DELTA")
        for _ in range(4)
    ]
    handled, desc = bucket_C_handled(results)
    assert handled is False
