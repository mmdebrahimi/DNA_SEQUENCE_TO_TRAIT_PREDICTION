"""Tests for the ProteinGym blosum62 rank-inverse generalization sweep (N=217).

Pins the coordinate gate + the reuse of the censored-assay/underpowered gates. Real-data tests skip
without the D: ProteinGym cache.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.forward_inverse_proteingym_sweep import (
    DMS_DIR,
    MATERIAL_MARGIN,
    REF,
    build_candidates,
    score_one,
)
from scripts.forward_inverse_deployable import TARGET_PCTS, random_null


def _write_dms(tmp_path, rows):
    p = tmp_path / "assay.csv"
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["mutant", "DMS_score"])
        w.writeheader()
        for m, s in rows:
            w.writerow({"mutant": m, "DMS_score": s})
    return p


def test_material_margin_matches_the_4_protein_sweep_not_retuned():
    assert MATERIAL_MARGIN == 0.25


# ---- the coordinate gate -------------------------------------------------------------------------------

def test_build_candidates_drops_a_wt_mismatch_rather_than_coercing(tmp_path):
    target = "MK"                                   # pos1=M, pos2=K
    dms = _write_dms(tmp_path, [("M1L", -1.0), ("A2L", -2.0)])   # A2L: pos2 is K, not A
    cands, dropped = build_candidates(target, dms)
    assert {c.mutant for c in cands} == {"M1L"}
    assert dropped == 1


def test_build_candidates_skips_multi_mutants(tmp_path):
    dms = _write_dms(tmp_path, [("M1L:K2R", -1.0), ("M1L", -1.0)])
    cands, _ = build_candidates("MK", dms)
    assert {c.mutant for c in cands} == {"M1L"}


def test_build_candidates_skips_out_of_range_and_nonstandard(tmp_path):
    dms = _write_dms(tmp_path, [("M9L", -1.0), ("M1X", -1.0), ("M1L", -1.0)])
    cands, dropped = build_candidates("MK", dms)
    assert {c.mutant for c in cands} == {"M1L"}      # M9L out of range, M1X nonstandard alt
    assert dropped >= 1


def test_build_candidates_requires_a_measured_value(tmp_path):
    dms = _write_dms(tmp_path, [("M1L", ""), ("M1V", -1.0)])
    cands, _ = build_candidates("MK", dms)
    assert {c.mutant for c in cands} == {"M1V"}


# ---- score_one gating ----------------------------------------------------------------------------------

def _null():
    return random_null(TARGET_PCTS, 5)


def test_score_one_flags_a_censored_assay_rather_than_scoring_it(tmp_path):
    # 300 variants across 25 positions, but 90% tied at one value -> censored
    rows = []
    for i in range(1, 26):
        for alt in "ACDEFGHIKLM":
            rows.append((f"A{i}{alt}", -2.0 if (i + ord(alt)) % 10 else -5.0))
    # force heavy ties: make most values identical
    rows = [(m, -2.0) for m, _ in rows[:270]] + [(m, v) for m, v in
            [(f"A{i}N", -float(i % 8)) for i in range(1, 26)]]
    target = "A" * 26
    dms = _write_dms(tmp_path, rows)
    row = {"DMS_id": "x", "taxon": "t", "source_organism": "o", "seq_len": "26", "target_seq": target}
    # point DMS_DIR resolution by writing the file where score_one looks
    import scripts.forward_inverse_proteingym_sweep as mod
    orig = mod.DMS_DIR
    mod.DMS_DIR = tmp_path
    try:
        (tmp_path / "x.csv").write_text(dms.read_text(encoding="utf-8"), encoding="utf-8")
        res = mod.score_one(row, _null())
    finally:
        mod.DMS_DIR = orig
    assert res["status"] == "DEGENERATE_CENSORED_ASSAY"


def test_score_one_underpowered_when_too_few_candidates(tmp_path):
    target = "MKAA"
    dms = _write_dms(tmp_path, [("M1L", -1.0), ("K2R", -2.0)])
    row = {"DMS_id": "x", "taxon": "", "source_organism": "", "seq_len": "4", "target_seq": target}
    import scripts.forward_inverse_proteingym_sweep as mod
    orig = mod.DMS_DIR
    mod.DMS_DIR = tmp_path
    try:
        (tmp_path / "x.csv").write_text(dms.read_text(encoding="utf-8"), encoding="utf-8")
        res = mod.score_one(row, _null())
    finally:
        mod.DMS_DIR = orig
    assert res["status"] == "UNDERPOWERED"


# ---- real data (skips without D:) ----------------------------------------------------------------------

@pytest.mark.skipif(not REF.exists(), reason="ProteinGym cache (D:) not mounted")
def test_reference_has_the_expected_217_assays():
    with REF.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 217
    assert all("target_seq" in r and "DMS_id" in r for r in rows)


@pytest.mark.skipif(not REF.exists() or not DMS_DIR.exists(), reason="ProteinGym cache (D:) not mounted")
def test_a_real_assay_scores_and_is_graded_in_percentile_points():
    with REF.open(encoding="utf-8") as fh:
        row = next(r for r in csv.DictReader(fh) if (DMS_DIR / f"{r['DMS_id']}.csv").exists())
    res = score_one(row, _null())
    assert res["status"] in ("SCORED", "DEGENERATE_CENSORED_ASSAY", "UNDERPOWERED")
    if res["status"] == "SCORED":
        assert 0.0 <= res["blosum_pct_err_best_of_k"] <= 1.0
        assert 0.0 <= res["null_pct_err_best_of_k"] <= 1.0
