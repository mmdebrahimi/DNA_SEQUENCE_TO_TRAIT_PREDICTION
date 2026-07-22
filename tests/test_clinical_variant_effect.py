"""Offline tests for the clinical-significance validation cell (scripts/clinical_variant_effect_validate.py).

Pure helpers only — no MaveDB/ClinVar network. Pins: AUROC correctness (incl. ties + perfect/inverse
separation), label-free DMS orientation, single-class AUROC-inapplicable gating, and the join logic.
"""
from __future__ import annotations

import math

import pytest

from scripts.clinical_variant_effect_validate import (
    auroc, _spearman_sign, score_from_tables, MIN_PER_CLASS,
)

META = {"uniprot": "P00000", "urn": "urn:mavedb:test", "assay": "test assay"}


def test_auroc_perfect_separation():
    # positives all score higher than negatives -> AUROC 1.0
    labels = [True, True, False, False]
    scores = [0.9, 0.8, 0.2, 0.1]
    assert auroc(labels, scores) == 1.0


def test_auroc_inverse_separation():
    labels = [True, True, False, False]
    scores = [0.1, 0.2, 0.8, 0.9]
    assert auroc(labels, scores) == 0.0


def test_auroc_chance_and_ties():
    # interleaved -> 0.5; all-tied -> 0.5 (mid-rank tie correction)
    assert auroc([True, False, True, False], [1.0, 1.0, 1.0, 1.0]) == 0.5
    # one tie straddling the boundary
    a = auroc([True, True, False, False], [0.5, 0.5, 0.5, 0.1])
    assert 0.5 <= a <= 1.0


def test_auroc_single_class_raises():
    with pytest.raises(ValueError):
        auroc([True, True, True], [1.0, 2.0, 3.0])


def test_spearman_sign_labelfree_orientation():
    assert _spearman_sign([1, 2, 3, 4], [10, 20, 30, 40]) == 1.0
    assert _spearman_sign([1, 2, 3, 4], [40, 30, 20, 10]) == -1.0
    assert _spearman_sign([1, 1, 1], [5, 6, 7]) == 0.0


def test_score_single_class_is_auroc_inapplicable():
    # 20 pathogenic, 1 benign -> below MIN_PER_CLASS on benign -> a FINDING, not a fake number
    dms = {("A", i, "V"): float(i) for i in range(1, 22)}
    clin = {("A", i, "V"): "PATH" for i in range(1, 21)}
    clin[("A", 21, "V")] = "BENIGN"
    rec = score_from_tables("GENEX", META, dms, clin, set())
    assert rec["auroc_applicable"] is False
    assert "dms_auroc" not in rec
    assert rec["n_path"] == 20 and rec["n_benign"] == 1
    assert "single-class" in rec["note"]


def test_score_orientation_recovers_high_auroc_when_dms_anticorrelates_blosum():
    # Build a joined set where LOW DMS = pathogenic, and DMS is oriented so higher=MORE damaging (anti-BLOSUM).
    # The label-free BLOSUM-anchored flip must still recover a high pathogenic-AUROC.
    dms, clin = {}, {}
    # pathogenic: radical substitutions (low BLOSUM), assigned HIGH raw dms (anti-oriented)
    path_muts = [("W", i, "R") for i in range(1, 21)]      # W->R radical
    ben_muts = [("I", i, "V") for i in range(21, 41)]       # I->V conservative (high BLOSUM)
    for k in path_muts:
        dms[k] = 100.0 + k[1]     # high raw score
        clin[k] = "PATH"
    for k in ben_muts:
        dms[k] = 1.0 + k[1]        # low raw score
        clin[k] = "BENIGN"
    rec = score_from_tables("GENEY", META, dms, clin, {"GENEY"})
    assert rec["auroc_applicable"] is True
    assert rec["in_proteingym"] is True
    # radical=pathogenic, conservative=benign -> BLOSUM alone separates well
    assert rec["blosum_auroc"] >= 0.9
    # DMS raw is anti-oriented (high=pathogenic) but BLOSUM-anchored flip fixes it -> high pathogenic-AUROC
    assert rec["dms_auroc"] >= 0.9
    assert rec["dms_orientation"] == "flipped-to-agree-with-blosum"


def test_score_join_is_intersection_only():
    dms = {("A", 1, "V"): 0.1, ("A", 2, "V"): 0.2, ("A", 3, "V"): 0.3}
    clin = {("A", 2, "V"): "PATH", ("A", 99, "L"): "BENIGN"}   # only (A,2,V) is shared
    rec = score_from_tables("GENEZ", META, dms, clin, set())
    assert rec["n_joined"] == 1


def test_min_per_class_constant_sane():
    assert MIN_PER_CLASS >= 10
