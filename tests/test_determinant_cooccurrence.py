"""Tests for the determinant co-occurrence / linkage world model (`scripts/determinant_cooccurrence.py`,
Family C). Pure-function tests run offline on synthetic data (numpy + sklearn are CORE -> CI-safe). The
end-to-end real-data test SKIPS when the cached AMRFinder runs (gitignored) are absent.
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import determinant_cooccurrence as dc  # noqa: E402


def test_organism_of_strips_drug_and_qualifier_suffixes():
    assert dc.organism_of("/x/data/raw/klebsiella_provdisjoint_ceftriaxone/amrfinder_runs/GCA_1/main.tsv") == "klebsiella"
    assert dc.organism_of("/x/data/raw/escherichia_coli_shigella_ciprofloxacin/amrfinder_runs/G/main.tsv") == "escherichia_coli_shigella"
    assert dc.organism_of("/x/data/raw/campylobacter_indep_ciprofloxacin/amrfinder_runs/G/main.tsv") == "campylobacter"


def test_auc_perfect_random_anti():
    lab = np.array([True] * 5 + [False] * 5)
    assert dc._auc(np.array([9, 8, 7, 6, 5, 4, 3, 2, 1, 0.0]), lab) == pytest.approx(1.0)   # perfect
    assert dc._auc(np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9.0]), lab) == pytest.approx(0.0)   # anti
    assert dc._auc(np.array([5.0] * 10), lab) == pytest.approx(0.5)                          # ties -> 0.5


def test_matrix_and_testable_both_class_filter():
    # det 'A' present in 12 of 30; 'B' present in 2 (below MIN_SUPPORT); 'C' present in 29 (absent<support)
    acc_org, acc_dets = {}, {}
    for i in range(30):
        acc_org[f"g{i}"] = "org"
        s = set()
        if i < 12: s.add("A")
        if i < 2: s.add("B")
        if i < 29: s.add("C")
        acc_dets[f"g{i}"] = s
    X, dets, didx, accs, testable = dc.matrix_for_organism(acc_org, acc_dets, "org", min_support=10)
    assert set(dets) == {"A", "B", "C"}
    assert "A" in testable                    # 12 present, 18 absent -> both-class
    assert "B" not in testable and "C" not in testable   # below both-class support


def test_dedup_profiles_keeps_one_per_unique_profile():
    X = np.array([[1, 0], [1, 0], [0, 1], [1, 1]], float)   # rows 0,1 identical
    keep = dc.dedup_profiles(X, ["a", "b", "c", "d"])
    assert len(keep) == 3                       # 3 unique profiles
    assert 0 in keep and 1 not in keep


def test_lift_table_computes_conditional_lift():
    # 20 genomes: target T present in 10; det X present in exactly those 10 (perfect co-occurrence) -> high lift
    n = 20
    X = np.zeros((n, 2))
    X[:10, 0] = 1.0   # T
    X[:10, 1] = 1.0   # X co-occurs perfectly
    lt = dc.lift_table(X, ["T", "X"], {"T": 0, "X": 1}, ["T"], min_cooc=5)
    row = lt["T"][0]
    assert row["det"] == "X"
    assert row["p_given_target"] == pytest.approx(1.0)      # P(X|T)=1
    assert row["lift"] == pytest.approx(2.0)                # P(X)=0.5 -> lift 1/0.5 = 2


def test_linkage_detects_perfect_linkage_and_independence():
    rng = np.random.default_rng(0)
    n = 120
    d1 = (rng.random(n) < 0.4).astype(float)
    X = np.column_stack([d1, d1.copy(), (rng.random(n) < 0.4).astype(float)])  # col0==col1; col2 independent
    linked = dc.linkage_for_determinant(X, 0, seed=0, boot=200)   # predict col0 from col1(==) + col2
    assert linked["linked"] is True and linked["auc"] > 0.9       # perfectly predictable
    indep = dc.linkage_for_determinant(X, 2, seed=0, boot=200)    # col2 independent of others
    assert indep["auc"] < 0.75                                    # near chance (not strongly linked)


def test_drug_inversion_groups_by_class_and_ranks():
    acc_org = {"g0": "o", "g1": "o", "g2": "o"}
    acc_dets = {"g0": {"blaX", "sul1"}, "g1": {"blaX"}, "g2": {"sul1"}}
    det_class = {"blaX": "CEPHALOSPORIN", "sul1": "SULFONAMIDE"}
    inv = dc.drug_inversion(acc_org, acc_dets, det_class, "o")
    assert inv["CEPHALOSPORIN"][0] == {"det": "blaX", "n": 2}     # blaX in 2 genomes, ranked first
    assert inv["SULFONAMIDE"][0]["det"] == "sul1"


def test_prereg_constants_frozen():
    assert dc.PASS_FRACTION == 0.5 and dc.LINKED_AUC == 0.5
    assert dc.MIN_GENOMES == 60 and dc.MIN_SUPPORT == 10


def test_module_does_not_import_frozen_amr_surface():
    src = (REPO / "scripts" / "determinant_cooccurrence.py").read_text(encoding="utf-8")
    assert "amr_rules" not in src and "calibrated_amr_rules" not in src and "shipped_decoder_surface" not in src


# --- end-to-end on real cached AMRFinder runs (skips when absent) ---
_HAS_RUNS = bool(glob.glob(str(REPO / "data" / "raw" / "*" / "amrfinder_runs" / "*" / "main.tsv")))


@pytest.mark.skipif(not _HAS_RUNS, reason="cached AMRFinder runs absent (gitignored)")
def test_harvest_and_run_organism_real_smoke():
    acc_org, acc_dets, det_class = dc.harvest(REPO)
    assert len(acc_org) > 100                        # many genomes cached
    from collections import Counter
    big = Counter(acc_org.values()).most_common(1)[0][0]
    res = dc.run_organism(acc_org, acc_dets, det_class, big, boot=100)
    assert res["n_testable"] >= 1
    assert 0.0 <= res["fraction_linked"] <= 1.0
    assert res["lift_table"] and res["phenotype_to_genotype_inversion"]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
