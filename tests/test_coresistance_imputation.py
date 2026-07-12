"""Tests for co-resistance imputation (`scripts/coresistance_imputation.py`, Family C deep-dive).
Pure-function tests offline (numpy+sklearn core -> CI-safe); real-data smoke SKIPS when cached AMRFinder
runs are absent.
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import coresistance_imputation as ci  # noqa: E402


def test_genome_classes_splits_compound_classes():
    acc_dets = {"g": {"blaX", "qnrB", "aacY"}}
    det_class = {"blaX": "BETA-LACTAM", "qnrB": "QUINOLONE", "aacY": "AMINOGLYCOSIDE/QUINOLONE"}
    cls = ci.genome_classes(acc_dets, det_class, "g")
    assert cls == {"BETA-LACTAM", "QUINOLONE", "AMINOGLYCOSIDE"}   # compound split on '/'


def test_class_matrix_and_testable():
    # 30 genomes: class A in 12 (both-class), B in 2 (below support), C in 30 (no absent)
    acc_org, acc_dets, det_class = {}, {}, {"dA": "A", "dB": "B", "dC": "C"}
    for i in range(30):
        acc_org[f"g{i}"] = "o"
        s = set()
        if i < 12: s.add("dA")
        if i < 2: s.add("dB")
        s.add("dC")
        acc_dets[f"g{i}"] = s
    X, classes, cidx, accs, testable = ci.class_matrix_for_organism(acc_org, acc_dets, det_class, "o")
    assert set(classes) == {"A", "B", "C"}
    assert testable == ["A"]        # only A has both-class support >= MIN_SUPPORT(10)


def test_run_organism_imputable_when_classes_linked():
    # class A present iff class B present (perfect co-resistance) -> A imputable from B
    rng = np.random.default_rng(0)
    acc_org, acc_dets, det_class = {}, {}, {"dA": "A", "dB": "B", "dC": "C"}
    for i in range(120):
        acc_org[f"g{i}"] = "o"
        has = rng.random() < 0.4
        s = set()
        if has: s.update({"dA", "dB"})           # A and B always together
        if rng.random() < 0.4: s.add("dC")        # C independent
        acc_dets[f"g{i}"] = s
    r = ci.run_organism(acc_org, acc_dets, det_class, "o", boot=200)
    assert r["per_class_imputation"]["A"]["imputable"] is True
    assert r["per_class_imputation"]["A"]["auc"] > 0.9


def test_prereg_constants():
    assert ci.PASS_FRACTION == 0.5 and ci.MIN_GENOMES == 60


def test_module_does_not_import_frozen_amr_surface():
    src = (REPO / "scripts" / "coresistance_imputation.py").read_text(encoding="utf-8")
    assert "amr_rules" not in src and "calibrated_amr_rules" not in src and "shipped_decoder_surface" not in src


_HAS_RUNS = bool(glob.glob(str(REPO / "data" / "raw" / "*" / "amrfinder_runs" / "*" / "main.tsv")))


@pytest.mark.skipif(not _HAS_RUNS, reason="cached AMRFinder runs absent (gitignored)")
def test_run_all_real_smoke():
    res = ci.run_all(REPO, boot=100)
    assert res["total_testable"] >= 1
    assert res["verdict"] in ("PASS_CORESISTANCE_IMPUTABLE", "FAIL_CLASSES_INDEPENDENT")
    # at least one organism should have a co-resistance network entry
    assert any(r["coresistance_network"] for r in res["per_organism"].values())


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
