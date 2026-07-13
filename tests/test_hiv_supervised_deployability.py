"""Tests for the deployability stress-test helpers (offline; real-data run skip-guarded)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))
pytest.importorskip("sklearn")

spec = importlib.util.spec_from_file_location(
    "hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_blindspot_metrics_single_class_guard():
    # a blind-spot test fold with only S (no R) must NOT crash or fake an AUROC
    y = [0, 0, 0, 1]; scores = [0.1, 0.2, 0.3, 0.9]; burden = [1, 2, 1, 3]
    neg_all_S = [True, True, False, False]      # the two catalog-negative isolates are both S
    m = mod._blindspot_metrics(y, scores, neg_all_S, burden)
    assert m["auroc"] is None and m["pass"] is False and "degenerate" in m["note"]


def test_blindspot_metrics_scores_directionally():
    y = [0, 0, 1, 1, 0, 1]; scores = [0.1, 0.2, 0.8, 0.9, 0.15, 0.85]; burden = [1, 1, 1, 1, 1, 1]
    neg = [True] * 6
    m = mod._blindspot_metrics(y, scores, neg, burden)
    assert m["auroc"] == pytest.approx(1.0) and m["n"] == 6 and m["R"] == 3


def test_config_matches_prior_experiment():
    assert mod.DRUG == "EFV" and mod.CUTOFF == 3.0
