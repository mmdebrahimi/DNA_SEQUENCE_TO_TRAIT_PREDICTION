"""Tests for the TB supervised-vs-catalog bacterial transfer test (offline pure helpers; real run skip-guarded)."""
from __future__ import annotations

import importlib.util, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))
pytest.importorskip("sklearn")

spec = importlib.util.spec_from_file_location("tb_supervised_vs_catalog", REPO / "scripts" / "tb_supervised_vs_catalog.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


def test_auroc_directional():
    assert mod._auroc([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]) == pytest.approx(1.0)
    assert mod._auroc([0, 0, 1, 1], [0.9, 0.8, 0.2, 0.1]) == pytest.approx(0.0)
    assert mod._auroc([1, 1, 1], [0.1, 0.2, 0.3]) is None      # single-class


def test_blind_single_class_guard():
    m = mod._blind([0, 0, 1], [0.1, 0.2, 0.9], [True, True, False], [1, 1, 1])
    assert m["auroc"] is None and m["pass"] is False and "single-class" in m["note"]


def test_blind_pass_requires_beat_burden_and_null():
    n = 40
    y = [0] * (n // 2) + [1] * (n // 2)
    sc = [0.1 + 0.001 * i for i in range(n // 2)] + [0.8 + 0.001 * i for i in range(n // 2)]
    neg = [True] * n
    burden = [1] * n                                   # constant burden -> burden AUROC 0.5
    m = mod._blind(y, sc, neg, burden)
    assert m["auroc"] == pytest.approx(1.0) and m["null"] < 0.55 and m["pass"] is True


def test_rpob_span_and_drug():
    assert mod.DRUG == "rifampicin" and mod.CODE == "RIF"
    assert mod.RPOB.start == 759807 and mod.RPOB.stop == 763326   # H37Rv rpoB span
