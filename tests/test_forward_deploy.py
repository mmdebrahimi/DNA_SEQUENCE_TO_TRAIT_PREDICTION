"""Offline tests for the one-call novel-protein orchestrator (heavy scorers monkeypatched)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import predict_hybrid_from_sequence  # noqa: E402
import dna_decode.forward.deploy as dep  # noqa: E402


def test_uses_precomputed_tables_without_touching_heavy_scorers(monkeypatch):
    # both tables supplied -> no esm_scorer / prosst_scorer import at all
    def boom(*a, **k):
        raise AssertionError("heavy scorer called despite precomputed tables")
    monkeypatch.setattr("dna_decode.forward.esm_scorer.esm2_logp_table", boom, raising=False)
    esm = {"W1L": 2.0, "W1P": -4.0, "W1V": 0.0}
    prosst = {"W1L": 1.5, "W1P": -3.0, "W1V": 0.1}
    p = predict_hybrid_from_sequence("W", "W1P", esm_table=esm, prosst_table=prosst)
    assert p.predicted_effect == "damaging" and p.method == "hybrid_esm2_prosst"


def test_computes_esm_when_absent(monkeypatch):
    # esm_table absent -> esm2_logp_table is invoked (stubbed); prosst supplied
    called = {"esm": False}
    def fake_esm(seq, model_name=None):
        called["esm"] = True
        return {1: {"W": 0.0, "L": 1.0, "P": -4.0, "V": -0.5}}
    monkeypatch.setattr("dna_decode.forward.esm_scorer.esm2_logp_table", fake_esm)
    prosst = {"W1L": 1.0, "W1P": -3.0, "W1V": 0.0}
    p = predict_hybrid_from_sequence("W", "W1L", prosst_table=prosst)
    assert called["esm"] and 0.0 <= p.raw_score <= 1.0


def test_needs_a_prosst_source():
    with pytest.raises(ValueError, match="need one of prosst_table"):
        predict_hybrid_from_sequence("W", "W1L", esm_table={"W1L": 1.0, "W1P": -1.0})


def test_all_single_mutants_saturation():
    m = dep._all_single_mutants("MA")
    assert len(m) == 2 * 19 and "M1L" in m and "A2G" in m and "M1M" not in m
