"""Tests for the expression-gated GPR runner (`scripts/fba_expression_gated_gpr.py`).

The run consumes real PRECISE-1K data off `D:`, so these cover the PURE decision logic — the places a
silent bug would have turned a negative result into a fake positive:

  * an UNMEASURED gene must never be gated off (absence of evidence != evidence of absence);
  * a gene absent from the model must never be gated;
  * the percentile must be computed over the MEASURED-AND-IN-MODEL set, not the whole matrix;
  * the frozen target set must refuse to run if it ever changes size;
  * the confusion counts that feed the guardrail must be right, since the guardrail is what
    actually decided this experiment.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.fba_expression_gated_gpr import confusion, frozen_target_set, gated_genes  # noqa: E402


# ------------------------------------------------------------------------------- the gate

def test_gates_only_genes_below_the_percentile():
    expr = {f"g{i}": float(i) for i in range(1, 11)}       # 1..10
    model = set(expr)
    off = gated_genes(expr, model, 30.0)                    # 30th pctl of 1..10 == 3.7
    assert off == {"g1", "g2", "g3"}


def test_an_unmeasured_gene_is_never_gated():
    """Load-bearing: mirrors the bridge's `eval_gpr` returning None for an unmeasured gene. Gating an
    unmeasured gene would knock it out on no evidence and manufacture essentiality."""
    expr = {"g1": 1.0, "g2": 100.0}
    model = {"g1", "g2", "g_unmeasured"}
    assert "g_unmeasured" not in gated_genes(expr, model, 90.0)


def test_a_gene_absent_from_the_model_is_never_gated():
    expr = {"g1": 1.0, "not_in_model": 0.0}
    assert gated_genes(expr, {"g1"}, 99.0) <= {"g1"}


def test_the_percentile_is_computed_over_the_in_model_set_only():
    """If the cutoff were taken over the whole expression matrix, out-of-model genes would drag it and
    silently change how many in-model genes get gated."""
    expr = {"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0, **{f"x{i}": 0.0 for i in range(100)}}
    model = {"a", "b", "c", "d"}
    off = gated_genes(expr, model, 50.0)                    # median of 1,2,3,4 == 2.5
    assert off == {"a", "b"}


def test_empty_measured_set_gates_nothing():
    assert gated_genes({"z": 1.0}, {"other"}, 50.0) == set()


# --------------------------------------------------------------------------- the frozen target set

def test_frozen_target_set_is_the_committed_eight():
    genes = frozen_target_set()
    assert len(genes) == 8
    assert genes == sorted(genes)
    src = json.loads(Path("wiki/fba_orphan_protection_2026-08-21.json").read_text(encoding="utf-8"))
    assert set(genes) == set(src["impact_on_experimental_deficit"]["genes_isozyme_masked"])


def test_frozen_target_set_refuses_a_changed_artifact(tmp_path, monkeypatch):
    """The pre-registration froze EIGHT genes. If the upstream artifact is ever regenerated with a
    different set, the run must abort rather than silently score a different endpoint."""
    import scripts.fba_expression_gated_gpr as mod

    p = tmp_path / "changed.json"
    p.write_text(json.dumps(
        {"impact_on_experimental_deficit": {"genes_isozyme_masked": ["b1", "b2"]}}), encoding="utf-8")
    monkeypatch.setattr(mod, "SCREEN_ARTIFACT", p)
    with pytest.raises(SystemExit):
        mod.frozen_target_set()


# ------------------------------------------------------------------------------ confusion counting

def test_confusion_counts_every_quadrant():
    genes, conds = ["g1", "g2"], ["c1"]
    calls = {"c1": {"g1": True, "g2": True}}
    truth = {("g1", "c1"): True, ("g2", "c1"): False}
    assert confusion(calls, truth, genes, conds) == (1, 1, 0, 0)


def test_a_missing_call_counts_as_not_essential_not_as_skipped():
    """The grid is fixed, so an absent prediction is a NEGATIVE call and must still land in the
    denominator -- silently dropping it would flatter the false-positive guardrail."""
    genes, conds = ["g1"], ["c1"]
    assert confusion({}, {("g1", "c1"): True}, genes, conds) == (0, 0, 1, 0)
    assert confusion({}, {}, genes, conds) == (0, 0, 0, 1)
