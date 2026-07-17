"""Tests for the FLC-route test — the flowering cell's distinctive claim."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.flowering_flc_route_test import FLC_CSV, TABLE_S3, _confusion, _predict, score

ART = Path(__file__).resolve().parent.parent / "wiki" / "flowering_flc_route_2026-07-17.json"


# ---- the FRI-only collapse (why the distinctive claim went untested) --------------------------------------

def test_unobserved_flc_collapses_the_cell_to_the_naive_fri_only_rule():
    """This IS what the Table S3 run had to do, and why it could not test the FLC route: with FLC
    unobserved the two-locus AND degenerates to 'late iff functional FRI'."""
    assert _predict(fri_lof=False, flc_strong=None) == "late"      # FRI-only says LATE for every one
    assert _predict(fri_lof=True, flc_strong=None) == "early"


def test_the_flc_route_changes_the_call_where_fri_only_cannot():
    """The Da(1)-12 class: functional FRI + weak FLC. A FRI-only rule says LATE; the cell says EARLY."""
    assert _predict(fri_lof=False, flc_strong=None) == "late"       # naive
    assert _predict(fri_lof=False, flc_strong=False) == "early"     # the cell's distinctive call


def test_flc_is_downstream_so_weak_flc_calls_early_regardless_of_fri():
    for lof in (True, False):
        assert _predict(fri_lof=lof, flc_strong=False) == "early"


def test_only_functional_fri_and_strong_flc_calls_late():
    assert _predict(False, True) == "late"
    assert _predict(True, True) == "early"


# ---- the harness -----------------------------------------------------------------------------------------

def test_score_reports_both_rules_on_the_same_rows():
    rows = [{"fri_lof": False, "flc_expr": 2.0, "ft16": 90.0},
            {"fri_lof": False, "flc_expr": 0.1, "ft16": 40.0},
            {"fri_lof": True, "flc_expr": 0.1, "ft16": 40.0}]
    s = score(rows, flc_cut=1.0, ft_cut=60.0)
    assert s["two_locus"]["accuracy"] == 1.0          # FLC rescues row 2
    assert s["fri_only"]["accuracy"] == pytest.approx(2 / 3)   # ...which FRI-only mis-calls LATE


def test_confusion_flags_a_single_class_set_degenerate():
    assert _confusion([("late", "late"), ("early", "late")])["degenerate"] is True


def test_confusion_reports_the_constant_predictor_null():
    m = _confusion([("late", "late")] * 8 + [("late", "early")] * 2)
    assert m["null_accuracy"] == 0.8


# ---- real data -------------------------------------------------------------------------------------------

@pytest.mark.skipif(not ART.exists(), reason="FLC-route artifact not present")
def test_all_four_rule_cells_are_called_correctly_on_real_data():
    """The cell's whole shape: each of the AND's 4 cells must call its majority correctly."""
    r = json.loads(ART.read_text(encoding="utf-8"))
    assert r["verdict"] == "FLC_ROUTE_VALIDATED_ADDS_OVER_FRI_ONLY"
    for c in r["rule_cells"]:
        assert c["call_correct_for_majority"], c


@pytest.mark.skipif(not ART.exists(), reason="FLC-route artifact not present")
def test_the_da112_class_is_where_the_flc_route_earns_its_keep():
    """functional FRI + weak FLC: a FRI-only rule calls these LATE and they are mostly EARLY."""
    r = json.loads(ART.read_text(encoding="utf-8"))
    da = next(c for c in r["rule_cells"] if c["is_da112_class"])
    strong = next(c for c in r["rule_cells"] if c["fri"] == "functional" and c["flc"] == "strong")
    assert da["pct_late"] < 0.5                       # mostly early -> FRI-only is wrong on them
    assert strong["pct_late"] > 0.8                   # ...while strong-FLC siblings really are late
    assert strong["pct_late"] - da["pct_late"] > 0.35  # FLC separates them materially


@pytest.mark.skipif(not ART.exists(), reason="FLC-route artifact not present")
def test_flc_adds_over_fri_only_within_ancestry_not_just_pooled():
    """The pooled gain could be ancestry confounding (the S3 run measured exactly that collapse)."""
    gw = json.loads(ART.read_text(encoding="utf-8"))["population_structure"]["group_weighted"]
    assert gw["mean_acc_two_locus"] > gw["mean_acc_fri_only"] > gw["mean_null"]
    assert gw["mean_delta"] > 0


@pytest.mark.skipif(not ART.exists(), reason="FLC-route artifact not present")
def test_the_threshold_dependence_is_recorded_not_buried():
    """The gain is NOT threshold-robust -- it dies at q60 and reverses at q70. A reader must see that."""
    sens = json.loads(ART.read_text(encoding="utf-8"))["threshold_sensitivity"]
    assert len(sens) >= 5
    assert any(x["delta"] > 0.02 for x in sens)       # it helps in the plausible (low-quantile) range...
    assert any(x["delta"] <= 0 for x in sens)         # ...and hurts if you over-call weak FLC


@pytest.mark.skipif(not ART.exists(), reason="FLC-route artifact not present")
def test_the_expression_proxy_mapping_is_declared():
    """FLC EXPRESSION is not an allele call; the artifact must say so rather than imply equivalence."""
    s = json.loads(ART.read_text(encoding="utf-8"))["substrate"]
    assert "PROXY" in s["honest_mapping"]
    assert "Atwell" in s["citations"]["flc_expression"]
