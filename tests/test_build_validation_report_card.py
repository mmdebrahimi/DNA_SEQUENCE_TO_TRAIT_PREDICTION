"""Pin the report-card cell-state machine (scripts/build_validation_report_card.py).

The 6-state classifier is the load-bearing honesty surface of Anchor-4: a mis-classified cell would let
"validated" drift (e.g. an underpowered cell rendering as scored, or an other-kingdom decoder claiming a
phenotype source it doesn't have). These tests pin each state + the precedence order.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "build_validation_report_card",
    Path(__file__).resolve().parent.parent / "scripts" / "build_validation_report_card.py",
)
mod = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(mod)


def _scored_cell():
    return {"metrics": {"acc": 0.97, "sens": 0.97, "spec": 0.97, "n_scored": 60,
                        "tp": 29, "fp": 1, "tn": 29, "fn": 1, "abstain": 0},
            "independence_tier": "provenance-disjoint ...", "_file": "x.json"}


def test_scored_state():
    key = ("klebsiella", "ciprofloxacin")
    c = mod.classify(key, {key: _scored_cell()}, {}, {})
    assert c["state"] == "SCORED" and c["acc"] == 0.97 and c["n"] == 60


def test_powered_unscored_state():
    key = ("klebsiella", "ceftriaxone")
    census = {key: {"organism": "Klebsiella", "drug": "ceftriaxone", "other_R": 505, "other_S": 410, "powered": True}}
    c = mod.classify(key, {}, census, {})
    assert c["state"] == "POWERED_UNSCORED" and "505R/410S" in c["note"]


def test_underpowered_state():
    key = ("salmonella", "ciprofloxacin")
    census = {key: {"organism": "Salmonella", "drug": "ciprofloxacin", "other_R": 4, "other_S": 87, "powered": False}}
    c = mod.classify(key, {}, census, {})
    assert c["state"] == "UNDERPOWERED"


def test_abstains_by_design_state():
    key = ("acinetobacter", "meropenem")
    registry = {key: {"verdict": "EXPRESSION_FLOOR", "counter": "broad", "threshold": 1}}
    c = mod.classify(key, {}, {}, registry)
    assert c["state"] == "ABSTAINS_BY_DESIGN"


def test_not_censused_state():
    c = mod.classify(("morganella", "ciprofloxacin"), {}, {}, {})
    assert c["state"] == "NOT_CENSUSED"


def test_scored_takes_precedence_over_census_and_registry():
    """A scored JSON must win even if census/registry also have the key (scored is ground truth)."""
    key = ("klebsiella", "ciprofloxacin")
    census = {key: {"organism": "K", "drug": "c", "other_R": 4, "other_S": 4, "powered": False}}
    registry = {key: {"verdict": "EXPRESSION_FLOOR", "counter": "broad", "threshold": 1}}
    c = mod.classify(key, {key: _scored_cell()}, census, registry)
    assert c["state"] == "SCORED"


def test_abstains_precedence_over_census():
    """EXPRESSION_FLOOR abstention outranks a powered census — an abstaining rule isn't 'unscored', it's a no-op by design."""
    key = ("acinetobacter", "meropenem")
    census = {key: {"organism": "A", "drug": "m", "other_R": 99, "other_S": 99, "powered": True}}
    registry = {key: {"verdict": "EXPRESSION_FLOOR", "counter": "broad", "threshold": 1}}
    c = mod.classify(key, {}, census, registry)
    assert c["state"] == "ABSTAINS_BY_DESIGN"
