"""Offline test for the AMR Portal independent report-card builder (pure build())."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.build_amr_portal_report_card import build, render_experimental_md  # noqa: E402


def test_build_tiers_and_routing():
    scores = {
        "Salmonella enterica|ciprofloxacin": {"organism": "Salmonella enterica", "drug": "ciprofloxacin",
            "n_R": 2434, "n_S": 22538, "sens": 0.907, "spec": 0.964, "accuracy": 0.959,
            "sens_wilson95": [0.895, 0.918], "spec_wilson95": [0.962, 0.967], "powered": True, "n_indeterminate": 0},
        "Salmonella enterica|meropenem": {"organism": "Salmonella enterica", "drug": "meropenem",
            "n_R": 0, "n_S": 19464, "sens": None, "spec": 1.0, "accuracy": 1.0,
            "sens_wilson95": None, "spec_wilson95": [1.0, 1.0], "powered": False, "n_indeterminate": 0},
    }
    card = build(scores)
    assert card["n_cells"] == 2 and card["n_scored_independent"] == 1 and card["n_underpowered"] == 1
    by = {(c["organism"], c["drug"]): c for c in card["cells"]}
    assert by[("Salmonella enterica", "ciprofloxacin")]["tier"] == "SCORED_INDEPENDENT"
    assert by[("Salmonella enterica", "ciprofloxacin")]["routing"] == "calibrated_registry"   # opt-in cipro
    assert by[("Salmonella enterica", "meropenem")]["tier"] == "UNDERPOWERED"                  # 0 R
    assert by[("Salmonella enterica", "meropenem")]["routing"] == "drug_rule_default"


def test_experimental_cells_are_namespace_separate():
    """The TMP-SMX experimental overlay must NOT leak into the deployed-surface cells/counts
    (the shared-key silent-overwrite trap). build() knows only deployed scores; the experimental
    block is a separate render + a separate top-level key wired in main()."""
    scores = {"Escherichia coli|ciprofloxacin": {"organism": "Escherichia coli", "drug": "ciprofloxacin",
        "n_R": 60, "n_S": 60, "sens": 0.9, "spec": 0.95, "accuracy": 0.92,
        "sens_wilson95": [0.8, 0.95], "spec_wilson95": [0.9, 0.98], "powered": True, "n_indeterminate": 0}}
    card = build(scores)
    assert card["n_cells"] == 1 and card["n_scored_independent"] == 1
    assert "experimental" not in card            # build() never mixes experimental in
    # the experimental renderer brands its own section + never touches deployed counts
    exp = {"n_scored": 1, "cells": {"Escherichia coli|trimethoprim-sulfamethoxazole": {
        "organism": "Escherichia coli", "drug": "trimethoprim-sulfamethoxazole", "headline": "SCORED",
        "n_R": 2619, "n_S": 9269, "strata_reproduced": True,
        "binary": {"sens": 0.727, "spec": 0.983, "acc": 0.926}}}}
    md = "\n".join(render_experimental_md(exp))
    assert "EXPERIMENTAL_SCORED" in md and "NOT counted in the deployed-surface totals" in md
    assert "trimethoprim-sulfamethoxazole" in md


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
