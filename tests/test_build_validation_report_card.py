"""Pin the report-card cell-state machine (scripts/build_validation_report_card.py).

The 6-state classifier is the load-bearing honesty surface of Anchor-4: a mis-classified cell would let
"validated" drift (e.g. an underpowered cell rendering as scored, or an other-kingdom decoder claiming a
phenotype source it doesn't have). These tests pin each state + the precedence order.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

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


def test_surface_no_free_source_state():
    """A surface cell flagged no_free_source classifies NO_FREE_PHENOTYPE_SOURCE (structural non-cell)."""
    key = ("candida_auris", "fluconazole")
    surface = {"phenotype_source_status": "no_free_source", "engine": "fungal_erg11"}
    c = mod.classify(key, {}, {}, {}, surface)
    assert c["state"] == "NO_FREE_PHENOTYPE_SOURCE"


def test_surface_label_confounded_state():
    """oxacillin/S. aureus -> LABEL_CONFOUNDED (M2), distinct from NOT_CENSUSED."""
    key = ("staphylococcus_aureus", "oxacillin")
    surface = {"phenotype_source_status": "label_confounded"}
    c = mod.classify(key, {}, {}, {}, surface)
    assert c["state"] == "LABEL_CONFOUNDED"


def test_label_confounded_precedence_over_scored():
    """A confounded label must NOT be presented as a clean SCORED number — structural property wins."""
    key = ("staphylococcus_aureus", "oxacillin")
    surface = {"phenotype_source_status": "label_confounded"}
    c = mod.classify(key, {key: _scored_cell()}, {}, {}, surface)
    assert c["state"] == "LABEL_CONFOUNDED"


# --- main() end-to-end emit (read-only roll-up; redirect WIKI/ROOT to tmp so real artifacts aren't clobbered) ---

import json  # noqa: E402


def _redirect_io(monkeypatch, tmp_path):
    """Point the module's WIKI (outputs + scored/census reads) + ROOT (registry read) at empty tmp dirs."""
    wiki = tmp_path / "wiki"
    wiki.mkdir()
    data = tmp_path / "dna_decode" / "data"
    data.mkdir(parents=True)
    monkeypatch.setattr(mod, "WIKI", wiki)
    monkeypatch.setattr(mod, "ROOT", tmp_path)
    return wiki


def test_main_emits_json_and_md_with_no_observations(monkeypatch, tmp_path):
    """With NO scored/census/registry files on disk, main() still emits both artifacts and every row comes
    from the shipped surface (a new decoder cannot ship invisibly)."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    rc = mod.main()
    assert rc == 0
    j = wiki / "decoder_validation_report_card.json"
    md = wiki / "decoder_validation_report_card.md"
    assert j.exists() and md.exists()
    doc = json.loads(j.read_text(encoding="utf-8"))
    assert doc["_schema"] == "decoder-validation-report-card-v0"
    assert doc["no_aggregate_headline"] is True
    # surface-only run: structural-label cells classify without any observation files
    cells = {(c["organism"], c["drug"]): c for c in doc["cells"]}
    assert cells[("candida_auris", "fluconazole")]["state"] == "NO_FREE_PHENOTYPE_SOURCE"
    assert cells[("staphylococcus_aureus", "oxacillin")]["state"] == "LABEL_CONFOUNDED"
    # an ncbi_pd surface cell with no census renders NOT_CENSUSED, never silently dropped
    assert cells[("escherichia_coli_shigella", "ciprofloxacin")]["state"] == "NOT_CENSUSED"
    assert sum(doc["state_counts"].values()) == len(doc["cells"])


def test_main_scored_json_renders_in_grid(monkeypatch, tmp_path):
    """A provenance_disjoint_validation_*.json on disk surfaces as a SCORED row in the emitted markdown."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provenance_disjoint_validation_kleb_cipro_2026-06-10.json").write_text(json.dumps({
        "organism": "Klebsiella", "drug": "ciprofloxacin",
        "metrics": {"acc": 0.95, "sens": 0.93, "spec": 0.97, "n_scored": 60,
                    "tp": 28, "fp": 1, "tn": 29, "fn": 2},
        "independence_tier": "provenance-disjoint ...",
    }), encoding="utf-8")
    rc = mod.main()
    assert rc == 0
    doc = json.loads((wiki / "decoder_validation_report_card.json").read_text(encoding="utf-8"))
    kleb = next(c for c in doc["cells"] if (c["organism"], c["drug"]) == ("klebsiella", "ciprofloxacin"))
    assert kleb["state"] == "SCORED" and kleb["acc"] == 0.95 and kleb["n"] == 60
    md_text = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "`SCORED`" in md_text and "0.95" in md_text


# --- lineage-disclosure layer (Step 3) ---


def _scored_lineage_cell(grade="clonal (<3 effective lineages)"):
    return {
        "organism": "Klebsiella", "drug": "ciprofloxacin", "raw_N": 60,
        "lineage_tier_emitted": True, "lineage_grade": grade,
        "thresholds": {
            "0.001": {"effective_lineage_N_R": 5, "effective_lineage_N_S": 12,
                      "cluster_weighted": {"sens": 0.8, "sens_ci": [0.3, 0.99], "sens_eff_n": 5,
                                           "spec": 1.0, "spec_ci": [0.7, 1.0], "spec_eff_n": 12,
                                           "n_discordant": 1}},
            "0.005": {"effective_lineage_N_R": 2, "effective_lineage_N_S": 8,
                      "cluster_weighted": {"sens": 0.5, "sens_ci": [0.09, 0.91], "sens_eff_n": 2,
                                           "spec": 1.0, "spec_ci": [0.6, 1.0], "spec_eff_n": 8,
                                           "n_discordant": 2}},
        },
    }


def test_c3_emitter_guard_refuses_weighted_without_ci():
    """A cluster-weighted point estimate with no Wilson CI is a honesty inversion — must raise (C3)."""
    with pytest.raises(AssertionError):
        mod._assert_weighted_renderable({"sens": 0.5, "sens_eff_n": 2})  # no sens_ci
    with pytest.raises(AssertionError):
        mod._assert_weighted_renderable({"sens": 0.5, "sens_ci": [0.1, 0.9]})  # no eff_n
    # a None metric needs no CI (nothing to render)
    mod._assert_weighted_renderable({"sens": None, "spec": None})


def test_build_lineage_block_states():
    assert mod.build_lineage_block(None)["status"] == "not_computed"
    inc = mod.build_lineage_block({"partial": True, "n_genomes_missing": 6, "raw_N": 54,
                                   "lineage_tier_emitted": False})
    assert inc["status"] == "incomplete" and inc["n_genomes_missing"] == 6
    sc = mod.build_lineage_block(_scored_lineage_cell())
    assert sc["status"] == "scored" and sc["effective_lineage_N"]["0.005"] == {"R": 2, "S": 8}


def test_main_renders_lineage_columns_with_ci(monkeypatch, tmp_path):
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provenance_disjoint_validation_kleb_cipro_2026-06-10.json").write_text(json.dumps({
        "organism": "Klebsiella", "drug": "ciprofloxacin",
        "metrics": {"acc": 0.967, "sens": 0.967, "spec": 0.967, "n_scored": 60,
                    "tp": 29, "fp": 1, "tn": 29, "fn": 1},
        "independence_tier": "x",
    }), encoding="utf-8")
    (wiki / "provdisjoint_lineage_metrics.json").write_text(json.dumps({
        "_schema": "provdisjoint-lineage-metrics-v1", "cells": [_scored_lineage_cell()],
    }), encoding="utf-8")
    assert mod.main() == 0
    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "Lineage disclosure" in md
    assert "0.5 [0.09–0.91] (n=2)" in md  # weighted sens @0.005 with CI + eff-N
    assert "clonal (<3 effective lineages)" in md
    doc = json.loads((wiki / "decoder_validation_report_card.json").read_text(encoding="utf-8"))
    kleb = next(c for c in doc["cells"] if (c["organism"], c["drug"]) == ("klebsiella", "ciprofloxacin"))
    assert kleb["state"] == "SCORED" and kleb["lineage"]["status"] == "scored"  # SCORED not removed


def test_main_scored_without_lineage_renders_not_computed(monkeypatch, tmp_path):
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provenance_disjoint_validation_kleb_cipro_2026-06-10.json").write_text(json.dumps({
        "organism": "Klebsiella", "drug": "ciprofloxacin",
        "metrics": {"acc": 0.95, "sens": 0.93, "spec": 0.97, "n_scored": 60,
                    "tp": 28, "fp": 1, "tn": 29, "fn": 2},
        "independence_tier": "x",
    }), encoding="utf-8")
    # NO lineage sidecar on disk
    assert mod.main() == 0
    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "lineage: not computed" in md  # never silently blank


def test_main_partial_lineage_renders_incomplete(monkeypatch, tmp_path):
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provenance_disjoint_validation_kleb_tetra_2026-06-10.json").write_text(json.dumps({
        "organism": "Klebsiella", "drug": "tetracycline",
        "metrics": {"acc": 0.9, "sens": 0.9, "spec": 0.9, "n_scored": 33,
                    "tp": 15, "fp": 2, "tn": 14, "fn": 2},
        "independence_tier": "x",
    }), encoding="utf-8")
    (wiki / "provdisjoint_lineage_metrics.json").write_text(json.dumps({
        "_schema": "provdisjoint-lineage-metrics-v1",
        "cells": [{"organism": "Klebsiella", "drug": "tetracycline", "raw_N": 33,
                   "partial": True, "n_genomes_missing": 27, "lineage_tier_emitted": False}],
    }), encoding="utf-8")
    assert mod.main() == 0
    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "lineage: incomplete (27 genomes missing)" in md
