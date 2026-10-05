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


def test_invisible_fraction_from_metrics():
    # fn / (tp + fn) = 1 - sens
    assert mod.invisible_fraction_from_metrics({"tp": 29, "fn": 1}) == round(1 / 30, 3)
    assert mod.invisible_fraction_from_metrics({"tp": 11, "fn": 23}) == round(23 / 34, 3)  # gono-tet shape
    assert mod.invisible_fraction_from_metrics({"tp": 20, "fn": 0}) == 0.0  # fully visible
    assert mod.invisible_fraction_from_metrics({"tp": 0, "fn": 0}) is None  # no measured-R scored
    assert mod.invisible_fraction_from_metrics({"tp": None, "fn": None}) is None  # missing counts


def test_scored_cell_carries_invisible_fraction():
    key = ("klebsiella", "ciprofloxacin")
    c = mod.classify(key, {key: _scored_cell()}, {}, {})  # tp=29 fn=1 -> 1/30
    assert c["invisible_fraction"] == round(1 / 30, 3)


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


def test_build_lineage_block_unreconciled_not_partial_is_not_computed():
    """A cell that didn't emit a tier but is NOT partial (e.g. reconcile failed) -> not_computed,
    never 'incomplete' (incomplete is reserved for genome-completeness gaps)."""
    blk = mod.build_lineage_block({"partial": False, "lineage_tier_emitted": False, "raw_N": 60,
                                   "n_genomes_missing": 0})
    assert blk["status"] == "not_computed" and blk["raw_N"] == 60


def test_load_lineage_metrics_reads_and_keys_by_canonical(tmp_path, monkeypatch):
    wiki = tmp_path / "wiki"
    wiki.mkdir()
    (wiki / mod.LINEAGE_SIDECAR).write_text(json.dumps({
        "_schema": "provdisjoint-lineage-metrics-v1",
        "cells": [{"organism": "Klebsiella", "drug": "Ciprofloxacin", "raw_N": 60}],
    }), encoding="utf-8")
    monkeypatch.setattr(mod, "WIKI", wiki)
    got = mod.load_lineage_metrics()
    assert got[("klebsiella", "ciprofloxacin")]["raw_N"] == 60  # canonical-keyed


def test_load_lineage_metrics_absent_is_empty(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "WIKI", tmp_path)  # no sidecar on disk
    assert mod.load_lineage_metrics() == {}


def test_load_lineage_metrics_malformed_is_empty(tmp_path, monkeypatch):
    wiki = tmp_path / "wiki"
    wiki.mkdir()
    (wiki / mod.LINEAGE_SIDECAR).write_text("{ not json", encoding="utf-8")
    monkeypatch.setattr(mod, "WIKI", wiki)
    assert mod.load_lineage_metrics() == {}  # malformed must not break the read-only roll-up


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


def test_naive_value_add_loader_and_section(monkeypatch, tmp_path):
    """The curated-vs-naive value-add (wrapper-vs-tool rail) must surface on the standing card; reconciled
    cells render in the table, RECONCILE_MISMATCH cells are excluded."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provdisjoint_naive_comparator_2026-06-27.json").write_text(json.dumps({
        "_schema": "provdisjoint-naive-comparator-v1",
        "cells": {
            "Klebsiella:ciprofloxacin": {"frozen_balacc": 0.967, "naive_balacc": 0.7,
                                         "delta_balacc": 0.267, "value_add_verdict": "CURATED_LAYER_ADDS_VALUE"},
            "Klebsiella:gentamicin": {"value_add_verdict": "RECONCILE_MISMATCH"},
        },
    }), encoding="utf-8")
    rows = mod.load_naive_value_add()
    keys = {(r["organism"], r["drug"], r["verdict"]) for r in rows}
    assert ("Klebsiella", "ciprofloxacin", "CURATED_LAYER_ADDS_VALUE") in keys
    assert mod.main() == 0
    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "Curated layer value-over-naive-baseline" in md
    assert "CURATED_LAYER_ADDS_VALUE" in md
    assert "RECONCILE_MISMATCH" not in md  # mismatch cells excluded from the rendered table


# --- prospective-lock disclosure layer (2026-08-24) ---


def _live_surface_hashes():
    """The current frozen-surface hashes, so a fixture artifact is one the pipeline could really emit."""
    from dna_decode.eval.prospective_lock import surface_hashes
    return surface_hashes()


def _prospective_artifact(org="Klebsiella", drug="ciprofloxacin", *, acc=0.60, sens=0.40,
                          spec=0.90, n=50, generated="2026-08-24", verified=True):
    return {
        "artifact": "prospective_lock_validation", "generated": generated,
        "organism": org, "drug": drug, "prospective_lock_verified": verified,
        # A REAL artifact stamps the manifest's whole `surface_sha256` (see build_artifact in
        # scripts/prospective_lock_validate.py). The fixture carried no hashes at all, which made it an
        # input the pipeline cannot produce -- and once the re-verification stopped blessing unpinnable
        # artifacts (2026-09-22) it correctly withheld these, exercising the wrong property. Stamping
        # the LIVE surface keeps each test aimed at what it was written for.
        "lock_manifest": {"lock_date": "2026-06-13",
                          "surface_sha256": _live_surface_hashes()},
        "confusion": {"n_scored": n, "acc": acc, "sens": sens, "spec": spec, "abstain": 0},
        "powering": {"status": "POWERED", "n_scored": n, "scored_R": 20, "scored_S": 30},
    }


def test_prospective_augments_a_scored_cell_without_overwriting_its_provdisjoint_numbers(
        monkeypatch, tmp_path):
    """THE shared-key trap: a prospective cell shares (organism, drug) with a provdisjoint cell.

    Merging them would silently replace one number with the other. This pins that the provdisjoint
    acc/n survive UNCHANGED while the prospective figures live in their own block -- deliberately using
    DIFFERENT values so an overwrite could not pass by coincidence.
    """
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provenance_disjoint_validation_kleb_cipro_2026-06-10.json").write_text(json.dumps({
        "organism": "Klebsiella", "drug": "ciprofloxacin",
        "metrics": {"acc": 0.95, "sens": 0.93, "spec": 0.97, "n_scored": 60,
                    "tp": 28, "fp": 1, "tn": 29, "fn": 2},
    }), encoding="utf-8")
    (wiki / "prospective_lock_validation_Klebsiella_ciprofloxacin_2026-08-24.json").write_text(
        json.dumps(_prospective_artifact()), encoding="utf-8")

    assert mod.main() == 0
    doc = json.loads((wiki / "decoder_validation_report_card.json").read_text(encoding="utf-8"))
    cell = next(c for c in doc["cells"] if (c["organism"], c["drug"]) == ("klebsiella", "ciprofloxacin"))

    assert cell["state"] == "SCORED"            # prospective AUGMENTS; it never demotes a state
    assert cell["acc"] == 0.95 and cell["n"] == 60          # provdisjoint numbers untouched
    p = cell["prospective"]
    assert p["status"] == "scored" and p["acc"] == 0.60 and p["n_scored"] == 50
    assert p["lock_date"] == "2026-06-13"

    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "Prospective-lock disclosure" in md
    assert "leakage-free BY CONSTRUCTION" in md
    assert "0.95" in md and "0.600" in md       # BOTH arms rendered, neither replaced


def test_prospective_section_absent_when_nothing_has_accrued(monkeypatch, tmp_path):
    """Non-vacuity: the section must not render on an empty accrual, or it would imply a number exists."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    assert mod.main() == 0
    assert "Prospective-lock disclosure" not in (
        wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")


def test_load_prospective_keeps_the_newest_artifact_per_cell(monkeypatch, tmp_path):
    """Cells are RE-scored as the cohort accrues, so several dated artifacts coexist."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    for gen, acc in (("2026-08-01", 0.10), ("2026-09-01", 0.90), ("2026-07-01", 0.50)):
        (wiki / f"prospective_lock_validation_Klebsiella_ciprofloxacin_{gen}.json").write_text(
            json.dumps(_prospective_artifact(acc=acc, generated=gen)), encoding="utf-8")
    got = mod.load_prospective()
    assert len(got) == 1
    assert next(iter(got.values()))["confusion"]["acc"] == 0.90     # newest wins, not last-globbed


def test_prospective_block_refuses_to_render_an_unverified_lock():
    """A stale artifact from a DRIFTED decoder must never read as validating the current one."""
    blk = mod.build_prospective_block(_prospective_artifact(verified=False))
    assert blk["status"] == "lock_unverified"
    assert "acc" not in blk                                        # no number leaks out of an unverified lock
    assert mod.build_prospective_block(None)["status"] == "not_accrued"


def test_prospective_prose_names_the_ACTIVE_lock_and_never_a_retired_one(monkeypatch, tmp_path):
    """A hardcoded manifest filename in the prose goes stale SILENTLY when the surface is revised.

    Retiring a lock IS the act of making its manifest stop verifying, so nothing re-syncs a restated
    filename. This card named the retired 2026-06-22 v1 manifest while rendering a cell scored against
    the 2026-08-31 v2 lock -- the standing trust surface citing the wrong pinning authority. The name
    must be DERIVED from resolve_active_lock(), and no other manifest may appear in the prose.
    """
    import re

    from dna_decode.eval.prospective_lock import resolve_active_lock

    active = resolve_active_lock()[0].name

    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "prospective_lock_validation_Klebsiella_ciprofloxacin_2026-08-24.json").write_text(
        json.dumps(_prospective_artifact()), encoding="utf-8")
    assert mod.main() == 0
    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")

    assert active in md, f"prose must name the active lock {active}"
    # Non-vacuity: any manifest it mentions must BE the active one, so a future hardcoded
    # filename fails here instead of quietly co-existing with the derived name.
    mentioned = set(re.findall(r"prospective_lock_manifest_[0-9-]+\.json", md))
    assert mentioned == {active}, f"prose mentions non-active manifest(s): {mentioned - {active}}"


def test_a_powered_prospective_regression_raises_a_TOP_LEVEL_flag(monkeypatch, tmp_path):
    """A consumer filtering `state == SCORED` must not be able to miss a contradicting prospective result.

    The state is deliberately NOT demoted -- the provenance-disjoint result is still what it was. The flag
    says something different: this cell has a prospective result that contradicts its standing.
    """
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provenance_disjoint_validation_kleb_gent_2026-06-10.json").write_text(json.dumps({
        "organism": "Klebsiella", "drug": "gentamicin",
        "metrics": {"acc": 0.95, "sens": 0.93, "spec": 0.97, "n_scored": 60,
                    "tp": 28, "fp": 1, "tn": 29, "fn": 2},
    }), encoding="utf-8")
    (wiki / "prospective_lock_validation_Klebsiella_gentamicin_2026-08-24.json").write_text(
        json.dumps(_prospective_artifact(org="Klebsiella", drug="gentamicin",
                                         acc=0.53, sens=0.429, spec=0.92, n=62)), encoding="utf-8")
    assert mod.main() == 0
    doc = json.loads((wiki / "decoder_validation_report_card.json").read_text(encoding="utf-8"))
    cell = next(c for c in doc["cells"] if (c["organism"], c["drug"]) == ("klebsiella", "gentamicin"))

    assert cell["state"] == "SCORED"                 # NOT demoted -- augment, never demote
    assert cell["prospective_regression"] is True
    assert "under-calls" in cell["deployment_caveat"]
    assert cell["prospective"]["regression"] is True
    md = (wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "**REGRESSION**" in md


def test_a_healthy_or_underpowered_prospective_cell_raises_no_flag(monkeypatch, tmp_path):
    """Non-vacuity in both directions: a good result and an unpowered one must both stay silent."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "prospective_lock_validation_Klebsiella_ciprofloxacin_2026-08-24.json").write_text(
        json.dumps(_prospective_artifact(sens=0.917, acc=0.967)), encoding="utf-8")
    assert mod.main() == 0
    doc = json.loads((wiki / "decoder_validation_report_card.json").read_text(encoding="utf-8"))
    cell = next(c for c in doc["cells"] if c.get("prospective"))
    assert cell.get("prospective_regression") is None and cell["prospective"]["regression"] is False

    # an UNDERPOWERED prospective cell makes no claim either way, even at a terrible sens
    blk = mod.build_prospective_block({
        "prospective_lock_verified": True, "generated": "2026-08-24",
        "lock_manifest": {"lock_date": "2026-06-13",
                          "surface_sha256": _live_surface_hashes()},
        "confusion": {"n_scored": 4, "acc": 0.1, "sens": 0.1, "spec": 0.1, "abstain": 0},
        "powering": {"status": "UNDERPOWERED", "scored_R": 2, "scored_S": 2}})
    assert blk["regression"] is False


# --- source-concentration disclosure layer ---

def test_source_concentration_never_overwrites_a_provdisjoint_metric():
    """ANTI-OVERWRITE, with deliberately DIFFERENT values so a merge bug cannot pass by coincidence.

    A source-concentration row shares its (organism, drug) key with a provdisjoint cell. Feeding it into
    load_scored() would silently replace one number with the other -- the documented shared-key trap. This
    asserts the block ATTACHES and the metrics survive untouched.
    """
    from scripts.build_validation_report_card import build_source_block
    blk = build_source_block({"organism": "x", "drug": "y", "n_cohort": 60, "spec": 0.111,
                              "bioproject": {"distinct": 2, "largest": ["PRJNA1", 58],
                                             "largest_share": 0.967, "n_unknown": 0},
                              "sra_center": {"distinct": 2}})
    assert blk["status"] == "measured"
    assert blk["single_source"] is True
    # the block must not carry a competing metric that a consumer could mistake for THE cell metric
    assert "spec" not in blk and "sens" not in blk and "acc" not in blk


def test_single_source_flag_is_a_disclosure_not_a_demotion():
    """A flagged cell keeps its state. The flag says the estimate rests on one source, which is a
    different fact from the estimate being wrong -- and the error is not even directional."""
    import json
    from pathlib import Path as _P
    card = _P(__file__).resolve().parent.parent / "wiki" / "decoder_validation_report_card.json"
    if not card.exists():
        import pytest
        pytest.skip("card absent")
    cells = json.loads(card.read_text(encoding="utf-8"))["cells"]
    flagged = [c for c in cells if c.get("source_concentration", {}).get("single_source")]
    if not flagged:
        import pytest
        pytest.skip("no single-source cells in the current card")
    for c in flagged:
        assert c["state"] == "SCORED", "the disclosure must not demote a cell"
        assert c.get("sens") is not None, "metrics must survive the disclosure"


def test_an_incomplete_provenance_sweep_renders_nothing():
    """A partial sweep cannot support a concentration claim. Rendering a floor that reads like a
    measurement is the failure this refuses -- same rule as the fail-closed leakage manifest."""
    import json
    from unittest import mock
    from pathlib import Path as _P
    from scripts import build_validation_report_card as B
    payload = json.dumps({"complete": False, "cells": [{"organism": "x", "drug": "y"}]})
    with mock.patch.object(_P, "exists", lambda self: True),          mock.patch.object(_P, "read_text", lambda self, **k: payload):
        assert B.load_source_concentration() == {}


def test_unknown_provenance_is_surfaced_not_hidden():
    """n_unknown must reach the block: a cell whose provenance is mostly missing looks diverse only
    because the metadata is absent, and the reader has to be able to see that."""
    from scripts.build_validation_report_card import build_source_block
    blk = build_source_block({"organism": "x", "drug": "y", "n_cohort": 60,
                              "bioproject": {"distinct": 3, "largest": ["P", 5],
                                             "largest_share": 0.5, "n_unknown": 50},
                              "sra_center": {"distinct": 1}})
    assert blk["n_unknown_provenance"] == 50


# --- L2 doubt disclosure layer (added 2026-08-31) -------------------------------------------------

def test_doubt_block_never_reports_absence_as_clean():
    """An unscreened drug must be distinguishable from a screened-and-clean one."""
    build_doubt_block = mod.build_doubt_block
    b = build_doubt_block(None)
    assert b["status"] == "not_measured"
    assert "unassessable" in b["note"]


def test_doubt_block_stamps_its_drug_level_scope():
    """The screen is index-wide, not per-cohort. Without the scope stamp a drug-level result would
    read as a statement about one cell's cohort — which it is not."""
    build_doubt_block = mod.build_doubt_block
    b = build_doubt_block({"drug": "gentamicin", "status": "scored", "n_families_uncounted": 131,
                           "n_strong": 1, "n_raw_signature": 1, "known_gap_recovered": True,
                           "strong": [{"evidence": {"symbol": "rmtE1"}}]})
    assert b["status"] == "measured"
    assert "NOT this cell's cohort" in b["scope"]
    assert b["strong_families"] == ["rmtE1"] and b["known_gap_recovered"] is True


def test_doubt_layer_is_augment_only_on_the_committed_card():
    """AUGMENT-ONLY, checked against the artifact: the block lives under its own key and carries no
    state, tier or metric. Merging it into a cell's own fields is the shared-key overwrite trap."""
    import json
    from pathlib import Path
    p = Path("wiki/decoder_validation_report_card.json")
    if not p.exists():
        pytest.skip("report card artifact absent")
    cells = json.loads(p.read_text(encoding="utf-8"))["cells"]
    with_doubt = [c for c in cells if "doubt_layer" in c]
    assert with_doubt, "no cell carries a doubt_layer — this guard would be vacuous"
    for c in with_doubt:
        d = c["doubt_layer"]
        forbidden = {"state", "tier", "acc", "sens", "spec", "n_scored", "cell_state"}
        assert not (forbidden & set(d)), f"{c['organism']}x{c['drug']}: doubt block carries {forbidden & set(d)}"


def test_the_strong_family_set_tracks_the_deployed_rule():
    """Pinned `{"rmtE1"}` until 2026-08-31; now EMPTY because the v2 lock deployed the fix.

    If a family ever goes STRONG that is a real finding needing adjudication, so this must fail loudly
    rather than pass unnoticed. Expressed against the deployed rule so it cannot drift: while the
    gentamicin rescue ships, `rmtE1` is counted and must NOT appear as an uncounted gap.
    """
    import json
    from pathlib import Path

    from dna_decode.eval.amr_rules import rule_for
    p = Path("wiki/decoder_validation_report_card.json")
    if not p.exists():
        pytest.skip("report card artifact absent")
    cells = json.loads(p.read_text(encoding="utf-8"))["cells"]
    fams = {f for c in cells for f in (c.get("doubt_layer", {}).get("strong_families") or [])}
    if rule_for("gentamicin").get("symbol_rescue"):
        assert "rmtE1" not in fams, "the rescue ships but rmtE1 still reads as an uncounted gap"
        assert fams == set(), f"a NEW strong completeness family appeared -- adjudicate it: {sorted(fams)}"
    else:
        assert fams == {"rmtE1"}, f"the set of STRONG completeness families changed: {sorted(fams)}"


# --- the prospective re-verification must not bless what it cannot check (found in review 2026-09-22)

def _stamped(hashes):
    return {"prospective_lock_verified": True, "generated": "2026-09-01",
            "lock_manifest": {"lock_date": "2026-08-31", "surface_sha256": hashes},
            "confusion": {"tp": 10, "fp": 0, "tn": 10, "fn": 0},
            "powering": {"powered": True}}


def test_an_artifact_that_pins_nothing_is_withheld_not_rendered():
    """REGRESSION. `if stamped:` was falsy for an absent/empty surface_sha256, so the ENTIRE
    re-verification was skipped and the cell rendered as a LIVE prospective score -- the strongest
    tier in the project, granted to an artifact that pinned nothing."""
    from scripts.build_validation_report_card import build_prospective_block
    assert build_prospective_block(_stamped({}))["status"] == "lock_unverifiable"


def test_a_partially_pinned_artifact_whose_hashes_all_match_is_withheld():
    """A partial stamp was checked only over the files it happened to carry, so the decoder could have
    changed in any UNPINNED file while the cell still read as verified.

    Scoped to the all-match case on purpose: if a pinned hash MISmatches we have positive proof of a
    revision, and `superseded_by_surface_change` is the truer, more informative status. Both withhold."""
    from dna_decode.eval.prospective_lock import FROZEN_SURFACE_FILES, surface_hashes
    from scripts.build_validation_report_card import build_prospective_block
    live = surface_hashes()
    partial = {FROZEN_SURFACE_FILES[0]: live[FROZEN_SURFACE_FILES[0]]}
    blk = build_prospective_block(_stamped(partial))
    assert blk["status"] == "lock_unverifiable"
    assert blk["unpinned_files"], "the reason must be reported, not just the refusal"


def test_a_fully_pinned_matching_artifact_still_scores():
    """NON-VACUITY: the tightening must reject the vacuous forms WITHOUT rejecting a genuine one. A
    guard that withholds everything is not fail-closed, it is broken."""
    from dna_decode.eval.prospective_lock import surface_hashes
    from scripts.build_validation_report_card import build_prospective_block
    blk = build_prospective_block(_stamped(surface_hashes()))
    assert blk["status"] not in ("lock_unverifiable", "superseded_by_surface_change")


def test_drift_is_still_detected_on_a_complete_pin():
    from dna_decode.eval.prospective_lock import FROZEN_SURFACE_FILES, surface_hashes
    from scripts.build_validation_report_card import build_prospective_block
    h = dict(surface_hashes()); h[FROZEN_SURFACE_FILES[0]] = "0" * 64
    assert build_prospective_block(_stamped(h))["status"] == "superseded_by_surface_change"


# ---------------------------------------------------------------------------
# The COMMITTED card vs the LIVE registry.
#
# Every test above builds through _redirect_io, which points the module's WIKI and ROOT at
# EMPTY tmp dirs -- so they pin the classifier's logic on synthetic inputs and say nothing
# about the artifact that actually SHIPS. That is the same structural blind spot that let
# wiki/certification_capstone.json sit in git reporting independent_measured=31 against a
# registry returning 33 (fixed in 1ff6020): its tests rebuilt before asserting, so drift was
# invisible and silently repaired.
#
# This card is exposed the same way -- it is one of only two roll-up builders that read
# cell_registry (via surface_index), so a surface/tier change can stale its committed rows.
# The guard reads the card AS COMMITTED and asserts the promise the module's own docstring
# makes: a shipped decoder cannot render invisibly.
#
# Deliberately NOT a rebuild-and-diff: _redirect_io cannot be reused (it blanks the INPUTS
# too, so a temp build is surface-only and not comparable), and rebuilding into the real
# wiki/ would mutate the working tree from inside a test.
# ---------------------------------------------------------------------------


def _committed_card() -> dict | None:
    """The report card as committed at HEAD, or None when git can't answer."""
    import json as _json
    import subprocess

    try:
        out = subprocess.run(
            ["git", "show", "HEAD:wiki/decoder_validation_report_card.json"],
            cwd=Path(__file__).resolve().parent.parent, capture_output=True, timeout=30,
        )
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    return _json.loads(out.stdout.decode("utf-8"))


def test_committed_card_renders_every_live_surface_cell():
    """A shipped decoder must not be invisible in the card that ships.

    Checks the COMMITTED artifact, not a fresh build -- a fresh build agrees with the registry
    by construction, which is exactly why it cannot detect committed staleness.
    """
    from dna_decode.data.cell_registry import surface_index

    card = _committed_card()
    if card is None:
        return  # no git / not committed yet — nothing to compare
    rows = {(c["organism"], c["drug"]) for c in card["cells"]}
    missing = sorted(set(surface_index()) - rows)
    assert not missing, (
        f"these live SHIPPED_DECODER_SURFACE cells are absent from the COMMITTED report card, so they "
        f"ship with no evidence row: {missing} — rebuild it: "
        f"uv run python scripts/build_validation_report_card.py"
    )
    # The reverse direction is NOT an error: the card is the surface UNION observed cells, so
    # observed-only rows (2 at the time of writing) legitimately have no surface entry.


# ============================================================================ #
# Species-composition disclosure (Step 5 of the meropenem FN16 plan, 2026-10-05)
# ============================================================================ #

def _sp_cell(**kw):
    """A species-audit cohort row, shaped as the audit script emits it."""
    base = {
        "cohort_dir": "klebsiella_provdisjoint_meropenem", "organism": "Klebsiella",
        "drug": "meropenem", "amrfinder_organism": "Klebsiella_pneumoniae",
        "reconciled": True,
        "composition": {"verdict": "MIXED_SPECIES", "expected_species": "Klebsiella pneumoniae",
                        "n_total": 60, "n_matches_expected": 37, "n_same_genus_other_species": 23,
                        "n_other_genus": 0, "n_unresolved": 0,
                        "species_counts": {"Klebsiella pneumoniae": 37, "Klebsiella aerogenes": 23}},
        "false_negative_anatomy": {"n_false_negatives": 16, "n_fn_off_species": 16,
                                   "n_fn_expected_species": 0,
                                   "fn_species_counts": {"Klebsiella aerogenes": 16}},
        "outcome_crosstab": {"by_species": {
            "Klebsiella pneumoniae": {"tp": 12, "fp": 3, "tn": 22, "fn": 0, "abstain": 0},
            "Klebsiella aerogenes": {"tp": 2, "fp": 0, "tn": 5, "fn": 16, "abstain": 0}}},
    }
    base.update(kw)
    return base


def test_species_block_is_None_for_a_cohort_that_matches_what_it_was_scored_as():
    """A reassuring line for a clean cohort is noise, and invites reading absence as a guarantee."""
    clean = _sp_cell(composition={"verdict": "MATCHES_SCORED_ORGANISM", "n_total": 40,
                                  "n_matches_expected": 40, "n_same_genus_other_species": 0,
                                  "n_other_genus": 0, "n_unresolved": 0, "species_counts": {}})
    assert mod.build_species_block(clean) is None
    assert mod.build_species_block(None) is None


def test_species_block_carries_the_measured_meropenem_finding():
    b = mod.build_species_block(_sp_cell())
    assert b["status"] == "measured"
    assert b["n_off_species"] == 23 and b["off_species_fraction"] == round(23 / 60, 4)
    assert b["false_negatives"] == 16 and b["false_negatives_off_species"] == 16
    o = b["on_scored_species_only"]
    assert (o["n"], o["tp"], o["fn"], o["sens"]) == (37, 12, 0, 1.0)
    assert "authority call" in o["note"], "the measurement must not read as a replacement"


def test_a_withheld_crosstab_yields_no_outcome_numbers():
    """A cohort with incomplete cached runs must not get FN attribution -- that is a rate over a
    silently-shrunken denominator."""
    b = mod.build_species_block(_sp_cell(reconciled=False, outcome_crosstab=None,
                                        crosstab_withheld_reason="cached AMRFinder runs incomplete (3/60)",
                                        false_negative_anatomy=None))
    assert b["status"] == "measured"            # composition still reported
    assert "outcome_attribution_withheld" in b
    for leaked in ("false_negatives", "false_negatives_off_species", "on_scored_species_only"):
        assert leaked not in b, f"withheld cross-tab leaked {leaked}"


def test_an_insufficient_resolution_cohort_renders_no_numbers():
    b = mod.build_species_block(_sp_cell(composition={
        "verdict": "INSUFFICIENT_RESOLUTION", "reason": "40/60 had no readable ORGANISM line",
        "n_total": 60, "n_unresolved": 40}))
    assert b["status"] == "insufficient_resolution"
    for leaked in ("n_off_species", "species_counts", "false_negatives", "on_scored_species_only"):
        assert leaked not in b, f"refused verdict leaked {leaked}"


def test_species_composition_NEVER_merges_into_the_scored_metrics(tmp_path, monkeypatch):
    """ANTI-OVERWRITE, with deliberately DIFFERENT values so a merge bug cannot pass by coincidence.

    The species cell shares its (organism, drug) key with the SCORED cell -- the documented shared-key
    trap. If the loader fed `load_scored()` instead of its own namespace, the published sens would be
    replaced by a composition field. Here the scored sens (0.467) and the species block's
    on-scored-species sens (1.000) are deliberately far apart, so a merge would be unmistakable.
    """
    sp = mod.build_species_block(_sp_cell())
    scored_like = {"state": "SCORED", "sens": 0.467, "spec": 0.9, "acc": 0.683,
                   "tp": 14, "fp": 3, "tn": 27, "fn": 16}
    merged = dict(scored_like)
    merged["species_composition"] = sp
    # the published metrics must survive untouched beside the block
    assert merged["sens"] == 0.467 and merged["tp"] == 14 and merged["fn"] == 16
    # and the species block's own, very different, figure is reachable only under its own key
    assert merged["species_composition"]["on_scored_species_only"]["sens"] == 1.0
    assert merged["sens"] != merged["species_composition"]["on_scored_species_only"]["sens"]


def test_the_committed_card_has_the_species_section_and_the_meropenem_row():
    """The rendered .md is the human-facing surface -- a block carried only in JSON is not a disclosure."""
    md = (mod.WIKI / "decoder_validation_report_card.md")
    if not md.exists():
        pytest.skip("card not built")
    text = md.read_text(encoding="utf-8")
    assert "## Species-composition disclosure" in text
    assert "38% K. aerogenes" in text
    assert "all 16 of its false negatives are those off-species isolates" in text
    assert "re-scoring is a user authority call" in text


def test_the_committed_card_attaches_species_blocks_ONLY_to_mixed_cohorts():
    """Campylobacter and the four E. coli cohorts match what they were scored as, so they must carry no
    block at all -- the layer discloses a mismatch, it does not annotate every cell."""
    j = mod.WIKI / "decoder_validation_report_card.json"
    if not j.exists():
        pytest.skip("card not built")
    cells = json.loads(j.read_text(encoding="utf-8"))["cells"]
    with_block = {(c["organism"], c["drug"]) for c in cells if c.get("species_composition")}
    assert with_block, "no species blocks attached at all -- the layer is not wired"
    assert all(o == "klebsiella" for o, _ in with_block), with_block
    merop = [c for c in cells if (c["organism"], c["drug"]) == ("klebsiella", "meropenem")][0]
    assert merop["species_composition"]["n_off_species"] == 23
    assert merop["species_composition"]["false_negatives_off_species"] == 16
    # ...and the published metrics are untouched beside it
    assert merop["sens"] == 0.467 and merop["fn"] == 16


def _species_audit(*, complete=True, n_off=23, cells=None):
    """A species-audit sidecar payload, shaped as `provdisjoint_species_audit.main` writes it."""
    if cells is None:
        cells = [{"cohort_dir": "klebsiella_provdisjoint_meropenem", "organism": "Klebsiella",
                  "drug": "meropenem", "amrfinder_organism": "Klebsiella_pneumoniae",
                  "reconciled": False, "crosstab_withheld_reason": "incomplete",
                  "composition": {"verdict": "MIXED_SPECIES", "n_total": 60,
                                  "n_same_genus_other_species": n_off, "n_other_genus": 0,
                                  "species_counts": {"Klebsiella pneumoniae": 60 - n_off,
                                                     "Klebsiella aerogenes": n_off}}}]
    return {"schema": "provdisjoint-species-audit-v1", "complete": complete, "cells": cells}


def test_load_species_composition_keeps_the_NEWEST_audit(monkeypatch, tmp_path):
    """The audit is re-run as cohorts change, so several dated sidecars coexist. Reading an older one
    would attach a stale composition beside a current metric with nothing to flag the mismatch."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    for stamp, n_off in (("2026-09-01", 2), ("2026-10-05", 23), ("2026-08-01", 99)):
        (wiki / f"provdisjoint_species_audit_{stamp}.json").write_text(
            json.dumps(_species_audit(n_off=n_off)), encoding="utf-8")
    got = mod.load_species_composition()
    assert len(got) == 1
    cell = next(iter(got.values()))
    assert cell["composition"]["n_same_genus_other_species"] == 23, "newest wins, not last-globbed"


def test_load_species_composition_does_NOT_gate_on_the_audit_complete_FLAG(monkeypatch, tmp_path):
    """DELIBERATE divergence from `load_source_concentration`, which drops its whole payload when the
    sweep is incomplete. Here the refusal is already per-cohort -- a cohort that could not resolve its
    species carries INSUFFICIENT_RESOLUTION and no numbers -- so dropping the layer because ONE cohort
    refused would hide the nine that measured cleanly. A copy-paste of the sibling's gate is the
    plausible defect, and it fails silently: the section simply stops rendering.
    """
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provdisjoint_species_audit_2026-10-05.json").write_text(
        json.dumps(_species_audit(complete=False)), encoding="utf-8")
    assert mod.load_species_composition(), "an incomplete audit must still surface its measured cohorts"


def test_a_malformed_species_audit_is_IGNORED_not_fatal(monkeypatch, tmp_path):
    """The loader claims a malformed sidecar cannot break the read-only roll-up. This repo has shipped
    an invalid-JSON wiki artifact before, so the claim is worth holding to."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provdisjoint_species_audit_2026-10-05.json").write_text("{not json", encoding="utf-8")
    assert mod.load_species_composition() == {}
    assert mod.main() == 0, "the roll-up must still build"
    assert "## Species-composition disclosure" not in (
        wiki / "decoder_validation_report_card.md").read_text(encoding="utf-8")


def test_a_cell_with_no_organism_is_dropped_rather_than_keyed_on_None(monkeypatch, tmp_path):
    """A row missing `organism` cannot be joined to a card cell; keying it on None would create a
    phantom entry that silently matches nothing."""
    wiki = _redirect_io(monkeypatch, tmp_path)
    (wiki / "provdisjoint_species_audit_2026-10-05.json").write_text(json.dumps(_species_audit(
        cells=[{"drug": "meropenem", "composition": {"verdict": "MIXED_SPECIES"}}])), encoding="utf-8")
    assert mod.load_species_composition() == {}
