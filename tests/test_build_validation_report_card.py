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


def _prospective_artifact(org="Klebsiella", drug="ciprofloxacin", *, acc=0.60, sens=0.40,
                          spec=0.90, n=50, generated="2026-08-24", verified=True):
    return {
        "artifact": "prospective_lock_validation", "generated": generated,
        "organism": org, "drug": drug, "prospective_lock_verified": verified,
        "lock_manifest": {"lock_date": "2026-06-13"},
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
        "lock_manifest": {"lock_date": "2026-06-13"},
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
