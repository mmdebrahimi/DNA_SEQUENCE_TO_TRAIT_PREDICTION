"""Unit tests for the external-validation roll-up (no Docker — synthetic clusters)."""
from __future__ import annotations

import json

import scripts.build_external_validation_report as rep


# --------------------------------------------------------------------------- #
# load_external_artifacts
# --------------------------------------------------------------------------- #
def test_load_external_artifacts(tmp_path):
    (tmp_path / "external_validation_spain_ciprofloxacin_2026-06-15.json").write_text(
        json.dumps({"_schema": "external-validation-v1", "cohort": "spain", "drug": "ciprofloxacin"}))
    (tmp_path / "external_validation_report_card.json").write_text(
        json.dumps({"_schema": "external-validation-report-card-v1"}))  # the roll-up itself -> excluded
    (tmp_path / "provenance_disjoint_validation_x.json").write_text(
        json.dumps({"_schema": "provenance-disjoint-validation-v1"}))   # frozen surface -> not globbed
    arts = rep.load_external_artifacts(tmp_path, allow_unscoped_glob=True)
    assert len(arts) == 1
    assert arts[0]["cohort"] == "spain"


def test_load_ignores_wrong_schema(tmp_path):
    (tmp_path / "external_validation_bad.json").write_text(json.dumps({"_schema": "something-else"}))
    assert rep.load_external_artifacts(tmp_path, allow_unscoped_glob=True) == []


# --- Step 9: run-scoping + hard_fail/degraded filtering ---
def test_load_refuses_glob_all_by_default(tmp_path):
    import pytest
    (tmp_path / "external_validation_x_cipro_run1_2026.json").write_text(
        json.dumps({"_schema": "external-validation-v1", "cohort": "x"}))
    with pytest.raises(ValueError):
        rep.load_external_artifacts(tmp_path)   # no run_id / artifacts / unscoped flag


def test_load_run_scoped(tmp_path):
    (tmp_path / "external_validation_x_cipro_run1_2026.json").write_text(
        json.dumps({"_schema": "external-validation-v1", "cohort": "x", "run_id": "run1"}))
    (tmp_path / "external_validation_x_cipro_run2_2026.json").write_text(
        json.dumps({"_schema": "external-validation-v1", "cohort": "x", "run_id": "run2"}))
    arts = rep.load_external_artifacts(tmp_path, run_id="run1")
    assert [a["run_id"] for a in arts] == ["run1"]   # stale run2 excluded


def test_load_skips_hard_fail_and_degraded(tmp_path):
    (tmp_path / "external_validation_a_cipro_run1_2026.json").write_text(json.dumps(
        {"_schema": "external-validation-v1", "cohort": "a", "powering": {"hard_fail": True}}))
    (tmp_path / "external_validation_b_cipro_run1_2026.json").write_text(json.dumps(
        {"_schema": "external-validation-v1", "cohort": "b", "run_degraded": True}))
    (tmp_path / "external_validation_c_cipro_run1_2026.json").write_text(json.dumps(
        {"_schema": "external-validation-v1", "cohort": "c"}))
    clean = rep.load_external_artifacts(tmp_path, run_id="run1")
    assert {a["cohort"] for a in clean} == {"c"}          # hard_fail + degraded skipped
    with_deg = rep.load_external_artifacts(tmp_path, run_id="run1", allow_degraded=True)
    assert {a["cohort"] for a in with_deg} == {"b", "c"}  # degraded now included, hard_fail still out


def test_load_skips_corrupt_json(tmp_path):
    # A truncated/corrupt artifact is skipped, not raised; a good one alongside still loads.
    (tmp_path / "external_validation_corrupt.json").write_text("{not valid json")
    (tmp_path / "external_validation_ok_ciprofloxacin_2026-06-15.json").write_text(
        json.dumps({"_schema": "external-validation-v1", "cohort": "ok", "drug": "ciprofloxacin"}))
    arts = rep.load_external_artifacts(tmp_path, allow_unscoped_glob=True)
    assert [a["cohort"] for a in arts] == ["ok"]


# --- gap: _publishable predicate (tested directly, not only via load) ------ #
def test_publishable_clean_cell():
    assert rep._publishable({}, allow_degraded=False) is True


def test_publishable_hard_fail_never_published():
    assert rep._publishable({"powering": {"hard_fail": True}}, allow_degraded=True) is False


def test_publishable_degraded_gated_by_flag():
    assert rep._publishable({"run_degraded": True}, allow_degraded=False) is False
    assert rep._publishable({"run_degraded": True}, allow_degraded=True) is True


# --- gap: explicit `artifacts` list takes precedence over run_id/glob ------- #
def test_load_explicit_artifacts_precedence(tmp_path):
    p = tmp_path / "external_validation_z_cipro_runZ_2026.json"
    p.write_text(json.dumps({"_schema": "external-validation-v1", "cohort": "z"}))
    arts = rep.load_external_artifacts(tmp_path, artifacts=[str(p)])
    assert [a["cohort"] for a in arts] == ["z"]
    assert arts[0]["_path"] == str(p)


# --------------------------------------------------------------------------- #
# _fmt_ci fallback (malformed / missing CI -> em-dash, never a crash)
# --------------------------------------------------------------------------- #
def test_fmt_ci_wellformed():
    assert rep._fmt_ci((0.5, 0.95)) == "[0.5, 0.95]"


def test_fmt_ci_malformed_falls_back():
    assert rep._fmt_ci(None) == "—"
    assert rep._fmt_ci((0.5,)) == "—"            # wrong length
    assert rep._fmt_ci("nonsense") == "—"        # wrong type


# --------------------------------------------------------------------------- #
# cluster_weighted_with_ci (synthetic clusters; reuses clonality math)
# --------------------------------------------------------------------------- #
def test_cluster_weighted_with_ci():
    # 3 clusters: one R-clone (correct), one S-clone (correct), one DISCORDANT.
    clusters = {"g1": 0, "g2": 0, "g3": 1, "g4": 2, "g5": 2}
    labels = {"g1": "R", "g2": "R", "g3": "S", "g4": "R", "g5": "S"}  # cluster 2 mixed -> DISCORDANT
    preds = {"g1": "R", "g2": "R", "g3": "S", "g4": "R", "g5": "S"}
    out = rep.cluster_weighted_with_ci(preds, labels, clusters)
    assert out["n_discordant"] == 1                 # cluster 2 excluded
    assert out["tp"] == 1 and out["tn"] == 1        # one R-clone, one S-clone
    assert isinstance(out["sens_ci"], tuple) and len(out["sens_ci"]) == 2
    assert out["effective_lineage_n_R"] == 1 and out["effective_lineage_n_S"] == 1


def test_cluster_weighted_with_ci_one_class_only():
    # All-R clusters: spec denominator is 0 -> spec_ci spans (0.0, 1.0), no crash;
    # eff-N S is 0. Guards the zero-denominator Wilson path the report card relies on.
    clusters = {"g1": 0, "g2": 1}
    labels = {"g1": "R", "g2": "R"}
    preds = {"g1": "R", "g2": "R"}
    out = rep.cluster_weighted_with_ci(preds, labels, clusters)
    assert out["tp"] == 2 and out["fp"] == 0 and out["tn"] == 0 and out["fn"] == 0
    assert out["spec_ci"] == (0.0, 1.0)
    assert out["effective_lineage_n_S"] == 0 and out["effective_lineage_n_R"] == 2


# --------------------------------------------------------------------------- #
# build_cell + rendering
# --------------------------------------------------------------------------- #
def _artifact():
    return {
        "cohort": "spain_probac", "organism": "Escherichia_coli_Shigella", "drug": "ciprofloxacin",
        "evidence_tier": "external_clinical",
        "strict": {"n_scored": 20, "acc": 0.9, "sens": 0.91, "spec": 0.89},
        "relaxed": {"n_scored": 28, "acc": 0.86, "sens": 0.87, "spec": 0.85},
        "independence_tier": "external clinical cohort ...",
    }


def test_build_cell_with_lineage():
    cell = rep.build_cell(_artifact(), {"status": "ok", "sens": 0.8, "spec": 0.9,
                                        "sens_ci": (0.5, 0.95), "spec_ci": (0.6, 0.99),
                                        "effective_lineage_n_R": 3, "effective_lineage_n_S": 4})
    assert cell["cohort"] == "spain_probac"
    assert cell["lineage"]["status"] == "ok"


def test_build_cell_lineage_unavailable_default():
    cell = rep.build_cell(_artifact(), None)
    assert cell["lineage"]["status"] == "unavailable"


def test_render_json_and_md():
    cell = rep.build_cell(_artifact(), {"status": "ok", "sens": 0.8, "spec": 0.9,
                                        "sens_ci": (0.5, 0.95), "spec_ci": (0.6, 0.99),
                                        "effective_lineage_n_R": 3, "effective_lineage_n_S": 4})
    j = rep.render_json([cell])
    assert j["_schema"] == "external-validation-report-card-v1"
    assert j["n_cells"] == 1
    md = rep.render_md([cell])
    assert "External-validation report card" in md
    assert "spain_probac" in md
    assert "[0.5, 0.95]" in md                       # Wilson CI rendered


def test_render_md_unavailable_lineage():
    md = rep.render_md([rep.build_cell(_artifact(), None)])
    assert "n/a" in md                               # degraded lineage shows n/a, raw still rendered
    assert "0.91" in md                              # strict raw sens still present


# --------------------------------------------------------------------------- #
# Fix C: our output namespace is DISTINCT from the frozen report card
# --------------------------------------------------------------------------- #
def test_output_namespace_distinct_from_frozen():
    assert rep.OUR_REPORT_CARD != rep.FROZEN_REPORT_CARD
    assert "external" in rep.OUR_REPORT_CARD


def test_main_writes_separate_namespace(tmp_path):
    (tmp_path / "external_validation_spain_ciprofloxacin_2026-06-15.json").write_text(
        json.dumps(_artifact() | {"_schema": "external-validation-v1"}))
    import sys
    argv = sys.argv
    sys.argv = ["prog", "--wiki-dir", str(tmp_path), "--no-clonality", "--allow-unscoped-glob"]
    try:
        assert rep.main() == 0
    finally:
        sys.argv = argv
    assert (tmp_path / "external_validation_report_card.json").exists()
    assert (tmp_path / "external_validation_report_card.md").exists()
    # the frozen decoder report card must NOT be created by our roll-up
    assert not (tmp_path / "decoder_validation_report_card.json").exists()


def test_build_cell_passes_through_experimental_branding():
    art = {"cohort": "sci234", "drug": "trimethoprim-sulfamethoxazole",
           "strict": {}, "binary": {"sens": 0.97, "spec": 0.99, "n_scored": 210},
           "rule_status": "EXPERIMENTAL_SCORED", "rule_scope": "scorer_local",
           "not_in_shipped_surface": True, "strata_reproduced": True}
    c = rep.build_cell(art, None)
    assert c["rule_status"] == "EXPERIMENTAL_SCORED"
    assert c["rule_scope"] == "scorer_local"
    assert c["not_in_shipped_surface"] is True
    assert c["binary"]["n_scored"] == 210


def test_build_cell_frozen_decoder_cell_unbranded():
    c = rep.build_cell({"cohort": "oxford", "drug": "ciprofloxacin",
                        "strict": {"sens": 0.94, "spec": 0.96, "n_scored": 1005}}, None)
    assert c["rule_status"] is None and c["rule_scope"] is None


def test_render_md_both_strict_and_binary_empty_defaults_to_strict():
    # When strict is not usable AND binary has no scored rows, the `not b.get("n_scored")`
    # branch keeps metric="strict" (all-None) rather than crashing or emitting "binary".
    cell = rep.build_cell({"cohort": "empty", "drug": "trimethoprim-sulfamethoxazole",
                           "strict": {}, "binary": {}}, None)
    md = rep.render_md([cell])
    assert "| strict |" in md          # falls back to strict label, not binary
    assert "| binary |" not in md
    assert "| empty |" in md           # row still renders (no crash on all-None metrics)


def test_render_md_marks_experimental_and_falls_back_to_binary():
    exp = rep.build_cell({"cohort": "sci234", "drug": "trimethoprim-sulfamethoxazole",
                          "strict": {}, "binary": {"sens": 0.97, "spec": 0.99, "n_scored": 210},
                          "rule_status": "EXPERIMENTAL_SCORED", "rule_scope": "scorer_local"}, None)
    dep = rep.build_cell({"cohort": "oxford", "drug": "ciprofloxacin",
                          "strict": {"sens": 0.94, "spec": 0.96, "n_scored": 1005}}, None)
    md = rep.render_md([exp, dep])
    assert "EXPERIMENTAL_SCORED (scorer_local)" in md   # experimental cell marked
    assert "deployed-decoder" in md                      # frozen-decoder cell marked
    assert "| binary |" in md                            # strict-degenerate cell falls back to binary
    assert "| strict |" in md                            # the deployed cell shows strict
