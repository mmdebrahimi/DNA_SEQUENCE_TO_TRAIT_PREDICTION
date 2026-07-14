"""Tests for the g8 residual-signal detector (pure, offline, CI-safe)."""
from __future__ import annotations

import pytest

from dna_decode.eval import residual_detector as R


def test_classify_residual_three_tiers():
    assert R.classify_residual(None) == R.UNTESTED
    assert R.classify_residual({"generalizes_beyond_lineage": True, "clade_concentrated": False}) == R.GENERALIZES
    assert R.classify_residual({"generalizes_beyond_lineage": False, "clade_concentrated": True}) == R.LINEAGE_MEDIATED
    # tested but neither flag -> conservative LINEAGE_MEDIATED (signal did not survive de-confounding)
    assert R.classify_residual({"generalizes_beyond_lineage": False, "clade_concentrated": False}) == R.LINEAGE_MEDIATED
    # clade_concentrated wins even if generalizes flag is somehow also set
    assert R.classify_residual({"generalizes_beyond_lineage": True, "clade_concentrated": True}) == R.LINEAGE_MEDIATED


def _artifact():
    return {
        "axis_label": "AMR determinant", "organism": "escherichia_coli_shigella",
        "verdict": "CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE",
        "median_auc_naive": 0.975, "median_auc_clade_grouped": 0.908,
        "mash": {"threshold": 0.005, "n_clades": 80, "largest_clade_frac": 0.353},
        "honest_caveats": ["source caveat X"],
        "per_gene": {
            "gyrA_S83L": {"n_present": 200, "auc_naive": 0.96, "auc_clade_grouped": 0.90, "drop": 0.06,
                          "generalizes_beyond_lineage": True, "clade_concentrated": False},
            "parC_E84V": {"n_present": 77, "auc_naive": 0.89, "auc_clade_grouped": 0.60, "drop": 0.29,
                          "generalizes_beyond_lineage": False, "clade_concentrated": True},
            "blaCTX-M-15": {"n_present": 90, "auc_naive": 0.80, "auc_clade_grouped": 0.72, "drop": 0.08,
                            "generalizes_beyond_lineage": True, "clade_concentrated": False},
        },
    }


def test_build_report_tier_counts_and_order():
    rep = R.build_residual_report(_artifact())
    assert rep["tier_counts"] == {R.GENERALIZES: 2, R.LINEAGE_MEDIATED: 1, R.UNTESTED: 0}
    # GENERALIZES first, strongest de-confounded (clade-grouped) AUC on top
    ids = [r["feature_id"] for r in rep["per_feature"]]
    assert ids[0] == "gyrA_S83L"      # generalizes, clade-grouped 0.90 (highest)
    assert ids[1] == "blaCTX-M-15"    # generalizes, clade-grouped 0.72
    assert ids[2] == "parC_E84V"      # lineage-mediated last
    assert rep["per_feature"][2]["tier"] == R.LINEAGE_MEDIATED


def test_family_rollup_and_meta():
    rep = R.build_residual_report(_artifact())
    fr = rep["family_rollup"]
    assert fr["QRDR"][R.GENERALIZES] == 1 and fr["QRDR"][R.LINEAGE_MEDIATED] == 1   # gyrA gen, parC lineage
    assert fr["beta-lactam"][R.GENERALIZES] == 1                                    # blaCTX-M
    assert rep["meta"]["mash_threshold"] == 0.005 and rep["meta"]["n_features"] == 3
    assert rep["meta"]["source_verdict"] == "CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE"


def test_honest_caveats_present_and_no_phenotype_claim():
    rep = R.build_residual_report(_artifact())
    joined = " ".join(rep["honest_caveats"]).lower()
    assert "not a phenotype prediction" in joined and "clonal" in joined
    assert "source caveat x" in joined       # source caveats appended
    # the report never emits an R/S field
    for r in rep["per_feature"]:
        assert "prediction" not in r and "R/S" not in str(r)


def test_viz_adapter_reuses_canonical_classification():
    # the network_adapter must re-export the SAME classification (no drift between the viz + the product)
    from dna_decode.viz import network_adapter as A
    assert A._lineage_status is R.classify_residual
    assert (A.GENERALIZES, A.LINEAGE_MEDIATED, A.UNTESTED) == (R.GENERALIZES, R.LINEAGE_MEDIATED, R.UNTESTED)


def test_empty_artifact_is_safe():
    rep = R.build_residual_report({})
    assert rep["tier_counts"] == {R.GENERALIZES: 0, R.LINEAGE_MEDIATED: 0, R.UNTESTED: 0}
    assert rep["per_feature"] == [] and rep["meta"]["n_features"] == 0
