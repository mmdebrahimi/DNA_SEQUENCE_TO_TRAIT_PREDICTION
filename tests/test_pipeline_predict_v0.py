"""Tests for the v0 predict schema helpers in scripts/pipeline.py.

The v0 schema is defined in `wiki/decoder_v0_ux_and_success_criterion.md`.
Pure-function helpers are tested directly; the end-to-end XGBoost + ISM
path is orchestration (skipped — uses real cache + trained model).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# ---- _confidence_tier ----


def test_confidence_tier_high_above_threshold():
    from scripts.pipeline import _confidence_tier
    assert _confidence_tier(0.95) == "HIGH"
    assert _confidence_tier(0.9) == "HIGH"


def test_confidence_tier_high_below_threshold():
    from scripts.pipeline import _confidence_tier
    assert _confidence_tier(0.05) == "HIGH"
    assert _confidence_tier(0.1) == "HIGH"


def test_confidence_tier_medium():
    from scripts.pipeline import _confidence_tier
    assert _confidence_tier(0.8) == "MEDIUM"
    assert _confidence_tier(0.2) == "MEDIUM"
    assert _confidence_tier(0.7) == "MEDIUM"
    assert _confidence_tier(0.3) == "MEDIUM"


def test_confidence_tier_low_near_threshold():
    from scripts.pipeline import _confidence_tier
    assert _confidence_tier(0.6) == "LOW"
    assert _confidence_tier(0.5) == "LOW"
    assert _confidence_tier(0.4) == "LOW"


# ---- _load_audit_verdict ----


def test_load_audit_verdict_none_when_no_path():
    from scripts.pipeline import _load_audit_verdict
    assert _load_audit_verdict(None, "S1") is None


def test_load_audit_verdict_none_when_file_missing(tmp_path: Path):
    from scripts.pipeline import _load_audit_verdict
    assert _load_audit_verdict(tmp_path / "nope.json", "S1") is None


def test_load_audit_verdict_returns_strain_row(tmp_path: Path):
    from scripts.pipeline import _load_audit_verdict

    merge_json = tmp_path / "merge.json"
    merge_json.write_text(json.dumps({
        "gate_verdict": "SUSPEND_CONDITION_4",
        "per_strain": [
            {
                "strain_id": "S1",
                "noise_class": "OPAQUE_R_no_mechanism",
                "mic_tier": "HIGH_R",
                "mechanism_opacity_flag": True,
                "primary_mechanisms": [],
                "co_resistance_modifiers": ["efflux"],
            },
            {"strain_id": "S2", "noise_class": "CLEAN_R_primary_mechanism", "mic_tier": "HIGH_R"},
        ],
    }), encoding="utf-8")

    v = _load_audit_verdict(merge_json, "S1")
    assert v is not None
    assert v["noise_class"] == "OPAQUE_R_no_mechanism"
    assert v["mic_tier"] == "HIGH_R"
    assert v["mechanism_opacity_flag"] is True
    assert v["co_resistance_modifiers"] == ["efflux"]
    assert v["cohort_gate_verdict"] == "SUSPEND_CONDITION_4"
    # SUSPEND propagates explicitly
    assert v["suspend_gate_fired"] is True
    assert "informational only" in v["verdict_explanation"]


def test_load_audit_verdict_returns_none_for_strain_not_in_audit(tmp_path: Path):
    from scripts.pipeline import _load_audit_verdict

    merge_json = tmp_path / "merge.json"
    merge_json.write_text(json.dumps({
        "gate_verdict": "RUN_FULL_AND_CLEAN",
        "per_strain": [{"strain_id": "S1", "noise_class": "CLEAN_R_primary_mechanism"}],
    }), encoding="utf-8")

    assert _load_audit_verdict(merge_json, "S_UNKNOWN") is None


def test_load_audit_verdict_no_suspend_when_gate_clean(tmp_path: Path):
    from scripts.pipeline import _load_audit_verdict

    merge_json = tmp_path / "merge.json"
    merge_json.write_text(json.dumps({
        "gate_verdict": "RUN_FULL_AND_CLEAN",
        "per_strain": [{"strain_id": "S1", "noise_class": "CLEAN_R_primary_mechanism", "mic_tier": "HIGH_R"}],
    }), encoding="utf-8")

    v = _load_audit_verdict(merge_json, "S1")
    assert v["suspend_gate_fired"] is False
    assert "verdict_explanation" not in v


# ---- _render_predict_markdown ----


def _sample_result(with_attribution: bool = True, with_audit: bool = True) -> dict:
    return {
        "strain_id": "562.12345",
        "drug": "ciprofloxacin",
        "prediction": "R",
        "calibrated_probability": 0.87,
        "confidence_tier": "MEDIUM",
        "top_k_attribution": [
            {"gene_id": "gene-gyrA", "locus_tag": "b2231", "score": 0.42, "tier": "Tier 1 (textbook QRDR)"},
            {"gene_id": "gene-parC", "locus_tag": "b3019", "score": 0.18, "tier": "Tier 1 (textbook QRDR)"},
        ] if with_attribution else [],
        "audit_verdict": {
            "noise_class": "CLEAN_R_primary_mechanism",
            "mic_tier": "HIGH_R",
            "mechanism_opacity_flag": False,
            "primary_mechanisms": ["QRDR_target_alteration"],
            "co_resistance_modifiers": [],
            "cohort_gate_verdict": "RUN_FULL_AND_CLEAN",
            "suspend_gate_fired": False,
        } if with_audit else None,
        "provenance": {
            "model": "nucleotide_transformer + XGBoost (frozen)",
            "training_cohort": "stage2_n150_cipro_cohort",
            "loso_auroc": 0.78,
            "lomo_clade_out_auroc": None,
            "trained_on": "2026-05-18",
        },
    }


def test_markdown_includes_prediction_header_and_probability():
    from scripts.pipeline import _render_predict_markdown
    md = _render_predict_markdown(_sample_result())
    assert "strain `562.12345`" in md
    assert "drug `ciprofloxacin`" in md
    assert "**Prediction:** R" in md
    assert "0.870" in md
    assert "**Confidence tier:** MEDIUM" in md


def test_markdown_renders_attribution_table():
    from scripts.pipeline import _render_predict_markdown
    md = _render_predict_markdown(_sample_result())
    assert "Top-K gene attribution" in md
    assert "gene-gyrA" in md
    assert "Tier 1 (textbook QRDR)" in md
    # Both rows present
    assert "gene-parC" in md


def test_markdown_handles_no_attribution():
    from scripts.pipeline import _render_predict_markdown
    md = _render_predict_markdown(_sample_result(with_attribution=False))
    assert "no attribution run" in md
    assert "gene-gyrA" not in md


def test_markdown_handles_no_audit_verdict():
    from scripts.pipeline import _render_predict_markdown
    md = _render_predict_markdown(_sample_result(with_audit=False))
    assert "no audit data" in md


def test_markdown_surfaces_suspend_warning():
    from scripts.pipeline import _render_predict_markdown
    r = _sample_result()
    r["audit_verdict"] = {
        "noise_class": "OPAQUE_R_no_mechanism",
        "mic_tier": "HIGH_R",
        "mechanism_opacity_flag": True,
        "primary_mechanisms": [],
        "co_resistance_modifiers": ["efflux"],
        "cohort_gate_verdict": "SUSPEND_CONDITION_4",
        "suspend_gate_fired": True,
        "verdict_explanation": "Training cohort fired SUSPEND_CONDITION_4; categorical labels carry noise. "
                               "Prediction is informational only; do not deploy as clinical decision support.",
    }
    md = _render_predict_markdown(r)
    assert "**SUSPEND gate fired on training cohort.**" in md
    assert "informational only" in md
    assert "do not deploy as clinical decision support" in md


def test_markdown_includes_provenance_block():
    from scripts.pipeline import _render_predict_markdown
    md = _render_predict_markdown(_sample_result())
    assert "## Provenance" in md
    assert "nucleotide_transformer + XGBoost (frozen)" in md
    assert "stage2_n150_cipro_cohort" in md
    assert "0.78" in md


def test_markdown_disclaimer_present():
    """v0 success criterion #4: honest output — no overclaiming."""
    from scripts.pipeline import _render_predict_markdown
    md = _render_predict_markdown(_sample_result())
    assert "Not a clinical decision support tool" in md


# ---- _locus_tag_prefix ----


def test_locus_tag_prefix_extracts_alpha_prefix():
    from scripts.pipeline import _locus_tag_prefix
    assert _locus_tag_prefix("ELX_001") == "ELX"
    assert _locus_tag_prefix("ERS123456_00010") == "ERS"
    assert _locus_tag_prefix("b2231") == "B"  # E. coli K12-style: alpha prefix only
    assert _locus_tag_prefix("OLZ10_24785") == "OLZ"


def test_locus_tag_prefix_handles_edge_cases():
    from scripts.pipeline import _locus_tag_prefix
    assert _locus_tag_prefix("") == ""
    assert _locus_tag_prefix("123_no_alpha_first") == ""
    assert _locus_tag_prefix(None) == ""  # tolerant


# ---- _classify_attribution_scope ----


def test_attribution_scope_indeterminate_pre_falsifier():
    """Pre-falsifier (verdict=None) always returns INDETERMINATE."""
    from scripts.pipeline import _classify_attribution_scope
    assert _classify_attribution_scope("ERS", falsifier_verdict=None) == "INDETERMINATE"
    assert _classify_attribution_scope("ELX", falsifier_verdict=None) == "INDETERMINATE"
    assert _classify_attribution_scope("", falsifier_verdict=None) == "INDETERMINATE"


def test_attribution_scope_indeterminate_when_saturated():
    """Saturated classifier -> INDETERMINATE regardless of locus-tag prefix."""
    from scripts.pipeline import _classify_attribution_scope
    assert _classify_attribution_scope("ERS", saturated=True, falsifier_verdict="PASS") == "INDETERMINATE"
    assert _classify_attribution_scope("ELX", saturated=True, falsifier_verdict="PASS") == "INDETERMINATE"


def test_attribution_scope_indeterminate_when_all_negative_delta():
    """Bucket C falsifier pattern -> INDETERMINATE."""
    from scripts.pipeline import _classify_attribution_scope
    result = _classify_attribution_scope(
        "OLZ", all_negative_delta=True, falsifier_verdict="PASS"
    )
    assert result == "INDETERMINATE"


def test_attribution_scope_high_for_ers_post_pass():
    """ERS-prefix strain after PASS verdict -> HIGH."""
    from scripts.pipeline import _classify_attribution_scope
    assert _classify_attribution_scope("ERS", falsifier_verdict="PASS") == "HIGH"
    assert _classify_attribution_scope("ers", falsifier_verdict="PASS") == "HIGH"  # case insensitive


def test_attribution_scope_partial_for_elx_family_post_pass():
    """ELX-family strains -> PARTIAL even on PASS (batch-clade failure case)."""
    from scripts.pipeline import _classify_attribution_scope
    for prefix in ("ELX", "ELY", "ELV", "ELU", "ELT"):
        assert _classify_attribution_scope(prefix, falsifier_verdict="PASS") == "PARTIAL"


def test_attribution_scope_partial_fallback_for_unknown_prefix_post_pass():
    """Unknown prefix on PASS -> PARTIAL (no falsifier evidence)."""
    from scripts.pipeline import _classify_attribution_scope
    assert _classify_attribution_scope("XYZ", falsifier_verdict="PASS") == "PARTIAL"
    assert _classify_attribution_scope("", falsifier_verdict="PASS") == "PARTIAL"


def test_attribution_scope_partial_for_ers_post_fail():
    """ERS on FAIL -> PARTIAL (ERS-control passed but cohort-wide method failed)."""
    from scripts.pipeline import _classify_attribution_scope
    assert _classify_attribution_scope("ERS", falsifier_verdict="FAIL") == "PARTIAL"


def test_attribution_scope_high_via_cohort_wide_override():
    """Explicit cohort-wide HIGH flag -> HIGH for any non-INDETERMINATE strain."""
    from scripts.pipeline import _classify_attribution_scope
    assert _classify_attribution_scope(
        "XYZ", falsifier_verdict="PASS", falsifier_pass_passes_high=True
    ) == "HIGH"
    # But saturation/negative-delta still wins
    assert _classify_attribution_scope(
        "XYZ", saturated=True, falsifier_verdict="PASS", falsifier_pass_passes_high=True
    ) == "INDETERMINATE"


# ---- attribution scope confidence rendering ----


def test_markdown_includes_attribution_scope_confidence():
    """Field shows up in markdown sidecar."""
    from scripts.pipeline import _render_predict_markdown
    r = _sample_result()
    r["attribution_scope_confidence"] = "PARTIAL"
    md = _render_predict_markdown(r)
    assert "**Attribution scope confidence:** PARTIAL" in md


def test_markdown_defaults_attribution_scope_to_indeterminate_when_missing():
    """Backward-compat: result dict without the field renders as INDETERMINATE."""
    from scripts.pipeline import _render_predict_markdown
    r = _sample_result()
    r.pop("attribution_scope_confidence", None)
    md = _render_predict_markdown(r)
    assert "**Attribution scope confidence:** INDETERMINATE" in md


# ---- RELOCKED 2026-05-23 v0 spec provenance fields ----


def test_markdown_renders_cv_strategy_and_auroc_when_present():
    """RELOCKED spec adds cv_strategy + cv_auroc to provenance; both surface in markdown."""
    from scripts.pipeline import _render_predict_markdown
    r = _sample_result()
    r["provenance"] = {
        "model": "nucleotide_transformer + XGBoost",
        "training_cohort": "stage2_n150_cipro_cohort",
        "cv_strategy": "leave_one_accession_out",
        "cv_auroc": 0.8697,
        "reporting_mode": "canonical_audit_aware",
        "trained_on": "2026-05-22",
    }
    md = _render_predict_markdown(r)
    assert "CV strategy: leave_one_accession_out" in md
    assert "CV AUROC: 0.8697" in md
    assert "Reporting mode: canonical_audit_aware" in md


def test_markdown_backward_compat_loso_auroc_only():
    """Older bundles only have loso_auroc; markdown still renders cleanly."""
    from scripts.pipeline import _render_predict_markdown
    r = _sample_result()
    r["provenance"] = {
        "model": "nucleotide_transformer + XGBoost (frozen)",
        "training_cohort": "stage2_n150_cipro_cohort",
        "loso_auroc": 0.78,
        "trained_on": "2026-05-18",
    }
    md = _render_predict_markdown(r)
    assert "LOSO AUROC (legacy field): 0.78" in md
    # No phantom cv_strategy line when fields are absent
    assert "CV strategy:" not in md


def test_markdown_renders_both_cv_and_legacy_loso_when_present():
    """Backward-compat: bundles may carry BOTH cv_auroc (canonical) AND loso_auroc (legacy)."""
    from scripts.pipeline import _render_predict_markdown
    r = _sample_result()
    r["provenance"] = {
        "model": "nucleotide_transformer + XGBoost",
        "training_cohort": "stage2_n150_cipro_cohort",
        "cv_strategy": "leave_one_accession_out",
        "cv_auroc": 0.87,
        "loso_auroc": 0.78,
        "trained_on": "2026-05-22",
    }
    md = _render_predict_markdown(r)
    assert "CV strategy: leave_one_accession_out" in md
    assert "CV AUROC: 0.87" in md
    assert "LOSO AUROC (legacy field): 0.78" in md
