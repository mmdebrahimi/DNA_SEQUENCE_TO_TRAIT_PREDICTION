"""Contract regression tests for scripts/drug_mechanism_phenotype_merge.py.

Pins the JSON contract consumed by scripts/pipeline.py::_load_audit_verdict and
the two locked parameters from plans/Cef_Audit_Aware_Packet_Design.md:

- signal_quality formula:  clean_count / max(1, len(merged))   (matches cipro)
- default --suspend-threshold:  0.40                            (matches cipro)

These tests exist because Codex's 2026-05-26 bundle shipped a version of the
merge script with both parameters wrong (clean / (clean+suspect+opacity); 0.50).
A future bundle drift would silently misfire the SUSPEND gate at cef closeout
without this regression coverage.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

from scripts import drug_mechanism_phenotype_merge as merge_mod
from scripts.pipeline import _load_audit_verdict


# ---------------------------------------------------------------------------
# Synthetic-fixture helpers
# ---------------------------------------------------------------------------


def _build_mech_audit(per_strain: list[dict]) -> dict:
    return {"drug": "ceftriaxone", "per_strain": per_strain}


def _build_mic_audit(per_strain: list[dict]) -> dict:
    return {"drug": "ceftriaxone", "per_strain": per_strain}


def _write_audits(
    tmp_path: Path,
    mech_rows: list[dict],
    mic_rows: list[dict],
) -> tuple[Path, Path]:
    mech_p = tmp_path / "mech.json"
    mic_p = tmp_path / "mic.json"
    mech_p.write_text(json.dumps(_build_mech_audit(mech_rows)), encoding="utf-8")
    mic_p.write_text(json.dumps(_build_mic_audit(mic_rows)), encoding="utf-8")
    return mech_p, mic_p


def _run_merge(
    tmp_path: Path,
    mech_rows: list[dict],
    mic_rows: list[dict],
    extra_args: list[str] | None = None,
) -> dict:
    mech_p, mic_p = _write_audits(tmp_path, mech_rows, mic_rows)
    out_md = tmp_path / "merge.md"
    argv = [
        "--drug", "ceftriaxone",
        "--mech-audit", str(mech_p),
        "--mic-audit", str(mic_p),
        "--output", str(out_md),
    ]
    if extra_args:
        argv.extend(extra_args)
    rc = merge_mod.main(argv)
    assert rc == 0, f"merge main() returned {rc}"
    out_json = out_md.with_suffix(".json")
    return json.loads(out_json.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Schema contract — consumed by pipeline._load_audit_verdict
# ---------------------------------------------------------------------------


def test_merge_emits_top_level_gate_verdict_consumable_by_pipeline_loader(tmp_path: Path):
    """Schema test 1: top-level gate_verdict + per-strain row read by _load_audit_verdict."""
    mech_rows = [{
        "strain_id": "S1",
        "status": "OK",
        "mechanisms_present": ["acquired_beta_lactamase"],
        "primary_mechanism_class": "acquired_beta_lactamase",
        "n_hits": 1,
        "mech_hits": {},
        "accession": "GCF_000000001.1",
        "mlst": "11",
    }]
    mic_rows = [{
        "strain_id": "S1",
        "tier": "HIGH_R",
        "cohort_binary_label": 1,
        "n_mic_rows": 3,
        "n_ast_rows": 3,
        "detail": {"median_mic": 32.0, "clsi_call": "R", "eucast_call": "R"},
        "accession": "GCF_000000001.1",
        "mlst": "11",
    }]
    payload = _run_merge(tmp_path, mech_rows, mic_rows)

    assert "gate_verdict" in payload
    assert "per_strain" in payload
    assert payload["per_strain"][0]["strain_id"] == "S1"

    merge_json = tmp_path / "merge.json"
    verdict = _load_audit_verdict(merge_json, strain_id="S1", fallback_to_cohort_gate=True)
    assert verdict is not None
    assert verdict["cohort_gate_verdict"] == payload["gate_verdict"]
    assert verdict["noise_class"] == "CLEAN_R_primary_mechanism"
    assert verdict["mic_tier"] == "HIGH_R"


def test_merge_per_strain_primary_mechanisms_is_plural_list(tmp_path: Path):
    """Schema test 2: primary_mechanisms is a list, never a scalar string."""
    mech_rows = [{
        "strain_id": "S1",
        "status": "OK",
        "mechanisms_present": ["acquired_beta_lactamase", "ampC_hyperproduction"],
        "primary_mechanism_class": "acquired_beta_lactamase",
        "n_hits": 2,
        "mech_hits": {},
        "accession": "",
        "mlst": "",
    }]
    mic_rows = [{
        "strain_id": "S1",
        "tier": "HIGH_R",
        "cohort_binary_label": 1,
        "detail": {},
    }]
    payload = _run_merge(tmp_path, mech_rows, mic_rows)

    row = payload["per_strain"][0]
    assert isinstance(row["primary_mechanisms"], list), \
        "primary_mechanisms MUST be a list (pipeline._load_audit_verdict expects iterable)"
    assert isinstance(row["co_resistance_modifiers"], list)
    assert set(row["primary_mechanisms"]) == {"acquired_beta_lactamase", "ampC_hyperproduction"}


def test_suspend_verdict_propagates_through_load_audit_verdict(tmp_path: Path):
    """Schema test 3: SUSPEND_CONDITION_4 surfaces suspend_gate_fired=True downstream."""
    # 10 strains, 2 CLEAN_R, 8 NOISY/OTHER -> signal_quality = 2/10 = 0.20 -> SUSPEND.
    mech_rows = []
    mic_rows = []
    for i in range(10):
        sid = f"S{i}"
        if i < 2:
            mechs = ["acquired_beta_lactamase"]
            tier = "HIGH_R"
            label = 1
        else:
            mechs = []
            tier = "BORDERLINE"
            label = 1
        mech_rows.append({
            "strain_id": sid,
            "status": "OK",
            "mechanisms_present": mechs,
            "primary_mechanism_class": mechs[0] if mechs else "MISSING",
            "n_hits": len(mechs),
            "mech_hits": {},
            "accession": "",
            "mlst": "",
        })
        mic_rows.append({
            "strain_id": sid,
            "tier": tier,
            "cohort_binary_label": label,
            "detail": {},
        })
    payload = _run_merge(tmp_path, mech_rows, mic_rows)
    assert payload["gate_verdict"] == "SUSPEND_CONDITION_4", \
        f"expected SUSPEND_CONDITION_4 at signal_quality=0.20, got {payload['gate_verdict']}"

    merge_json = tmp_path / "merge.json"
    verdict = _load_audit_verdict(merge_json, strain_id="S0", fallback_to_cohort_gate=True)
    assert verdict is not None
    assert verdict["suspend_gate_fired"] is True
    assert "verdict_explanation" in verdict
    assert "SUSPEND" in verdict["verdict_explanation"].upper()


# ---------------------------------------------------------------------------
# Drift-catch — the two locked parameters
# ---------------------------------------------------------------------------


def test_signal_quality_uses_len_merged_denominator(tmp_path: Path):
    """Drift-catch 1: signal_quality MUST be clean / len(merged), NOT clean / (clean+suspect+opacity).

    The earlier-bundle formula `clean / (clean+suspect+opacity)` excludes NOISY/OTHER
    buckets from the denominator, which biases the gate above the K/N null baseline.
    This test pins the cipro-matching denominator.
    """
    # Build 10 strains: 3 CLEAN_R, 0 SUSPECT/OPAQUE, 7 NOISY_R_borderline.
    mech_rows = []
    mic_rows = []
    for i in range(10):
        sid = f"S{i}"
        if i < 3:
            mechs = ["acquired_beta_lactamase"]
            tier = "HIGH_R"
        else:
            mechs = []
            tier = "BORDERLINE"
        mech_rows.append({
            "strain_id": sid,
            "status": "OK",
            "mechanisms_present": mechs,
            "primary_mechanism_class": mechs[0] if mechs else "MISSING",
            "n_hits": len(mechs),
            "mech_hits": {},
            "accession": "",
            "mlst": "",
        })
        mic_rows.append({
            "strain_id": sid,
            "tier": tier,
            "cohort_binary_label": 1,
            "detail": {},
        })
    payload = _run_merge(tmp_path, mech_rows, mic_rows)

    # With locked formula: clean / len(merged) = 3 / 10 = 0.30 -> SUSPEND (below 0.40).
    # With old buggy formula: clean / (clean+suspect+opacity) = 3 / 3 = 1.00 -> RUN_FULL_AND_CLEAN.
    # If this test fails with verdict=RUN_FULL_AND_CLEAN, the old denominator regressed in.
    assert abs(payload["signal_quality"] - 0.30) < 1e-9, \
        f"signal_quality should be 3/10 = 0.30, got {payload['signal_quality']}"
    assert payload["gate_verdict"] == "SUSPEND_CONDITION_4", \
        "denominator regression — re-check scripts/drug_mechanism_phenotype_merge.py line containing signal_quality"


def test_default_suspend_threshold_matches_cipro_calibration():
    """Drift-catch 2: default --suspend-threshold MUST be 0.40 (NOT 0.50)."""
    # Force a fresh import + reflect the argparse default by parsing minimal args.
    importlib.reload(merge_mod)

    # Rebuild the parser directly (don't run main) to inspect the default.
    parser = argparse.ArgumentParser()
    parser.add_argument("--drug", required=False, default="ceftriaxone")
    parser.add_argument("--mech-audit", type=Path, default=Path("dummy"))
    parser.add_argument("--mic-audit", type=Path, default=Path("dummy"))
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--suspend-threshold", type=float, default=0.40)
    args = parser.parse_args([])
    assert args.suspend_threshold == 0.40, "smoke check — local parser builds at 0.40"

    # Now inspect the actual script's parser.
    src_text = Path("scripts/drug_mechanism_phenotype_merge.py").read_text(encoding="utf-8")
    assert "default=0.40" in src_text, (
        "default --suspend-threshold drifted away from 0.40 (cipro 2026-05-17 calibration). "
        "See plans/Cef_Audit_Aware_Packet_Design.md edit #2 lock."
    )
    assert "default=0.50" not in src_text, (
        "default=0.50 regressed — this was Codex's 2026-05-26 bundle bug. "
        "Cef calibration is 0.40; raising it weakens the SUSPEND gate."
    )
