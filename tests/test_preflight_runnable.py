"""Tests for scripts/preflight_runnable.py pure-logic helpers.

Read-only preflight; pure-function helpers tested with synthetic capability dicts.
Real-host integration (torch / docker / git) is sanity-only — these tests don't
require those tools to be present.
"""
from __future__ import annotations

from pathlib import Path

import pytest


def test_check_python_always_available():
    from scripts.preflight_runnable import check_python
    result = check_python()
    assert result["available"] is True
    assert "." in result["version"]


def test_check_path_exists_missing(tmp_path: Path):
    from scripts.preflight_runnable import check_path_exists
    missing = tmp_path / "does-not-exist.txt"
    result = check_path_exists(missing)
    assert result["available"] is False
    assert "not found" in result["reason"].lower()


def test_check_path_exists_present(tmp_path: Path):
    from scripts.preflight_runnable import check_path_exists
    f = tmp_path / "file.txt"
    f.write_text("hello")
    result = check_path_exists(f)
    assert result["available"] is True
    assert result["type"] == "file"
    assert result["size_bytes"] == 5


def test_check_path_exists_directory(tmp_path: Path):
    from scripts.preflight_runnable import check_path_exists
    d = tmp_path / "subdir"
    d.mkdir()
    result = check_path_exists(d)
    assert result["available"] is True
    assert result["type"] == "dir"


def test_evaluate_ep_passes_when_all_caps_satisfied():
    """Post-push checklist requires git + pytest only; if both available -> RUNNABLE."""
    from scripts.preflight_runnable import EP_PROFILES, evaluate_ep
    profile = EP_PROFILES["post_push_checklist"]
    caps = {
        "git": {"available": True},
        "pytest": {"available": True},
        "file_checks": {},
    }
    result = evaluate_ep("post_push_checklist", profile, caps)
    assert result.runnable is True
    assert result.missing_capabilities == []
    assert result.missing_files == []


def test_evaluate_ep_fails_when_gpu_missing():
    """EP-1A genome-input requires GPU; missing GPU -> NOT runnable."""
    from scripts.preflight_runnable import EP_PROFILES, evaluate_ep
    profile = EP_PROFILES["EP_1A_genome_input"]
    caps = {
        "gpu": {"available": False, "reason": "torch.cuda.is_available() False"},
        "refseq_cache": {"available": True},
        "model_pickles": {
            "data/processed/models/ciprofloxacin_nucleotide_transformer.pkl": {
                "available": True,
                "is_leakage_safe": True,
            },
        },
        "file_checks": {
            "data/processed/stage2_n150_cipro_cohort.parquet": {"available": True},
        },
    }
    result = evaluate_ep("EP_1A_genome_input", profile, caps)
    assert result.runnable is False
    assert "gpu" in result.missing_capabilities


def test_evaluate_ep_fails_when_gpu_vram_too_small():
    """GTX 860M-class hardware: GPU available but < 6 GiB -> NOT runnable for NT."""
    from scripts.preflight_runnable import EP_PROFILES, evaluate_ep
    profile = EP_PROFILES["EP_1A_genome_input"]
    caps = {
        "gpu": {"available": True, "fits_nt_v2_100m": False, "vram_gib": 4.0},
        "refseq_cache": {"available": True},
        "model_pickles": {
            "data/processed/models/ciprofloxacin_nucleotide_transformer.pkl": {
                "available": True,
                "is_leakage_safe": True,
            },
        },
        "file_checks": {
            "data/processed/stage2_n150_cipro_cohort.parquet": {"available": True},
        },
    }
    result = evaluate_ep("EP_1A_genome_input", profile, caps)
    assert result.runnable is False
    assert "gpu" in result.missing_capabilities
    assert "VRAM" in result.details["gpu"]


def test_evaluate_ep_fails_when_no_leakage_safe_pickle():
    """A pickle that exists but has cv_strategy != leave_one_accession_out -> missing model_pickle_real."""
    from scripts.preflight_runnable import EP_PROFILES, evaluate_ep
    profile = EP_PROFILES["EP_1A_genome_input"]
    caps = {
        "gpu": {"available": True, "fits_nt_v2_100m": True, "vram_gib": 12.0, "bitsandbytes_compatible": True},
        "refseq_cache": {"available": True},
        "model_pickles": {
            "data/processed/models/ciprofloxacin_nucleotide_transformer.pkl": {
                "available": True,
                "is_leakage_safe": False,  # e.g., cv_strategy="loso" or "strain_id"
                "cv_strategy": "loso",
            },
        },
        "file_checks": {
            "data/processed/stage2_n150_cipro_cohort.parquet": {"available": True},
        },
    }
    result = evaluate_ep("EP_1A_genome_input", profile, caps)
    assert result.runnable is False
    assert "model_pickle_real" in result.missing_capabilities


def test_evaluate_ep_passes_when_leakage_safe_pickle_present():
    """A leakage-safe pickle satisfies the model_pickle_real capability."""
    from scripts.preflight_runnable import EP_PROFILES, evaluate_ep
    profile = EP_PROFILES["EP_1A_genome_input"]
    caps = {
        "gpu": {"available": True, "fits_nt_v2_100m": True, "vram_gib": 12.0, "bitsandbytes_compatible": True},
        "refseq_cache": {"available": True},
        "model_pickles": {
            "data/processed/models/some_old.pkl": {"available": True, "is_leakage_safe": False, "cv_strategy": "loso"},
            "data/processed/models/ciprofloxacin_nucleotide_transformer.pkl": {
                "available": True,
                "is_leakage_safe": True,
                "cv_strategy": "leave_one_accession_out",
            },
        },
        "file_checks": {
            "data/processed/stage2_n150_cipro_cohort.parquet": {"available": True},
        },
    }
    result = evaluate_ep("EP_1A_genome_input", profile, caps)
    assert result.runnable is True
    assert result.missing_capabilities == []


def test_evaluate_ep_reports_missing_files():
    """Required cohort file absent -> NOT runnable, missing_files populated."""
    from scripts.preflight_runnable import EP_PROFILES, evaluate_ep
    profile = EP_PROFILES["EP_1_5_architecture_poc"]
    caps = {
        "gpu": {"available": True, "fits_nt_v2_100m": True, "vram_gib": 12.0, "bitsandbytes_compatible": True},
        "docker": {"available": True},
        "amrfinder_db": {"available": True, "latest_present": True},
        "file_checks": {
            "data/processed/gate_b_mini_cipro_cohort.parquet": {"available": False},
            "data/processed/gate_b_mini_cef_cohort.parquet": {"available": True},
            "data/processed/gate_b_mini_tet_cohort.parquet": {"available": True},
        },
    }
    result = evaluate_ep("EP_1_5_architecture_poc", profile, caps)
    assert result.runnable is False
    assert "data\\processed\\gate_b_mini_cipro_cohort.parquet" in result.missing_files \
        or "data/processed/gate_b_mini_cipro_cohort.parquet" in result.missing_files


def test_host_fingerprint_no_gpu():
    from scripts.preflight_runnable import host_fingerprint
    caps = {"gpu": {"available": False, "reason": "no CUDA"}}
    assert "CPU-only" in host_fingerprint(caps)


def test_host_fingerprint_with_gpu():
    from scripts.preflight_runnable import host_fingerprint
    caps = {"gpu": {"available": True, "device_name": "GTX 860M", "vram_gib": 4.0, "compute_capability": "5.0"}}
    fp = host_fingerprint(caps)
    assert "GTX 860M" in fp
    assert "4.0 GiB" in fp


def test_ep_profiles_pinned():
    """Pin the 5 EPs so accidental edits show up in code review."""
    from scripts.preflight_runnable import EP_PROFILES
    assert set(EP_PROFILES.keys()) == {
        "EP_1A_genome_input",
        "EP_1B_external_benchmark",
        "EP_1_5_architecture_poc",
        "cef_slice",
        "post_push_checklist",
    }
