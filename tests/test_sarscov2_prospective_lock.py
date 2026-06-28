"""Pins the SARS-CoV-2 Mpro prospective-lock: surface-hash pin + tamper-evidence + date-eligibility.

Mirrors the AMR prospective-lock test intent for the SARS cell. Reuses the generic
`is_prospective_eligible` (date-based) primitive; touches no AMR-surface file.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.prospective_lock import is_prospective_eligible  # noqa: E402
from scripts.sarscov2_prospective_lock import (  # noqa: E402
    SARS_MPRO_SURFACE_FILES, compute_manifest, surface_hashes, verify_lock,
)

LOCK = "2026-06-27"


def test_surface_files_exist_and_hash():
    h = surface_hashes()
    assert set(h) == set(SARS_MPRO_SURFACE_FILES)
    assert all(len(v) == 64 for v in h.values())          # full sha256


def test_manifest_schema_and_live_surface_intact():
    man = compute_manifest(LOCK)
    assert man["schema"] == "sarscov2-mpro-prospective-lock-v1"
    assert man["lock_date"] == LOCK
    ok, drifted = verify_lock(man)
    assert ok and drifted == []                           # the just-built manifest matches the live surface


def test_verify_detects_drift():
    man = compute_manifest(LOCK)
    rel = SARS_MPRO_SURFACE_FILES[0]
    man["surface_sha256"][rel] = "0" * 64                 # tamper one hash
    ok, drifted = verify_lock(man)
    assert not ok and rel in drifted                      # tamper-evident


def test_eligibility_is_strictly_after_lock():
    assert is_prospective_eligible("2026-07-01", LOCK).eligible        # after lock -> independent
    assert not is_prospective_eligible("2026-06-13", LOCK).eligible    # before lock -> leakage risk
    assert not is_prospective_eligible(LOCK, LOCK).eligible            # ON lock -> not strictly after
    assert not is_prospective_eligible(None, LOCK).eligible            # unknown date -> fail-closed
