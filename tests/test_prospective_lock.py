"""Tests for the prospective-lock layer (LA-3) — pure + offline (no network/Docker).

Pins: (1) the committed manifest is byte-true to the LIVE frozen decoder (the lock is valid right now);
(2) verify_lock detects tamper/drift; (3) eligibility is strictly-after + fail-closed; (4) the scorer's
pure parts (partition / powering / artifact stamp) behave; (5) the scorer HARD-FAILS on a broken lock and
runs the lock-only path clean.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.eval.prospective_lock import (  # noqa: E402
    FROZEN_SURFACE_FILES, SCORED_CELLS, compute_lock_manifest, is_prospective_eligible,
    surface_hashes, verify_lock,
)
from scripts.prospective_lock_validate import (  # noqa: E402
    build_artifact, partition_by_eligibility, powering_gate,
)
import scripts.prospective_lock_validate as plv  # noqa: E402

_REPO = Path(__file__).resolve().parent.parent
_MANIFEST = _REPO / "wiki" / "prospective_lock_manifest_2026-06-22.json"


def _manifest() -> dict:
    return json.loads(_MANIFEST.read_text(encoding="utf-8"))


# ---- the committed lock is valid against the live decoder ----

def test_committed_manifest_matches_live_decoder():
    """The committed manifest must hash-match the LIVE frozen surface — i.e. the decoder is still frozen."""
    assert _MANIFEST.exists(), "committed prospective-lock manifest missing"
    v = verify_lock(_manifest())
    assert v.ok, f"lock drifted={v.drifted} missing={v.missing} — decoder no longer matches the lock"


def test_manifest_pins_all_surface_files_and_cells():
    m = _manifest()
    assert set(m["surface_sha256"]) == set(FROZEN_SURFACE_FILES)
    assert len(m["scored_cells"]) == len(SCORED_CELLS)
    assert m["lock_date"] == "2026-06-22" and m["schema"] == "prospective-lock-manifest-v1"


def test_verify_lock_detects_tamper():
    m = _manifest()
    tampered = {**m, "surface_sha256": {**m["surface_sha256"],
                                        FROZEN_SURFACE_FILES[0]: "0" * 64}}
    v = verify_lock(tampered)
    assert not v.ok and FROZEN_SURFACE_FILES[0] in v.drifted


def test_verify_lock_detects_missing():
    v = verify_lock({"surface_sha256": {"dna_decode/data/does_not_exist.json": "x" * 64}})
    assert not v.ok and "dna_decode/data/does_not_exist.json" in v.missing


def test_surface_hashes_are_64hex():
    for h in surface_hashes().values():
        assert len(h) == 64 and all(c in "0123456789abcdef" for c in h)


# ---- temporal eligibility (the leakage-by-construction gate) ----

def test_eligibility_strictly_after_and_fail_closed():
    L = "2026-06-22"
    assert is_prospective_eligible("2026-06-23", L).eligible          # day after -> eligible
    assert is_prospective_eligible("2026-07", L).eligible             # next month (earliest=07-01) -> eligible
    assert not is_prospective_eligible("2026-06-22", L).eligible      # lock day itself -> NOT (not strictly after)
    assert not is_prospective_eligible("2026-06-01", L).eligible      # before -> NOT
    assert is_prospective_eligible("2026", L).reason == "undatable_fail_closed"   # year-only -> fail closed
    assert is_prospective_eligible(None, L).reason == "undatable_fail_closed"     # missing -> fail closed
    assert is_prospective_eligible("garbage", L).reason == "undatable_fail_closed"


def test_eligibility_requires_full_lock_date():
    try:
        is_prospective_eligible("2027-01-01", "2026")
    except ValueError:
        return
    raise AssertionError("expected ValueError on a non-full lock_date")


# ---- scorer pure parts ----

def test_partition_by_eligibility():
    cohort = [
        {"biosample": "A", "first_public_date": "2026-07-01"},   # post-lock
        {"biosample": "B", "first_public_date": "2026-05-01"},   # pre-lock
        {"biosample": "C", "first_public_date": None},           # undatable
    ]
    elig, excl = partition_by_eligibility(cohort, "2026-06-22")
    assert [r["biosample"] for r in elig] == ["A"]
    reasons = {r["biosample"]: r["exclusion_reason"] for r in excl}
    assert reasons == {"B": "pre_or_equal_lock", "C": "undatable_fail_closed"}


def test_powering_gate_states():
    assert powering_gate({"n_scored": 0})["status"] == "ACCRUING"
    assert powering_gate({"n_scored": 5, "tp": 2, "fn": 1, "tn": 1, "fp": 1})["status"] == "UNDERPOWERED"
    powered = powering_gate({"n_scored": 40, "tp": 10, "fn": 5, "tn": 12, "fp": 13})
    assert powered["status"] == "POWERED"


def test_build_artifact_stamps_lock_provenance():
    m = _manifest()
    art = build_artifact(m, "Escherichia_coli_Shigella", "ciprofloxacin", lock_ok=True,
                         eligible=[{"biosample": "A"}], excluded=[], conf={"n_scored": 0},
                         powering={"status": "ACCRUING"}, generated="2026-06-22")
    assert art["prospective_lock_verified"] is True
    assert art["lock_manifest"]["surface_sha256"] == m["surface_sha256"]   # provably the locked decoder
    assert art["lock_manifest"]["lock_date"] == "2026-06-22"


# ---- scorer CLI gates ----

def test_scorer_lock_only_path_passes(capsys):
    rc = plv.main(["--manifest", str(_MANIFEST)])   # no --cohort-tsv -> lock-only verification
    assert rc == 0
    assert "verify_lock OK" in capsys.readouterr().out


def test_scorer_hard_fails_on_broken_lock(tmp_path, capsys):
    m = _manifest()
    bad = {**m, "surface_sha256": {**m["surface_sha256"], FROZEN_SURFACE_FILES[0]: "0" * 64}}
    p = tmp_path / "broken_manifest.json"
    p.write_text(json.dumps(bad), encoding="utf-8")
    rc = plv.main(["--manifest", str(p)])
    assert rc == 2 and "LOCK BROKEN" in capsys.readouterr().err


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
