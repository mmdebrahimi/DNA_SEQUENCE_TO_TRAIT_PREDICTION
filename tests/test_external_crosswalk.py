"""Unit tests for the alias->BioSample crosswalk + conflict taxonomy."""
from __future__ import annotations

import json

from dna_decode.data import external_crosswalk as xw


def test_many_runs_one_biosample_ok():
    # two runs share one BioSample -> both run keys resolve, NO conflict
    records = [
        {"run_accession": "ERR1", "sample_accession": "SAMEA1", "sample_alias": "alphaA"},
        {"run_accession": "ERR2", "sample_accession": "SAMEA1", "sample_alias": "alphaA"},
    ]
    built = xw.build_crosswalk(records)
    assert built["conflicts"] == []
    assert built["crosswalk"]["ERR1"] == "SAMEA1"
    assert built["crosswalk"]["ERR2"] == "SAMEA1"
    assert built["crosswalk"]["alphaA"] == "SAMEA1"
    assert built["crosswalk"]["SAMEA1"] == "SAMEA1"   # BioSample maps to itself


def test_one_key_two_biosamples_is_conflict():
    # sample_alias "dup" appears on two different BioSamples -> HARD CONFLICT
    records = [
        {"run_accession": "ERR1", "sample_accession": "SAMEA1", "sample_alias": "dup"},
        {"run_accession": "ERR2", "sample_accession": "SAMEA2", "sample_alias": "dup"},
    ]
    built = xw.build_crosswalk(records)
    assert "dup" not in built["crosswalk"]
    conflict_keys = {c["mic_key"] for c in built["conflicts"]}
    assert "dup" in conflict_keys
    # provenance fields present
    c = next(c for c in built["conflicts"] if c["mic_key"] == "dup")
    assert set(c) >= {"mic_key", "candidate_field", "candidate_value", "resolved_biosample", "source_row_id"}


def test_cross_field_collision_is_conflict():
    # value "X" is a run_accession on one BioSample and a sample_alias on another
    records = [
        {"run_accession": "X", "sample_accession": "SAMEA1"},
        {"run_accession": "ERR9", "sample_accession": "SAMEA2", "sample_alias": "X"},
    ]
    built = xw.build_crosswalk(records)
    assert "X" not in built["crosswalk"]
    assert "X" in {c["mic_key"] for c in built["conflicts"]}


def test_resolve_keys_resolved_unresolved_conflicts():
    records = [
        {"run_accession": "ERR1", "sample_accession": "SAMEA1", "sample_alias": "good"},
        {"run_accession": "ERR2", "sample_accession": "SAMEA2", "sample_alias": "dup"},
        {"run_accession": "ERR3", "sample_accession": "SAMEA3", "sample_alias": "dup"},
    ]
    out = xw.resolve_keys(["good", "dup", "missing"], records)
    assert out["resolved"] == {"good": "SAMEA1"}
    assert out["unresolved"] == ["missing"]
    assert "dup" in {c["mic_key"] for c in out["conflicts"]}


def test_resolve_keys_conflicts_scoped_to_wanted():
    records = [
        {"run_accession": "ERR2", "sample_accession": "SAMEA2", "sample_alias": "dup"},
        {"run_accession": "ERR3", "sample_accession": "SAMEA3", "sample_alias": "dup"},
        {"run_accession": "ERR1", "sample_accession": "SAMEA1", "sample_alias": "good"},
    ]
    out = xw.resolve_keys(["good"], records)   # don't ask about "dup"
    assert out["resolved"] == {"good": "SAMEA1"}
    assert out["conflicts"] == []              # dup conflict not in scope


def test_records_missing_biosample_skipped():
    # A record with empty/missing sample_accession contributes NO crosswalk entries
    # (it has no resolution target); other records still resolve.
    records = [
        {"run_accession": "ERR1", "sample_accession": ""},        # no BioSample -> skipped
        {"run_accession": "ERR2"},                                # missing field -> skipped
        {"run_accession": "ERR3", "sample_accession": "SAMEA3"},
    ]
    built = xw.build_crosswalk(records)
    assert built["crosswalk"] == {"SAMEA3": "SAMEA3", "ERR3": "SAMEA3"}
    assert built["conflicts"] == []


def test_resolve_keys_strips_and_ignores_blank_native_keys():
    # Blank / whitespace native keys are dropped; surviving keys are stripped before lookup.
    records = [{"run_accession": "ERR1", "sample_accession": "SAMEA1"}]
    out = xw.resolve_keys([" ERR1 ", "", "   ", None], records)
    assert out["resolved"] == {"ERR1": "SAMEA1"}
    assert out["unresolved"] == []


def test_write_crosswalk_roundtrip(tmp_path):
    built = xw.build_crosswalk([{"run_accession": "ERR1", "sample_accession": "SAMEA1"}])
    p = tmp_path / "cw.json"
    xw.write_crosswalk(p, built)
    loaded = json.loads(p.read_text())
    assert loaded["crosswalk"]["ERR1"] == "SAMEA1"
