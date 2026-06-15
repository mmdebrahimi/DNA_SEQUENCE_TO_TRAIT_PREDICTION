"""Unit tests for the Oxford label-emission driver (pure helpers; no network)."""
from __future__ import annotations

from dna_decode.data.external_mic_labels import MicValue
import scripts.build_oxford_labels as bol


def test_rekey_merges_same_biosample():
    data = {
        "ERR1": {"ciprofloxacin": {"mics": [MicValue(8.0, "=", "8")], "calls": set()}},
        "ERR2": {"ciprofloxacin": {"mics": [MicValue(16.0, "=", "16")], "calls": {"R"}}},
        "aliasX": {"ciprofloxacin": {"mics": [MicValue(0.06, "=", "0.06")], "calls": set()}},
    }
    resolved = {"ERR1": "SAMEA1", "ERR2": "SAMEA1", "aliasX": "SAMEA2"}  # ERR1+ERR2 same BioSample
    bs_data, dropped = bol.rekey_to_biosample(data, resolved)
    assert set(bs_data) == {"SAMEA1", "SAMEA2"}
    assert len(bs_data["SAMEA1"]["ciprofloxacin"]["mics"]) == 2   # merged
    assert bs_data["SAMEA1"]["ciprofloxacin"]["calls"] == {"R"}
    assert dropped == []


def test_rekey_drops_unresolved():
    data = {"missing": {"ciprofloxacin": {"mics": [MicValue(8.0, "=", "8")], "calls": set()}}}
    bs_data, dropped = bol.rekey_to_biosample(data, {})
    assert bs_data == {}
    assert dropped == ["missing"]


def test_manifest_rows_tier_and_label():
    bs_data = {
        "SAMEA1": {"ciprofloxacin": {"mics": [MicValue(16.0, "=", "16")], "calls": set()}},   # HIGH_R
        "SAMEA2": {"ciprofloxacin": {"mics": [MicValue(0.75, "=", "0.75")], "calls": set()}},  # BORDERLINE
        "SAMEA3": {"ciprofloxacin": {"mics": [MicValue(8.0, ">", ">8")], "calls": set()}},     # CENSORED_HIGH_R
    }
    rows = bol.manifest_rows_for_drug(bs_data, "ciprofloxacin")
    by_bs = {r["biosample"]: r for r in rows}
    assert by_bs["SAMEA1"]["label"] == "R" and by_bs["SAMEA1"]["strict"] is True
    assert by_bs["SAMEA2"]["label"] == "EXCLUDED" and by_bs["SAMEA2"]["strict"] is False
    assert by_bs["SAMEA3"]["censor_meta"] is True and by_bs["SAMEA3"]["label"] == "R"


def test_manifest_rows_decisive_is_relaxed_only():
    # A DECISIVE_R isolate (MIC 5) -> relaxed True / strict False (the RELAXED_EXTRA branch).
    bs_data = {"SAMEA1": {"ciprofloxacin": {"mics": [MicValue(5.0, "=", "5")], "calls": set()}}}
    row = bol.manifest_rows_for_drug(bs_data, "ciprofloxacin")[0]
    assert row["tier"] == "DECISIVE_R"
    assert row["label"] == "R" and row["strict"] is False and row["relaxed"] is True
    assert row["censor_meta"] is False


def test_manifest_rows_skips_biosample_without_drug():
    # A BioSample carrying a DIFFERENT drug slot produces no row for the queried drug.
    bs_data = {"SAMEA1": {"gentamicin": {"mics": [MicValue(16.0, "=", "16")], "calls": set()}}}
    assert bol.manifest_rows_for_drug(bs_data, "ciprofloxacin") == []


def test_rekey_partial_resolution_drops_only_unresolved():
    # Mixed: one native key resolves, one doesn't -> resolved kept, other in dropped.
    data = {
        "ERR1": {"ciprofloxacin": {"mics": [MicValue(8.0, "=", "8")], "calls": set()}},
        "missing": {"ciprofloxacin": {"mics": [MicValue(16.0, "=", "16")], "calls": set()}},
    }
    bs_data, dropped = bol.rekey_to_biosample(data, {"ERR1": "SAMEA1"})
    assert set(bs_data) == {"SAMEA1"}
    assert dropped == ["missing"]


def test_drug_inputs_shape():
    bs_data = {"SAMEA1": {"ciprofloxacin": {"mics": [MicValue(8.0, "=", "8")], "calls": {"R"}}}}
    iso_mics, iso_calls = bol._drug_inputs(bs_data, "ciprofloxacin")
    assert iso_mics == {"SAMEA1": ["8"]}        # raw tokens for build_drug_labels
    assert iso_calls == {"SAMEA1": {"R"}}


def test_manifest_keys_are_biosamples():
    from dna_decode.data.external_cohort_genomes import is_biosample_key
    bs_data = {"SAMN12345": {"ciprofloxacin": {"mics": [MicValue(16.0, "=", "16")], "calls": set()}}}
    rows = bol.manifest_rows_for_drug(bs_data, "ciprofloxacin")
    assert all(is_biosample_key(r["biosample"]) for r in rows)
