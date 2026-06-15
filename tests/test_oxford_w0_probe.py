"""Offline unit tests for the W0 probe summarizers."""
from __future__ import annotations

import scripts.oxford_w0_probe as w0


def test_candidate_key_cardinality():
    records = [
        {"run_accession": "ERR1", "sample_accession": "SAMEA1", "sample_alias": "a", "secondary_sample_accession": "ERS1"},
        {"run_accession": "ERR2", "sample_accession": "SAMEA1", "sample_alias": "a", "secondary_sample_accession": "ERS1"},
    ]
    card = w0.candidate_key_cardinality(records)
    assert card["run_accession"] == 2     # two distinct runs
    assert card["sample_accession"] == 1  # one BioSample (two runs share it)


def test_mic_table_summary_censoring_and_dups():
    rows = [
        {"iso": "k1", "CIP": ">8"},
        {"iso": "k1", "CIP": "4"},       # duplicate (k1, cipro) pair
        {"iso": "k2", "CIP": "<=0.5"},
        {"iso": "k3", "XYZ": "1"},       # non-pilot column ignored
    ]
    s = w0.mic_table_summary(rows, "iso", {"CIP": "ciprofloxacin", "XYZ": "amikacin"})
    assert s["n_rows"] == 4
    assert s["n_unique_keys"] == 3
    assert "k1|ciprofloxacin" in s["duplicate_isolate_drug"]
    assert s["censoring_by_drug"]["ciprofloxacin"][">"] == 1
    assert s["censoring_by_drug"]["ciprofloxacin"]["<="] == 1
    assert s["censoring_by_drug"]["ciprofloxacin"]["="] == 1
    assert "amikacin" not in s["censoring_by_drug"]   # non-pilot skipped


def test_resolution_summary():
    records = [{"run_accession": "ERR1", "sample_accession": "SAMEA1", "sample_alias": "aliasX"}]
    res = w0.resolution_summary(["SAMEA1", "aliasX", "missing9"], records)
    assert res["n_mic_keys"] == 3
    assert res["n_resolved"] == 2          # SAMEA1 + aliasX resolve, missing9 doesn't
    assert res["resolution_rate"] == round(2 / 3, 4)
    assert "missing9" in res["unresolved_sample"]


def test_read_table_csv_and_tsv(tmp_path):
    csv_p = tmp_path / "t.csv"
    csv_p.write_text("iso,CIP\nk1,8\n")
    assert w0.read_table(csv_p) == [{"iso": "k1", "CIP": "8"}]
    tsv_p = tmp_path / "t.tsv"
    tsv_p.write_text("iso\tCIP\nk1\t8\n")
    assert w0.read_table(tsv_p) == [{"iso": "k1", "CIP": "8"}]
