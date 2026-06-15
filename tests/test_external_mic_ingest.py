"""Unit tests for the config-driven MIC-table ingester."""
from __future__ import annotations

import pytest

from dna_decode.data import external_mic_ingest as ing
from dna_decode.data.external_mic_labels import MicValue


def test_ingest_csv_drug_mapping_and_unmapped(tmp_path):
    p = tmp_path / "mic.csv"
    p.write_text("iso,CIP,AMK,note\nk1,>8,2,foo\nk2,0.06,4,bar\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin", "AMK": "amikacin"})
    # AMK is non-pilot -> skipped; CIP mapped
    assert "ciprofloxacin" in out["data"]["k1"]
    assert "amikacin" not in out["data"]["k1"]
    # "note" + "AMK" are unmapped-to-a-pilot... AMK IS in drug_cols (mapped col) so only "note" unmapped
    assert out["unmapped_columns"] == ["note"]


def test_ingest_preserves_micvalue_operator(tmp_path):
    p = tmp_path / "mic.csv"
    p.write_text("iso,CIP\nk1,>8\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"})
    mics = out["data"]["k1"]["ciprofloxacin"]["mics"]
    assert mics == [MicValue(8.0, ">", ">8")]


def test_ingest_tsv(tmp_path):
    p = tmp_path / "mic.tsv"
    p.write_text("iso\tCIP\nk1\t4\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"})
    assert out["data"]["k1"]["ciprofloxacin"]["mics"][0].value == 4.0


def test_ingest_collects_calls(tmp_path):
    p = tmp_path / "mic.csv"
    p.write_text("iso,CIP,CIP_SIR\nk1,8,R\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"},
                               call_cols={"CIP_SIR": "ciprofloxacin"})
    assert out["data"]["k1"]["ciprofloxacin"]["calls"] == {"R"}


def test_ingest_missing_key_col_raises(tmp_path):
    p = tmp_path / "mic.csv"
    p.write_text("wrongkey,CIP\nk1,8\n")
    with pytest.raises(ValueError):
        ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"})


def test_ingest_blank_key_rows_skipped(tmp_path):
    p = tmp_path / "mic.csv"
    p.write_text("iso,CIP\n,8\nk2,4\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"})
    assert list(out["data"].keys()) == ["k2"]


# --- gap coverage ---------------------------------------------------------- #
def test_ingest_empty_table_no_rows(tmp_path):
    # Header-only file: no rows -> empty data + no unmapped, and NO missing-key raise
    # (the key_col presence check is guarded on `rows` being non-empty).
    p = tmp_path / "empty.csv"
    p.write_text("iso,CIP\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"})
    assert out == {"data": {}, "unmapped_columns": []}


def test_ingest_call_col_non_pilot_skipped(tmp_path):
    # A call column mapped to a NON-pilot drug is skipped (canonical_drug -> None),
    # leaving the isolate present but with no drug slot from that call column.
    p = tmp_path / "c.csv"
    p.write_text("iso,AMK_SIR\nk1,R\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={},
                               call_cols={"AMK_SIR": "amikacin"})
    assert out["data"] == {"k1": {}}


def test_ingest_accumulates_multiple_rows_per_isolate(tmp_path):
    # Two AST rows for the same isolate/drug accumulate both MicValues (no overwrite).
    p = tmp_path / "m.csv"
    p.write_text("iso,CIP\nk1,8\nk1,16\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"CIP": "ciprofloxacin"})
    mics = out["data"]["k1"]["ciprofloxacin"]["mics"]
    assert [mv.value for mv in mics] == [8.0, 16.0]


def test_ingest_drug_col_non_pilot_creates_empty_slot(tmp_path):
    # A non-pilot DRUG column is skipped; an isolate with only that column is present
    # but carries no drug slot (loud: still in data, just empty).
    p = tmp_path / "d.csv"
    p.write_text("iso,AMK\nk1,2\n")
    out = ing.ingest_mic_table(p, key_col="iso", drug_cols={"AMK": "amikacin"})
    assert out["data"] == {"k1": {}}
