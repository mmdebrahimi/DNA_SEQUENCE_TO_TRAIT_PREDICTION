"""Offline tests for the TB Portals -> candidate-TSV adapter (pure mapping + the gate-G1 method filter)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.build_tbportals_candidates import (  # noqa: E402
    classify_drug, classify_method, detect_columns, normalize_result, pivot_dst, _accession_cols,
)


def test_normalize_result_maps_measured_calls():
    assert normalize_result("Resistant") == "R"
    assert normalize_result("Susceptible") == "S"
    assert normalize_result("Sensitive") == "S"
    assert normalize_result("R") == "R" and normalize_result("s") == "S"
    # not a clean measured call -> blank
    for v in ("Intermediate", "Indeterminate", "Not tested", "Contaminated", "", None, "unknown"):
        assert normalize_result(v) == ""


def test_classify_method_molecular_wins():
    assert classify_method("MGIT 960") == "phenotypic"
    assert classify_method("Lowenstein-Jensen proportion") == "phenotypic"
    assert classify_method("Hain LPA") == "molecular"
    assert classify_method("WGS") == "molecular"
    assert classify_method("Xpert MTB/RIF") == "molecular"
    # a row that mentions both -> molecular (genotype-derived is circular)
    assert classify_method("WGS-confirmed MGIT") == "molecular"
    assert classify_method("") == "unknown"
    assert classify_method("some novel assay") == "unknown"


def test_classify_drug():
    assert classify_drug("Rifampicin") == "RIF"
    assert classify_drug("Rifampin") == "RIF"
    assert classify_drug("Isoniazid") == "INH"
    assert classify_drug("INH") == "INH"
    assert classify_drug("Ethambutol") == ""


def test_detect_columns_and_accession_aliases():
    header = ["condition_id", "drug_name", "dst_result", "dst_method", "run_accession", "biosample"]
    cols = detect_columns(header)
    assert cols["id"] == "condition_id"
    assert cols["drug"] == "drug_name"
    assert cols["result"] == "dst_result"
    assert cols["method"] == "dst_method"
    acc = _accession_cols(header)
    assert acc["run_accession"] == "run_accession"
    assert acc["biosample_accession"] == "biosample"


def test_pivot_dst_long_format_drops_molecular_and_pivots():
    cols = {"id": "iso", "drug": "drug", "result": "res", "method": "method"}
    acc = {"run_accession": "run", "biosample_accession": "bios"}
    rows = [
        # isolate A: phenotypic RIF-R + INH-S (kept, pivoted to one row)
        {"iso": "A", "drug": "Rifampicin", "res": "Resistant", "method": "MGIT", "run": "SRR1", "bios": "SAMN1"},
        {"iso": "A", "drug": "Isoniazid", "res": "Susceptible", "method": "MGIT", "run": "SRR1", "bios": "SAMN1"},
        # isolate B: ONLY a molecular RIF call -> dropped entirely (gate G1)
        {"iso": "B", "drug": "Rifampicin", "res": "Resistant", "method": "Hain LPA", "run": "SRR2", "bios": ""},
        # isolate C: phenotypic INH-R only
        {"iso": "C", "drug": "Isoniazid", "res": "Resistant", "method": "LJ proportion", "run": "SRR3", "bios": ""},
    ]
    cands, stats = pivot_dst(rows, cols, acc, keep_methods=("phenotypic", "unknown"))
    by_id = {c["strain_id"]: c for c in cands}
    assert set(by_id) == {"A", "C"}                       # B dropped (molecular-only)
    assert by_id["A"]["rif_label"] == "R" and by_id["A"]["inh_label"] == "S"
    assert by_id["A"]["run_accession"] == "SRR1" and by_id["A"]["biosample_accession"] == "SAMN1"
    assert by_id["C"]["inh_label"] == "R" and by_id["C"]["rif_label"] == ""
    assert stats["dropped_method_molecular"] == 1
    assert stats["n_isolates_usable"] == 2


def test_pivot_dst_strict_phenotypic_drops_unknown_method():
    cols = {"id": "iso", "drug": "drug", "result": "res", "method": "method"}
    rows = [{"iso": "X", "drug": "Rifampicin", "res": "Resistant", "method": "mystery assay"}]
    # default keeps unknown
    cands, _ = pivot_dst(rows, cols, {}, keep_methods=("phenotypic", "unknown"))
    assert cands and cands[0]["rif_label"] == "R"
    # strict drops it
    cands2, stats2 = pivot_dst(rows, cols, {}, keep_methods=("phenotypic",))
    assert cands2 == []
    assert stats2["dropped_method_unknown"] == 1


def test_pivot_dst_conflicting_repeats_keep_first_and_flag():
    cols = {"id": "iso", "drug": "drug", "result": "res", "method": "method"}
    rows = [
        {"iso": "A", "drug": "Rifampicin", "res": "Resistant", "method": "MGIT"},
        {"iso": "A", "drug": "Rifampicin", "res": "Susceptible", "method": "MGIT"},
    ]
    cands, stats = pivot_dst(rows, cols, {}, keep_methods=("phenotypic",))
    assert cands[0]["rif_label"] == "R"        # first kept
    assert stats["conflict_rif"] == 1


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
