"""Offline tests for the TB Portals -> candidate-TSV adapter (pure mapping + the gate-G1 method filter)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.build_tbportals_candidates import (  # noqa: E402
    classify_drug, classify_method, detect_columns, normalize_result, pivot_dst, _accession_cols,
    is_tbportals_wide, parse_tbportals_result, pivot_tbportals_wide,
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


# --- TB Portals WIDE format (verified vs the 2026 data dictionary) ---------------------------------------
def test_parse_tbportals_result_handles_aggregates():
    assert parse_tbportals_result("R") == "R"
    assert parse_tbportals_result("S") == "S"
    assert parse_tbportals_result("{R}") == "R"
    assert parse_tbportals_result("{S, R}") == ""        # discordant aggregate -> blank
    assert parse_tbportals_result("{S, S}") == "S"
    for v in ("Not Reported", "I", "Ind", "", None, "nan"):
        assert parse_tbportals_result(v) == ""


def test_is_tbportals_wide_detection():
    assert is_tbportals_wide(["condition_id", "bactec_rifampicin", "le_isoniazid"])
    assert not is_tbportals_wide(["iso", "drug", "result", "method"])


def test_pivot_tbportals_wide_keeps_phenotypic_drops_molecular():
    header = ["condition_id", "bactec_rifampicin", "le_rifampicin", "bactec_isoniazid",
              "genexpert_rifampicin", "hain_isoniazid", "ncbi_biosample", "ncbi_sra"]
    rows = [
        # A: bactec RIF-R, le RIF blank (agree-by-single), bactec INH-S; molecular cols must be IGNORED.
        {"condition_id": "A", "bactec_rifampicin": "R", "le_rifampicin": "Not Reported",
         "bactec_isoniazid": "S", "genexpert_rifampicin": "S", "hain_isoniazid": "R",
         "ncbi_biosample": "SAMN100", "ncbi_sra": "22001"},
        # B: bactec RIF-R and le RIF-S DISAGREE -> RIF blank; INH-R kept.
        {"condition_id": "B", "bactec_rifampicin": "R", "le_rifampicin": "S",
         "bactec_isoniazid": "R", "genexpert_rifampicin": "", "hain_isoniazid": "",
         "ncbi_biosample": "SAMN101", "ncbi_sra": "22002"},
        # C: no phenotypic call at all (only molecular) -> dropped.
        {"condition_id": "C", "bactec_rifampicin": "Not Reported", "le_rifampicin": "",
         "bactec_isoniazid": "Not Reported", "genexpert_rifampicin": "R", "hain_isoniazid": "R",
         "ncbi_biosample": "SAMN102", "ncbi_sra": "22003"},
    ]
    cands, stats = pivot_tbportals_wide(rows, header)
    by_id = {c["strain_id"]: c for c in cands}
    assert set(by_id) == {"A", "B"}                          # C dropped (molecular-only)
    assert by_id["A"]["rif_label"] == "R" and by_id["A"]["inh_label"] == "S"
    assert by_id["A"]["biosample_accession"] == "SAMN100"    # ncbi_biosample -> alias
    assert by_id["A"]["sample_accession"] == "22001"         # ncbi_sra -> alias
    assert by_id["B"]["rif_label"] == "" and by_id["B"]["inh_label"] == "R"   # bactec/le conflict -> blank
    assert stats["conflict_rif"] == 1
    assert stats["molecular_cols_ignored"] == 2              # genexpert_rifampicin + hain_isoniazid
    # molecular hain_isoniazid='R' on isolate A must NOT have leaked into A's INH label (it stayed 'S')
    assert by_id["A"]["inh_label"] == "S"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
