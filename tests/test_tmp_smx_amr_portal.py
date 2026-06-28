"""Pins the TMP-SMX EBI-AMR-Portal experimental scorer (Tier-1 of the unscored-cell triage).

Tests the pure SIR-binning + strata-reproduction gate on synthetic isolates (no parquet / no network).
The strata gate = sul+dfr is the highest-R stratum AND sul-only R-rate < 0.5 (reproducing the
Sci234/Oxford `(sul AND dfr)` pattern). Real-data scoring needs the gitignored D: parquet (manual).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.tmp_smx_amr_portal_validate import _symbols, score_cell_sir  # noqa: E402


def _det(sym):
    return {"amr_element_symbol": sym, "gene_symbol": sym}


def test_symbols_extracts_gene_names():
    dets = [{"amr_element_symbol": "sul1", "gene_symbol": "x"}, {"amr_element_symbol": "", "gene_symbol": "dfrA1"}]
    assert _symbols(dets) == ["sul1", "dfrA1"]


def test_strata_reproduced_clean_and_rule():
    # sul+dfr -> all R; sul-only -> all S; neither -> all S. The AND rule separates cleanly.
    geno = {
        "g1": [_det("sul1"), _det("dfrA1")], "g2": [_det("sul2"), _det("dfrA17")],
        "g3": [_det("sul1")], "g4": [_det("sul2")],            # sul-only
        "g5": [],                                               # neither
    }
    isolates = [("g1", False, "R"), ("g2", False, "R"),        # sul+dfr R
                ("g3", False, "S"), ("g4", False, "S"),        # sul-only S
                ("g5", False, "S")]                            # neither S
    r = score_cell_sir(isolates, geno)
    assert r["strata_reproduced"] is True and r["headline"] == "SCORED"
    assert r["strata"]["sul+dfr"]["r_rate"] == 1.0
    assert r["strata"]["sul-only"]["r_rate"] == 0.0
    assert r["binary"]["tp"] == 2 and r["binary"]["tn"] == 3          # AND rule: 2 R-calls correct, 3 S correct


def test_strata_not_reproduced_flags_indeterminate():
    # sul-only stratum is mostly-R (R-rate 0.75 >= 0.5) -> the AND rule's premise fails -> INDETERMINATE.
    geno = {"a": [_det("sul1")], "b": [_det("sul1")], "c": [_det("sul2")], "d": [_det("sul2")],
            "e": [_det("sul1"), _det("dfrA1")]}
    isolates = [("a", False, "R"), ("b", False, "R"), ("c", False, "R"), ("d", False, "S"),  # sul-only 3R/1S
                ("e", False, "R")]                                                            # sul+dfr 1R
    r = score_cell_sir(isolates, geno)
    assert r["strata_reproduced"] is False and r["headline"] == "INDETERMINATE"


def test_leaked_isolates_excluded():
    geno = {"g1": [_det("sul1"), _det("dfrA1")]}
    isolates = [("g1", True, "R")]                              # leaked -> skipped
    r = score_cell_sir(isolates, geno)
    assert r["n_R"] == 0 and r["n_S"] == 0
