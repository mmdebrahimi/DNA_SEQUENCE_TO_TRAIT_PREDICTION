"""Offline pins for the fair HIV resistance test's pure helpers (no torch / no model)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.esm_hiv_resistance_test import (  # noqa: E402
    R_CUTOFF, _single_mutant_variants, _wt_protein,
)


def test_single_mutant_extraction_and_label():
    wt = "MKVLA" * 20                                       # 100-residue synthetic protein
    # isolate A: single mutant at pos 3 (V->I), high fold -> R
    # isolate B: single mutant at pos 3 (V->I), low fold -> S  (same variant, diff isolate)
    # isolate C: TWO mutations -> excluded
    # isolate D: single mutant but no fold -> excluded
    rows = [
        {"Method": "PhenoSense", "EFV": "10", **{f"P{i}": "-" for i in range(1, 101)}},
        {"Method": "PhenoSense", "EFV": "1.2", **{f"P{i}": "-" for i in range(1, 101)}},
        {"Method": "PhenoSense", "EFV": "50", **{f"P{i}": "-" for i in range(1, 101)}},
        {"Method": "PhenoSense", "EFV": "", **{f"P{i}": "-" for i in range(1, 101)}},
    ]
    rows[0]["P3"] = "I"
    rows[1]["P3"] = "I"
    rows[2]["P3"] = "I"; rows[2]["P5"] = "W"               # two muts -> excluded
    rows[3]["P3"] = "I"
    got = _single_mutant_variants(rows, wt, 100, "EFV")
    assert len(got) == 2                                   # A + B only (C two-mut, D no-fold)
    assert all(v["pos"] == 3 and v["wt"] == "V" and v["mut"] == "I" for v in got)
    labels = sorted(v["label_R"] for v in got)
    assert labels == [False, True]                         # one R (fold 10), one S (fold 1.2)


def test_cutoff_boundary():
    wt = "MKVLA" * 20
    row = {"Method": "PhenoSense", "EFV": str(R_CUTOFF), **{f"P{i}": "-" for i in range(1, 101)}}
    row["P2"] = "R"
    v = _single_mutant_variants([row], wt, 100, "EFV")
    assert len(v) == 1 and v[0]["label_R"] is True         # fold == cutoff -> R (>=)


def test_non_phenosense_excluded():
    wt = "MKVLA" * 20
    row = {"Method": "Antivirogram", "EFV": "50", **{f"P{i}": "-" for i in range(1, 101)}}
    row["P3"] = "I"
    assert _single_mutant_variants([row], wt, 100, "EFV") == []


def test_hiv_references_translate():
    ref_dir = Path(__file__).resolve().parent.parent / "data/hiv_ref"
    if not (ref_dir / "HIV1_RT_HXB2_cds.fna").exists():
        import pytest
        pytest.skip("HIV references not present")
    rt = _wt_protein(ref_dir, "RT")
    pr = _wt_protein(ref_dir, "PR")
    assert len(rt) >= 440 and rt[0] == "P"                 # HXB2 RT starts PISPIET...
    assert len(pr) == 99 and pr[0] == "P"                  # HXB2 protease starts PQITL...; 99 aa
