"""Pin the generic curated-cell AMR-Portal scorer mechanics on synthetic isolates (no parquet, no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.amr_portal_curated_score import SPEC_FLOOR, score_cell  # noqa: E402


def _rule(dets):
    present = any((d.get("amr_element_symbol") or "").startswith("tet(M)") for d in dets)
    return {"prediction": "R" if present else "S", "determinant_present": present}


def test_generic_score_cell():
    geno = {"mut": [{"amr_element_symbol": "tet(M)"}], "wt": [{"amr_element_symbol": "mtrR"}]}
    iso = ([("mut", False, "R")] * 12 + [("wt", False, "S")] * 12
           + [("wt", False, "R")]            # FN
           + [("mut", False, "S")]           # FP
           + [("mut", True, "R")])           # leaked -> excluded
    r = score_cell(iso, geno, _rule)
    b = r["binary"]
    assert b["tp"] == 12 and b["tn"] == 12 and b["fn"] == 1 and b["fp"] == 1
    assert r["n_R"] == 13 and r["n_S"] == 13 and r["powered"] is True
    assert r["strata"]["determinant_present"]["R"] == 12 and r["strata_reproduced"] is True
    assert b["spec"] >= SPEC_FLOOR and r["headline"] == "SCORED"


def test_low_spec_is_indeterminate():
    # rule fires on everything -> spec collapses -> INDETERMINATE (the over-call falsifier)
    geno = {}
    rule = lambda dets: {"prediction": "R", "determinant_present": True}
    r = score_cell([("x", False, "R")] * 12 + [("x", False, "S")] * 12, geno, rule)
    assert r["binary"]["spec"] == 0.0 and r["headline"] == "INDETERMINATE"
