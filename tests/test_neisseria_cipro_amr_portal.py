"""Pin the NG-cipro AMR-Portal scorer mechanics on synthetic isolates (no parquet, no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.neisseria_cipro_amr_portal_validate import SPEC_FLOOR, score  # noqa: E402


def test_confusion_strata_and_endorsement():
    # 12 R-gyrA-mut (TP), 12 S-no-mut (TN), 1 R-no-mut (FN), 1 S-with-mut (FP); 1 leaked (excluded)
    geno = {"mut": [{"amr_element_symbol": "gyrA_S91F"}], "wt": [{"amr_element_symbol": "mtrR"}]}
    iso = ([("mut", False, "R")] * 12 + [("wt", False, "S")] * 12
           + [("wt", False, "R")]            # FN (R without gyrA mut)
           + [("mut", False, "S")]           # FP (S with gyrA mut)
           + [("mut", True, "R")])           # leaked -> excluded
    r = score(iso, geno)
    b = r["binary"]
    assert b["tp"] == 12 and b["tn"] == 12 and b["fn"] == 1 and b["fp"] == 1
    assert r["n_R"] == 13 and r["n_S"] == 13            # leaked R excluded
    assert r["strata"]["gyrA_qrdr_present"]["R"] == 12 and r["strata"]["gyrA_qrdr_absent"]["R"] == 1
    assert r["strata_reproduced"] is True               # present r_rate > absent r_rate
    assert r["powered"] is True and b["spec"] >= SPEC_FLOOR
    assert r["headline"] == "SCORED"


def test_underpowered_when_few():
    geno = {"mut": [{"amr_element_symbol": "gyrA_S91F"}]}
    r = score([("mut", False, "R")] * 3 + [("mut", False, "S")] * 3, geno)
    assert r["headline"] == "UNDERPOWERED"
