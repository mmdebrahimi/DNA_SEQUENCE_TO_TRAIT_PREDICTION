"""Pin the S. aureus-cipro AMR-Portal scorer mechanics on synthetic isolates (no parquet, no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.staph_cipro_amr_portal_validate import SPEC_FLOOR, score  # noqa: E402


def test_confusion_strata_and_endorsement():
    geno = {"mut": [{"amr_element_symbol": "parC_S80F"}], "wt": [{"amr_element_symbol": "mecA"}]}
    iso = ([("mut", False, "R")] * 12 + [("wt", False, "S")] * 12
           + [("wt", False, "R")]            # FN
           + [("mut", False, "S")]           # FP
           + [("mut", True, "R")])           # leaked -> excluded
    r = score(iso, geno)
    b = r["binary"]
    assert b["tp"] == 12 and b["tn"] == 12 and b["fn"] == 1 and b["fp"] == 1
    assert r["n_R"] == 13 and r["n_S"] == 13
    assert r["strata_reproduced"] is True and r["powered"] is True
    assert b["spec"] >= SPEC_FLOOR and r["headline"] == "SCORED"


def test_underpowered_when_few():
    geno = {"mut": [{"amr_element_symbol": "gyrA_S84L"}]}
    assert score([("mut", False, "R")] * 3 + [("mut", False, "S")] * 3, geno)["headline"] == "UNDERPOWERED"
