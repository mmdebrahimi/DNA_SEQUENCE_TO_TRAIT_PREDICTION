"""Offline pins for the HIV v0.2 absolute-cutoff calibration (scripts/hiv_absolute_cutoff_validate.py).

No network: exercises the sourced-cutoff dict + the confusion math + the CUTOFF_UNAVAILABLE wall + a
synthetic end-to-end. The real numbers live in the wiki artifacts; these pin the CONTRACT (esp. that
DOR + all INSTI are NEVER assigned a fabricated cutoff)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.hiv_absolute_cutoff_validate import (  # noqa: E402
    DRMCV_LOWER_CUTOFF, _confusion_at_cutoff, run_class,
)


def test_cutoffs_sourced_pi_and_nnrti_only():
    # PI: all 8
    for d in ("fosamprenavir", "atazanavir", "indinavir", "lopinavir",
              "nelfinavir", "saquinavir", "tipranavir", "darunavir"):
        assert d in DRMCV_LOWER_CUTOFF
    # the differentiated PI cutoffs (the genuine v0.2 content) match DRMcv.R exactly
    assert DRMCV_LOWER_CUTOFF["lopinavir"] == 9.0
    assert DRMCV_LOWER_CUTOFF["tipranavir"] == 2.0
    assert DRMCV_LOWER_CUTOFF["darunavir"] == 10.0
    # NNRTI: 4 (all == 3)
    for d in ("efavirenz", "nevirapine", "etravirine", "rilpivirine"):
        assert DRMCV_LOWER_CUTOFF[d] == 3.0


def test_no_fabricated_cutoff_for_dor_or_insti():
    # DOR + every INSTI must be ABSENT (postdate DRMcv.R) — never a guessed value
    for d in ("doravirine", "raltegravir", "elvitegravir", "dolutegravir", "bictegravir", "cabotegravir"):
        assert d not in DRMCV_LOWER_CUTOFF


def test_confusion_at_cutoff_perfect_separation():
    pairs = [(50.0, True)] * 10 + [(1.0, False)] * 10   # called-R high fold, called-S low fold
    m = _confusion_at_cutoff(pairs, cutoff=3.0)
    assert m["sens"] == 1.0 and m["spec"] == 1.0 and m["balacc"] == 1.0
    # under-powered
    assert _confusion_at_cutoff([(1.0, False)] * 5, 3.0)["note"].startswith("under-powered")


def test_confusion_respects_cutoff_value():
    # a fold of 5 is R at cutoff 3 but S at cutoff 9
    pairs = [(5.0, True)] * 20
    assert _confusion_at_cutoff(pairs, 3.0)["confusion"]["tp"] == 20   # label_R
    assert _confusion_at_cutoff(pairs, 9.0)["confusion"]["fp"] == 20   # label_S, called R -> FP


def test_run_class_insti_all_unavailable(tmp_path):
    # minimal INI .Full header (no fold needed — every drug walls before scoring)
    cols = ["SeqID", "PtID", "Subtype", "Method", "RefID", "Type", "IsolateName", "SeqType",
            "RAL", "EVG", "DTG", "BIC", "CAB"] + [f"P{p}" for p in range(1, 289)]
    p = tmp_path / "INI_DataSet.Full.txt"
    p.write_text("\t".join(cols) + "\n", encoding="utf-8")
    res = run_class("INSTI", p)
    assert res["n_drugs_calibrated"] == 0
    assert all(m["status"] == "CUTOFF_UNAVAILABLE" for m in res["per_drug"].values())


def test_run_class_pi_end_to_end(tmp_path):
    cols = ["SeqID", "PtID", "Subtype", "Method", "RefID", "Type", "IsolateName", "SeqType",
            "FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"] + [f"P{p}" for p in range(1, 100)]
    lines = ["\t".join(cols)]
    for i in range(20):
        d = {c: "-" for c in cols}
        d.update({"SeqID": str(i), "Subtype": "B", "Method": "PhenoSense", "Type": "Clinical"})
        if i < 10:
            d["P84"] = "V"                                   # I84V major -> position-based R
            for code in ("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"):
                d[code] = "50"                               # high fold -> label R at any cutoff
        else:
            for code in ("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"):
                d[code] = "1"
        lines.append("\t".join(d[c] for c in cols))
    p = tmp_path / "PI_DataSet.Full.txt"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    res = run_class("PI", p)
    assert res["n_drugs_calibrated"] == 8
    drv = res["per_drug"]["darunavir"]
    assert drv["status"] == "CALIBRATED" and drv["cutoff"] == 10.0
    assert drv["all"]["sens"] == 1.0 and drv["all"]["spec"] == 1.0    # perfect synthetic separation
