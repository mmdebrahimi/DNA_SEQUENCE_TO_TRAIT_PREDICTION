"""Offline pins for the HIV within-subtype de-confounding check (scripts/hiv_within_subtype.py).

No network / no real dataset: exercises the pure helpers (subtype grouping, mutant vs position observed
extraction, the frozen verdict) + a tiny synthetic .Full TSV end-to-end. Real-data numbers live in the
wiki artifacts; these pin the CONTRACT."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.hiv_amr import call_hiv_observed  # noqa: E402
from scripts.hiv_within_subtype import (  # noqa: E402
    CLASSES, _class_verdict, _group, _group_metrics, _observed, run_class,
)


def test_group_b_vs_nonb():
    assert _group("B") == "B"
    assert _group("C") == "non-B" and _group("CRF01_AE") == "non-B"
    assert _group("") == "unknown" and _group("Unknown") == "unknown"


def test_observed_mutant_uses_real_wt():
    # NNRTI (mutant-level, gene RT): K103N at P103; consensus '-' yields nothing.
    row = {"P103": "N", "P181": "-", "P190": ""}
    obs = _observed(row, "RT", CLASSES["NNRTI"]["positions"], "mutant")
    assert obs == {"RT": {"K103N"}}
    # and it drives an R call for an NNRTI drug through the frozen dispatch
    assert call_hiv_observed("efavirenz", obs).prediction == "R"


def test_observed_position_placeholder_wt_still_calls_R():
    # PI (position-based, gene PR): a residue at major position 84 -> R regardless of the WT letter used.
    row = {f"P{p}": "-" for p in CLASSES["PI"]["positions"]}
    row["P84"] = "V"
    obs = _observed(row, "PR", CLASSES["PI"]["positions"], "position")
    assert obs == {"PR": {"X84V"}}
    assert call_hiv_observed("darunavir", obs).prediction == "R"
    # empty gene -> S
    empty = _observed({f"P{p}": "-" for p in CLASSES["PI"]["positions"]}, "PR",
                      CLASSES["PI"]["positions"], "position")
    assert call_hiv_observed("darunavir", empty).prediction == "S"


def test_group_metrics_underpowered_and_auc_direction():
    # under-powered
    assert _group_metrics([(1.0, False)] * 5)["note"].startswith("under-powered")
    # called-R isolates carry the high folds -> AUC ~ 1.0
    fc = [(100.0, True)] * 10 + [(1.0, False)] * 10
    m = _group_metrics(fc)
    assert m["auc"] is not None and m["auc"] > 0.95


def test_class_verdict_holds():
    per_drug = {d: {"all": {"auc": 0.90}, "B": {"auc": 0.89}} for d in ("a", "b", "c")}
    v = _class_verdict(per_drug)
    assert v["verdict"] == "HOLDS_WITHIN_SUBTYPE" and v["median_within_b_auc"] == 0.89


def test_class_verdict_subtype_inflated():
    # pooled materially above within-B -> the class-mixed number rode subtype structure
    per_drug = {d: {"all": {"auc": 0.88}, "B": {"auc": 0.60}} for d in ("a", "b")}
    assert _class_verdict(per_drug)["verdict"] == "SUBTYPE_INFLATED"


def test_class_verdict_underpowered():
    per_drug = {"a": {"all": {"auc": 0.9}, "B": {"note": "under-powered (<15)"}}}
    assert _class_verdict(per_drug)["verdict"] == "WITHIN_B_UNDERPOWERED"


def test_end_to_end_synthetic_tsv(tmp_path):
    # Build a tiny PI .Full-shaped TSV: header + rows, one major-position mutant carrying high fold.
    cols = ["SeqID", "PtID", "Subtype", "Method", "RefID", "Type", "IsolateName", "SeqType",
            "FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"] + \
           [f"P{p}" for p in range(1, 100)]
    lines = ["\t".join(cols)]
    for i in range(20):
        d = {c: "-" for c in cols}
        d.update({"SeqID": str(i), "Subtype": "B", "Method": "PhenoSense", "Type": "Clinical"})
        if i < 10:                       # resistant: I84V + high fold
            d["P84"] = "V"
            for code in ("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"):
                d[code] = "50"
        else:                            # susceptible: WT + low fold
            for code in ("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"):
                d[code] = "1"
        lines.append("\t".join(d[c] for c in cols))
    p = tmp_path / "PI_DataSet.Full.txt"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    res = run_class("PI", p)
    assert res["drug_class"] == "PI" and res["n_clinical_phenosense"] == 20
    # perfect separation -> HOLDS
    assert res["class_verdict"]["verdict"] == "HOLDS_WITHIN_SUBTYPE"
    assert res["per_drug"]["darunavir"]["B"]["auc"] == 1.0
