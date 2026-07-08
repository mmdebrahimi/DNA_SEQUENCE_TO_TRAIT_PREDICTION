"""Integrated multi-gene PGx interpretation (warfarin triad / statin pair / thiopurine pair)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import interpret as I  # noqa: E402


def test_warfarin_high_sensitivity_poor_metab_lowers_dose():
    r = I.interpret_warfarin("High sensitivity", "Poor Metabolizer", "Normal Function")
    assert r["status"] == "ok"
    assert r["dose_direction"] == "much_lower_dose_requirement"  # -2 (VKORC1) + -2 (CYP2C9) + 0


def test_warfarin_normal_all_is_standard():
    r = I.interpret_warfarin("Normal sensitivity", "Normal Metabolizer", "Normal Function")
    assert r["dose_direction"] == "standard_dose_requirement"


def test_warfarin_cyp4f2_star3_raises_dose():
    # Normal VKORC1 + Normal CYP2C9 + Reduced CYP4F2 (*3/*3) -> +1.0 -> higher
    r = I.interpret_warfarin("Normal sensitivity", "Normal Metabolizer", "Reduced Function")
    assert r["dose_direction"] == "higher_dose_requirement"


def test_warfarin_indeterminate_on_missing_gene():
    assert I.interpret_warfarin(None, "Normal Metabolizer", "Normal Function")["status"] == "indeterminate"


def test_statins_report_per_statin_risk():
    r = I.interpret_statins("Poor Function", "Normal Function")
    assert r["simvastatin_myopathy_risk"] == "high_risk"      # SLCO1B1-driven
    assert r["rosuvastatin_exposure_risk"] == "typical_risk"  # ABCG2 normal


def test_thiopurine_takes_worse_of_two_genes():
    # CPIC rule: the MORE-deficient gene governs. TPMT NM + NUDT15 PM -> high risk, governed by NUDT15.
    r = I.interpret_thiopurines("Normal Metabolizer", "Poor Metabolizer")
    assert r["toxicity_risk"] == "high_risk"
    assert r["governing_gene"] == "NUDT15"


def test_thiopurine_both_normal_is_normal_risk():
    r = I.interpret_thiopurines("Normal Metabolizer", "Normal Metabolizer")
    assert r["toxicity_risk"] == "normal_risk"
    assert r["governing_gene"] == "TPMT+NUDT15 (equal)"


def test_interpret_all_wires_realizer_shape():
    results = {
        "vkorc1": {"sensitivity": "High sensitivity"},
        "cyp2c9": {"phenotype": "Intermediate Metabolizer"},
        "cyp4f2": {"function": "Normal Function"},
        "slco1b1": {"function": "Decreased Function"},
        "abcg2": {"function": "Poor Function"},
        "tpmt": {"phenotype": "Intermediate Metabolizer"},
        "nudt15": {"phenotype": "Normal Metabolizer"},
    }
    out = I.interpret_all(results)
    assert out["warfarin"]["status"] == "ok"
    assert out["statins"]["rosuvastatin_exposure_risk"] == "high_risk"
    assert out["thiopurines"]["toxicity_risk"] == "increased_risk"  # IM in TPMT governs


def test_all_caveats_say_not_clinical():
    r = I.interpret_all({"vkorc1": {"sensitivity": "Normal sensitivity"},
                         "cyp2c9": {"phenotype": "Normal Metabolizer"},
                         "cyp4f2": {"function": "Normal Function"},
                         "slco1b1": {"function": "Normal Function"}, "abcg2": {"function": "Normal Function"},
                         "tpmt": {"phenotype": "Normal Metabolizer"}, "nudt15": {"phenotype": "Normal Metabolizer"}})
    for drug in ("warfarin", "statins", "thiopurines"):
        assert "NOT a clinical tool" in r[drug]["caveat"]
