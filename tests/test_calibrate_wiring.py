"""Tests for the calibrate_organism -> call_resistance wiring (opt-in per-organism path).

Verifies:
  - organism with a CALIBRATED registry entry -> uses the calibrated counter+threshold+intrinsic exclusion
  - organism with an EXPRESSION_FLOOR entry  -> ABSTAIN (refuses to predict)
  - organism=None (or unknown organism)      -> backward-compat DRUG_RULE path, unchanged
  - the committed registry loads and has the expected shape
"""
import csv
from pathlib import Path

import pytest

from dna_decode.eval.amr_rules import (
    call_resistance, calibrated_rule_for, load_calibrated_registry,
)

_COLS = ["Element symbol", "Element name", "Method", "Class", "Subclass", "% Identity to reference"]


def _write_main_tsv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in _COLS})


def _qrdr(sym):
    return {"Element symbol": sym, "Method": "POINTX", "Class": "QUINOLONE", "Subclass": "FLUOROQUINOLONE"}

def _gene(sym, cls="QUINOLONE", sub="FLUOROQUINOLONE"):
    return {"Element symbol": sym, "Method": "EXACTX", "Class": cls, "Subclass": sub}


# synthetic registry (independent of the committed one, so the test is hermetic)
_REG = {"rules": {
    "Campylobacter|ciprofloxacin": {"verdict": "CALIBRATED", "counter": "qrdr_point", "threshold": 1,
                                    "intrinsic_families_excluded": [], "loo_balanced_accuracy": 1.0, "n": 30},
    "Klebsiella|ciprofloxacin": {"verdict": "CALIBRATED", "counter": "qrdr_point", "threshold": 2,
                                 "intrinsic_families_excluded": ["oqxA", "oqxB"],
                                 "loo_balanced_accuracy": 1.0, "n": 30},
    "Salmonella|ciprofloxacin": {"verdict": "CALIBRATED", "counter": "broad", "threshold": 1,
                                 "intrinsic_families_excluded": [], "loo_balanced_accuracy": 1.0, "n": 30},
    "Acinetobacter|meropenem": {"verdict": "EXPRESSION_FLOOR", "counter": "broad", "threshold": 3,
                                "intrinsic_families_excluded": ["blaADC", "blaOXA-51-family"],
                                "loo_balanced_accuracy": 0.13, "n": 30},
}}


def test_calibrated_campylobacter_single_mutation_is_R(tmp_path):
    # single gyrA point mutation -> R under threshold 1 (Campylobacter calibration)
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_qrdr("gyrA_T86I")])
    r = call_resistance(mt, "ciprofloxacin", organism="Campylobacter", registry=_REG)
    assert r["prediction"] == "R"
    assert r["resistance_threshold"] == 1
    assert "calibrated_organism_v1" in r["rule"]


def test_default_rule_calls_single_mutation_S(tmp_path):
    # SAME genome, no organism -> DRUG_RULE default (cipro threshold 2) -> single mutation is S
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_qrdr("gyrA_T86I")])
    r = call_resistance(mt, "ciprofloxacin")
    assert r["prediction"] == "S"
    assert "calibrated_organism" not in r["rule"]   # backward-compat: default path


def test_calibrated_salmonella_qnr_is_R_via_broad_counter(tmp_path):
    # qnr gene, ZERO qrdr point -> default qrdr_point rule would miss it; Salmonella broad@1 catches it
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_gene("qnrB19")])
    cal = call_resistance(mt, "ciprofloxacin", organism="Salmonella", registry=_REG)
    assert cal["prediction"] == "R"
    assert cal["counter"] == "broad"
    default = call_resistance(mt, "ciprofloxacin")          # qrdr_point counter -> 0 -> S
    assert default["prediction"] == "S"


def test_calibrated_klebsiella_excludes_intrinsic_oqxab(tmp_path):
    # S strain carrying only intrinsic oqxAB -> calibrated rule excludes the family -> S
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_gene("oqxA"), _gene("oqxB")])
    r = call_resistance(mt, "ciprofloxacin", organism="Klebsiella", registry=_REG)
    assert r["prediction"] == "S"
    assert "oqxA" in r["intrinsic_families_excluded"]


def test_expression_floor_organism_abstains(tmp_path):
    # Acinetobacter meropenem -> EXPRESSION_FLOOR -> ABSTAIN regardless of determinants present
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_gene("blaOXA-23", "CARBAPENEM", "CARBAPENEM")])
    r = call_resistance(mt, "meropenem", organism="Acinetobacter", registry=_REG)
    assert r["prediction"] == "ABSTAIN"
    assert "ABSTAIN" in r["caveat"]


def test_unknown_organism_falls_back_to_default(tmp_path):
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_qrdr("gyrA_S83L"), _qrdr("parC_S80I")])
    r = call_resistance(mt, "ciprofloxacin", organism="Nonexistensia", registry=_REG)
    assert "calibrated_organism" not in r["rule"]          # no entry -> default path
    assert r["prediction"] == "R"                          # 2 qrdr points >= default threshold 2


def test_explicit_threshold_overrides_calibrated(tmp_path):
    # passing resistance_threshold disables the calibrated path (explicit override wins)
    mt = tmp_path / "main.tsv"; _write_main_tsv(mt, [_qrdr("gyrA_T86I")])
    r = call_resistance(mt, "ciprofloxacin", resistance_threshold=2, organism="Campylobacter", registry=_REG)
    assert "calibrated_organism" not in r["rule"]
    assert r["prediction"] == "S"


# ---- the committed registry ----
def test_committed_registry_loads_and_has_expected_rules():
    reg = load_calibrated_registry()
    rules = reg.get("rules", {})
    assert "Salmonella|ciprofloxacin" in rules
    assert rules["Salmonella|ciprofloxacin"]["counter"] == "broad"
    assert rules["Acinetobacter|meropenem"]["verdict"] == "EXPRESSION_FLOOR"
    assert "IN-SAMPLE" in reg.get("_provenance", "")

def test_calibrated_rule_for_case_insensitive():
    assert calibrated_rule_for("salmonella", "ciprofloxacin", registry=_REG)["counter"] == "broad"
    assert calibrated_rule_for("SALMONELLA", "ciprofloxacin", registry=_REG) is not None
    assert calibrated_rule_for("Campylobacter", "meropenem", registry=_REG) is None
