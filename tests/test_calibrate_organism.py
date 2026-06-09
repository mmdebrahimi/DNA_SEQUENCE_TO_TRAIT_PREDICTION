"""Tests for dna_decode.eval.calibrate_organism — auto-select per-organism AMR rule config.

Synthetic fixtures encode the four boundary regimes the wider-AMR validation thread found, so the test is
self-contained (no data/raw dependency, CI-portable):
  1. single-mutation organism (Campylobacter-like)  -> threshold 1
  2. double-mutation organism (E. coli/Klebsiella-like, + intrinsic efflux) -> qrdr_point @ threshold 2
  3. qnr-content organism (Salmonella-like)          -> broad counter @ threshold 1
  4. intrinsic-gene organism (Acinetobacter-like)    -> intrinsic family auto-excluded
"""
from dna_decode.eval.calibrate_organism import (
    CalibratedRule, balanced_accuracy, calibrate, family_token, intrinsic_families,
)


# ---------- helpers ----------
def _cohort(specs):
    """specs: list of (qrdr_symbols, broad_symbols, label). Returns (strains, labels)."""
    strains = [{"qrdr_point": list(q), "broad": list(b)} for q, b, _ in specs]
    labels = [l for _, _, l in specs]
    return strains, labels


# ---------- unit: family_token ----------
def test_family_token_oxa51_family():
    for allele in ("blaOXA-68", "blaOXA-69", "blaOXA-66", "blaOXA-100", "blaOXA-312"):
        assert family_token(allele) == "blaOXA-51-family"

def test_family_token_strips_allele_number():
    assert family_token("blaADC-76") == "blaADC"
    assert family_token("blaTEM-1") == "blaTEM"
    assert family_token("qnrB19") == "qnrB"

def test_family_token_strips_point_suffix():
    assert family_token("gyrA_S83L") == "gyrA"
    assert family_token("parC_S80I") == "parC"

def test_family_token_idempotent():
    assert family_token(family_token("blaOXA-68")) == "blaOXA-51-family"
    assert family_token(family_token("blaADC-76")) == "blaADC"


# ---------- unit: balanced_accuracy ----------
def test_balanced_accuracy_perfect_separation():
    cl = [(2, "R"), (2, "R"), (0, "S"), (0, "S")]
    assert balanced_accuracy(cl, 2) == 1.0
    assert balanced_accuracy(cl, 3) == 0.5   # all called S -> sens 0


# ---------- unit: intrinsic_families ----------
def test_intrinsic_families_flags_ubiquitous_both_classes():
    # blaOXA-51-family in every strain (R and S); strong acquired only in R
    psf = [({"blaOXA-51-family", "blaOXA-23"}, "R")] * 5 + [({"blaOXA-51-family"}, "S")] * 5
    out = intrinsic_families(psf)
    assert "blaOXA-51-family" in out
    assert "blaOXA-23" not in out            # discriminative -> not flagged

def test_intrinsic_families_empty_when_one_class_missing():
    psf = [({"x"}, "R")] * 3
    assert intrinsic_families(psf) == []


# ---------- regime 1: single-mutation (Campylobacter) -> threshold 1 ----------
def test_regime_single_mutation_picks_threshold_1():
    specs = [(["gyrA_T86I"], ["gyrA_T86I"], "R")] * 8 + [([], [], "S")] * 8
    strains, labels = _cohort(specs)
    rule = calibrate(strains, labels, "ciprofloxacin")
    assert rule.threshold == 1
    assert rule.full_balanced_accuracy == 1.0
    assert rule.loo_balanced_accuracy == 1.0


# ---------- regime 2: double-mutation + intrinsic efflux -> qrdr_point @ threshold 2 ----------
def test_regime_double_mutation_picks_qrdr_threshold_2():
    # R: gyrA+parC double (qrdr=2). S: a SINGLE sub-clinical gyrA (qrdr=1) -> threshold 1 over-calls (FP),
    # threshold 2 is required. BOTH classes carry intrinsic oqxAB efflux in the broad set.
    R = (["gyrA_S83L", "parC_S80I"], ["gyrA_S83L", "parC_S80I", "oqxA", "oqxB"], "R")
    S = (["gyrA_S83L"], ["gyrA_S83L", "oqxA", "oqxB"], "S")
    strains, labels = _cohort([R] * 8 + [S] * 8)
    rule = calibrate(strains, labels, "ciprofloxacin")
    assert rule.counter == "qrdr_point"
    assert rule.threshold == 2
    assert rule.full_balanced_accuracy == 1.0
    # qrdr_point@1 over-calls the single-mutant S strains; broad@1 over-calls via oqxAB -> neither is 1.0
    assert rule.per_config["qrdr_point@1"] < 1.0
    assert rule.per_config["broad@1"] < 1.0


# ---------- regime 3: qnr-content (Salmonella) -> broad counter @ threshold 1 ----------
def test_regime_qnr_content_picks_broad_threshold_1():
    # R: some single-gyrA, some qnr-only (qnr in broad, NOT in qrdr_point). S: nothing.
    single = (["gyrA_S83F"], ["gyrA_S83F"], "R")
    qnr_only = ([], ["qnrB19"], "R")          # zero QRDR point -> qrdr_point counter can't see it
    S = ([], [], "S")
    strains, labels = _cohort([single] * 4 + [qnr_only] * 4 + [S] * 8)
    rule = calibrate(strains, labels, "ciprofloxacin")
    assert rule.counter == "broad"            # must switch counter to catch qnr
    assert rule.threshold == 1
    assert rule.full_balanced_accuracy == 1.0
    # qrdr_point can never catch the qnr-only strains -> its best bal_acc < 1
    assert max(v for k, v in rule.per_config.items() if k.startswith("qrdr_point@")) < 1.0


# ---------- regime 4: intrinsic-gene (Acinetobacter) -> family auto-excluded ----------
def test_regime_intrinsic_gene_excluded_restores_separation():
    # Every strain carries a DIFFERENT OXA-51-family allele (intrinsic). R additionally carry OXA-23.
    import itertools
    alleles = itertools.cycle(["blaOXA-66", "blaOXA-68", "blaOXA-69", "blaOXA-100"])
    R = [([], [next(alleles), "blaOXA-23"], "R") for _ in range(8)]
    S = [([], [next(alleles)], "S") for _ in range(8)]
    strains, labels = _cohort(R + S)
    rule = calibrate(strains, labels, "meropenem")
    assert "blaOXA-51-family" in rule.intrinsic_families_excluded
    # after excluding the intrinsic family, broad@1 separates (only OXA-23 remains, R-only)
    assert rule.full_balanced_accuracy == 1.0


# ---------- regime 5: EXPRESSION floor -> abstain (no presence config works) ----------
def test_regime_expression_floor_abstains():
    # Resistance is expression-driven: R and S strains are determinant-indistinguishable (the same single
    # acquired gene appears in half of EACH class; the real driver is unobservable derepression). No
    # presence config can separate -> verdict EXPRESSION_FLOOR, abstain.
    R = [([], ["blaACT"], "R")] * 4 + [([], [], "R")] * 4      # half the R carry the gene, half don't
    S = [([], ["blaACT"], "S")] * 4 + [([], [], "S")] * 4      # ...same for S -> no signal
    strains, labels = _cohort(R + S)
    rule = calibrate(strains, labels, "ceftriaxone")
    assert rule.verdict == "EXPRESSION_FLOOR"
    assert rule.loo_balanced_accuracy < 0.70
    assert "ABSTAIN" in rule.note

def test_passing_regimes_have_calibrated_verdict():
    specs = [(["gyrA_T86I"], ["gyrA_T86I"], "R")] * 8 + [([], [], "S")] * 8
    strains, labels = _cohort(specs)
    assert calibrate(strains, labels, "ciprofloxacin").verdict == "CALIBRATED"


# ---------- predict() round-trips the calibrated config ----------
def test_predict_applies_calibrated_rule():
    specs = [(["gyrA_T86I"], ["gyrA_T86I"], "R")] * 6 + [([], [], "S")] * 6
    strains, labels = _cohort(specs)
    rule = calibrate(strains, labels, "ciprofloxacin")
    assert rule.predict({"qrdr_point": ["gyrA_T86I"], "broad": ["gyrA_T86I"]}) == "R"
    assert rule.predict({"qrdr_point": [], "broad": []}) == "S"


# ---------- guards ----------
def test_calibrate_rejects_length_mismatch():
    import pytest
    with pytest.raises(ValueError):
        calibrate([{"qrdr_point": [], "broad": []}], ["R", "S"], "ciprofloxacin")

def test_calibrate_rejects_empty_cohort():
    import pytest
    with pytest.raises(ValueError):
        calibrate([], [], "ciprofloxacin")
