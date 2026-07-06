"""Pins the CYP2D6 STRUCTURAL surface (v0.2) — read-depth copy-number caller.

The real-data validation (CN-class concordance on 39 1000G CRAMs) is a committed artifact + a data-driven
test (tests/../data/pgx_1000g/cyp2d6_structural_ratios.tsv). These tests pin the pure ratio->copy-number
logic + the load-bearing honesty rail: the surface NEVER resolves hybrid allele identity (*13/*36/*68).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx.cyp2d6_structural import (  # noqa: E402
    HYBRID_D7_THRESHOLD,
    NORMAL_BASELINE,
    copy_number_from_ratio,
    hybrid_detection,
    hybrid_suspected,
    structural_call,
    truth_copy_number_class,
)

_RATIOS = REPO / "tests" / "data" / "pgx_getrm" / "cyp2d6_structural_ratios.tsv"
_HYBRID_RATIOS = REPO / "tests" / "data" / "pgx_getrm" / "cyp2d6_hybrid_ratios.tsv"


# --- pure copy-number logic (calibrated so a 2-copy sample ~= NORMAL_BASELINE) ---
@pytest.mark.parametrize("ratio,cn", [
    (0.0, 0), (0.30, 0),                 # homozygous deletion
    (0.58, 1), (0.63, 1), (0.66, 1),     # het deletion (observed *X/*5 band)
    (1.26, 2), (1.0, 2), (1.5, 2),       # normal 2-copy
    (1.73, 3), (1.89, 3),                # duplication (3 copies)
    (2.29, 4), (2.5, 4),                 # higher duplication
])
def test_copy_number_from_ratio(ratio, cn):
    assert copy_number_from_ratio(ratio) == cn


def test_copy_number_clamped_and_validated():
    assert copy_number_from_ratio(99.0) == 6         # clamp
    with pytest.raises(ValueError):
        copy_number_from_ratio(-1.0)


@pytest.mark.parametrize("ratio,status", [
    (0.30, "homozygous_deletion"), (0.62, "deletion"),
    (1.26, "normal_copy_number"), (1.9, "duplication"), (2.3, "duplication"),
])
def test_structural_call_status(ratio, status):
    assert structural_call(ratio).status == status


def test_hybrid_identity_never_resolved():
    """LOAD-BEARING honesty rail: this surface NEVER claims a hybrid allele (*13/*36/*68) identity."""
    for ratio in (0.6, 1.26, 2.3):
        c = structural_call(ratio)
        assert c.hybrid_identity_unresolved is True
        assert "*13" not in c.star_consistent and "PSV" in c.note


def test_truth_copy_number_class():
    assert truth_copy_number_class("*4/*5") == "deletion"
    assert truth_copy_number_class("*1/*5") == "deletion"
    assert truth_copy_number_class("*5/*5") == "homozygous_deletion"
    assert truth_copy_number_class("*15/*17") == "normal_copy_number"
    assert truth_copy_number_class("*1/*4") == "normal_copy_number"
    assert truth_copy_number_class("*2x2/*71") == "duplication"
    assert truth_copy_number_class("*68+*4/*5") is None       # hybrid -> not CN-scored (copy contrib unresolved)
    assert truth_copy_number_class("*1/*36x2+*10") is None    # *36 hybrid present -> excluded from CN scoring
    assert truth_copy_number_class("*1/(*68)+*4") is None     # pure hybrid -> excluded


# --- data-driven validation on the real measured ratios (committed artifact) ---
@pytest.mark.skipif(not _RATIOS.exists(), reason="committed structural ratios TSV not present")
def test_real_cn_class_concordance():
    import csv
    rows = list(csv.DictReader(_RATIOS.open(encoding="utf-8"), delimiter="\t"))
    scored = correct = 0
    for r in rows:
        tclass = truth_copy_number_class(r["truth"])
        if tclass is None:
            continue
        pred = structural_call(float(r["ratio"])).status
        scored += 1
        # deletion / homozygous_deletion are one copy-loss class for scoring
        dels = {"deletion", "homozygous_deletion"}
        correct += (pred == tclass) or (pred in dels and tclass in dels)
    assert scored >= 25
    assert correct / scored >= 0.85          # the copy-number surface separates DEL/NORMAL/DUP cleanly


# --- HYBRID DETECTION (v0.3) — CYP2D7 paralog depth signal ---
@pytest.mark.parametrize("d7_ratio,suspected", [
    (1.10, False), (1.16, False), (1.24, False),   # non-hybrid band
    (1.25, True), (1.30, True), (1.62, True),       # hybrid band (*68 family)
])
def test_hybrid_suspected(d7_ratio, suspected):
    assert hybrid_suspected(d7_ratio) is suspected


def test_hybrid_detection_never_resolves_identity():
    """LOAD-BEARING: detects PRESENCE only — never claims the exact *13/*36/*68 identity."""
    for d7 in (1.1, 1.4):
        assert hybrid_detection(d7)["hybrid_identity_resolved"] is False
    # a POSITIVE call names PRESENCE + defers identity to PSV analysis (never a specific allele call)
    pos = hybrid_detection(1.4)
    assert pos["hybrid_suspected"] is True and "PRESENCE" in pos["note"] and "PSV" in pos["note"]


@pytest.mark.skipif(not _HYBRID_RATIOS.exists(), reason="committed hybrid ratios TSV not present")
def test_hybrid_detection_high_specificity():
    """The CYP2D7-depth detector must SEPARATE hybrids from non-hybrids incl. pure dups (the confound),
    at HIGH specificity (a positive is trustworthy) — validated on the real 38-sample structural set."""
    import csv
    from scripts.cyp2d6_hybrid_validate import validate
    rows = list(csv.DictReader(_HYBRID_RATIOS.open(encoding="utf-8"), delimiter="\t"))
    rep = validate(rows)
    assert rep["n_hybrid_truth"] >= 10
    assert rep["specificity"] == 1.0            # never fires on a pure dup / normal (the confound test)
    assert rep["sensitivity"] >= 0.5            # partial (the *68 family detects cleanly)
    assert rep["auroc"] >= 0.75


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
