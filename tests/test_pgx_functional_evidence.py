"""Pins the PGx independent functional-evidence layer (Unit A) -- the circularity-break cross-check."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx.functional_evidence import EVIDENCE, derive_missense_verdict, summary  # noqa: E402


@pytest.mark.parametrize("reduces,damaging,verdict", [
    (True, True, "AGREE"),     # CPIC reduced + predictor damaging -> agree
    (True, False, "DISAGREE"), # CPIC reduced but predictor benign -> disagree (the informative case)
    (False, False, "AGREE"),   # CPIC normal + predictor benign -> agree
    (False, True, "DISAGREE"),
    (True, None, "NO_SIGNAL"), # no usable prediction
])
def test_derive_missense_verdict(reduces, damaging, verdict):
    assert derive_missense_verdict(reduces, damaging) == verdict


def test_evidence_covers_the_warfarin_pair_plus_cyp2c19():
    keys = {(e.gene, e.allele) for e in EVIDENCE}
    assert ("CYP2C19", "*2") in keys and ("CYP2C19", "*3") in keys and ("CYP2C19", "*17") in keys
    assert ("CYP2C9", "*2") in keys and ("CYP2C9", "*3") in keys
    assert ("VKORC1", "-1639A") in keys


def test_cyp2c9_star3_is_the_disagree_finding():
    """The load-bearing honesty finding: predictors call *3 (I359L) benign but CPIC = no function."""
    e = next(e for e in EVIDENCE if e.gene == "CYP2C9" and e.allele == "*3")
    assert e.verdict == "DISAGREE"
    assert "BENIGN" in e.independent_signal.upper()


def test_cyp2c19_star2_is_flagged_synonymous_splice():
    e = next(e for e in EVIDENCE if e.gene == "CYP2C19" and e.allele == "*2")
    assert e.verdict == "FLAG" and "splice" in e.note.lower()


def test_cyp2c19_star3_stopgain_agrees():
    e = next(e for e in EVIDENCE if e.gene == "CYP2C19" and e.allele == "*3")
    assert e.verdict == "AGREE" and "stop" in e.variant_class


def test_summary_counts():
    s = summary()
    assert s["n"] == len(EVIDENCE)
    assert s["AGREE"] + s["DISAGREE"] + s["FLAG"] + s["NO_SIGNAL"] == s["n"]
    assert s["DISAGREE"] >= 1 and s["FLAG"] >= 1   # the layer surfaces real non-agreement


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
