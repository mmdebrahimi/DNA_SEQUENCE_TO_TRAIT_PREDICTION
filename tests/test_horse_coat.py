"""Pin the horse base coat-colour Mendelian rule (MC1R-E epistatic to ASIP-A).

Definitional cases from Rieder et al. 2001 / UC Davis VGL (the deployed rule):
e/e -> chestnut (any agouti); E_ A_ -> bay; E_ a/a -> black.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.horse_coat import INDETERMINATE, call_horse_base_colour  # noqa: E402


def test_chestnut_is_epistatic():
    # e/e -> chestnut regardless of agouti (including when agouti is unknown/absent)
    assert call_horse_base_colour("ee", "AA") == "chestnut"
    assert call_horse_base_colour("ee", "aa") == "chestnut"
    assert call_horse_base_colour("ee", "") == "chestnut"     # agouti irrelevant under e/e


def test_bay_needs_E_and_A():
    assert call_horse_base_colour("EE", "AA") == "bay"
    assert call_horse_base_colour("Ee", "Aa") == "bay"
    assert call_horse_base_colour("EE", "Aa") == "bay"
    assert call_horse_base_colour("Ee", "AA") == "bay"


def test_black_is_E_with_aa():
    assert call_horse_base_colour("EE", "aa") == "black"
    assert call_horse_base_colour("Ee", "aa") == "black"


def test_indeterminate_on_missing_or_bad():
    assert call_horse_base_colour("", "AA") == INDETERMINATE       # MC1R missing
    assert call_horse_base_colour("E", "AA") == INDETERMINATE      # single allele
    assert call_horse_base_colour("EE", "") == INDETERMINATE       # need agouti when E present
    assert call_horse_base_colour("EE", "A") == INDETERMINATE      # agouti single allele


def test_allele_separators_tolerated():
    assert call_horse_base_colour("E/e", "a/a") == "black"
    assert call_horse_base_colour("e e", "A A") == "chestnut"
