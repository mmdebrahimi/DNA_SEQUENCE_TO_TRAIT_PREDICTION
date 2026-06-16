"""Pin the non-frozen TMP-SMX overlay rule + co-trimoxazole tiering."""
import pytest

from dna_decode.data.experimental_drug_rules import (
    COTRIMOXAZOLE_BREAKPOINTS, cotrimoxazole_tier, is_dfr, is_sul, tmp_smx_call,
)


@pytest.mark.parametrize("genes, expected", [
    (["sul1", "dfrA17"], "R"),          # both families -> R
    (["sul2", "dfrA1", "blaCTX-M-15"], "R"),
    (["sul2"], "S"),                    # sul-only -> S
    (["dfrA1"], "S"),                   # dfr-only -> S
    (["blaCTX-M-15", "aac(3)-IIa"], "S"),  # neither -> S
    ([], "S"),
])
def test_and_truth_table(genes, expected):
    assert tmp_smx_call(genes)["prediction"] == expected


def test_matched_lists_reported():
    r = tmp_smx_call(["sul1", "sul2", "dfrA17", "blaTEM-1B"])
    assert r["matched_sul"] == ["sul1", "sul2"] and r["matched_dfr"] == ["dfrA17"]


@pytest.mark.parametrize("sym", ["sul1", "sul2", "sul3", "sul4", "sul2_1"])
def test_is_sul_accepts_acquired_families(sym):
    assert is_sul(sym)


@pytest.mark.parametrize("sym", ["sulR", "sulP", "sul", "sul5", "folP", "blaSUL"])
def test_is_sul_rejects_lookalikes(sym):
    assert not is_sul(sym)


@pytest.mark.parametrize("sym", ["dfrA17", "dfrA1", "dfrB1", "dfrA14b", "dfrA17_2"])
def test_is_dfr_accepts_acquired_families(sym):
    assert is_dfr(sym)


@pytest.mark.parametrize("sym", ["dfrA", "dfrB", "dfr", "folA", "dfrD"])
def test_is_dfr_rejects_lookalikes(sym):
    # bare family w/o allele number, target gene folA, and dfrD (not A/B) are excluded
    assert not is_dfr(sym)


def test_breakpoints_tmp_component():
    assert COTRIMOXAZOLE_BREAKPOINTS == {"clsi_r": 4.0, "clsi_s": 2.0, "eucast_r": 4.0, "eucast_s": 2.0}


@pytest.mark.parametrize("mic, tier", [
    (32, "HIGH_R"),     # >= 4*clsi_r (16)
    (16, "HIGH_R"),
    (0.25, "HIGH_S"),   # <= clsi_s/4 (0.5)
    (0.5, "HIGH_S"),
    (1, "BORDERLINE"),  # clean measured MIC the strict tier over-excludes
    (4, "BORDERLINE"),
])
def test_cotrimoxazole_tier(mic, tier):
    assert cotrimoxazole_tier([mic]) == tier


def test_tier_never_ambiguous_eucast_equals_clsi():
    # EUCAST == CLSI for co-trimoxazole here, so classify_tier must not return AMBIGUOUS
    assert all(cotrimoxazole_tier([m]) != "AMBIGUOUS" for m in (0.25, 0.5, 1, 2, 4, 8, 16, 32))
