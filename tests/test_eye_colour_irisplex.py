"""Pin the IrisPlex 6-SNP v0.1 model (the deployed Walsh/HIrisPlex rule applied to independent data).

THE load-bearing test is the heterozygote RESCUE: rs12913832 AG (which v0 abstains on as 'intermediate')
resolves to a definite colour under the 6-SNP model — that is the entire reason v0.1 exists. Plus
strand-agnosticism (non-palindromic SNPs), the palindromic literal-strand handling, complete-case gating,
and the sourced-coefficient sanity (all-blue / all-brown homozygotes).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.eye_colour_irisplex import (  # noqa: E402
    effect_allele_count, predict_irisplex,
)

# a fully blue 6-SNP genotype: rs12913832 GG (0 brown-A); every other SNP homozygous NON-effect
BLUE_GT = {"rs12913832": "GG", "rs1800407": "CC", "rs12896399": "GG",
           "rs16891982": "GG", "rs1393350": "GG", "rs12203592": "CC"}
# fully brown: rs12913832 AA + brown-pushing effect alleles
BROWN_GT = {"rs12913832": "AA", "rs1800407": "TT", "rs12896399": "TT",
            "rs16891982": "CC", "rs1393350": "AA", "rs12203592": "TT"}


def test_homozygous_blue_predicts_blue():
    r = predict_irisplex(BLUE_GT)
    assert r["status"] == "PREDICTED"
    assert r["prediction"] == "blue"
    assert r["p_blue"] > 0.7


def test_homozygous_brown_predicts_brown():
    r = predict_irisplex(BROWN_GT)
    assert r["prediction"] == "brown"
    assert r["p_brown"] > 0.7


def test_heterozygote_rescue_AG_resolves_to_brown():
    """v0 abstains on rs12913832 AG (calls it 'intermediate'); v0.1 resolves it. The brown coefficient
    (5.41) dominates the intermediate one (3.16), so a lone AG -> brown — recovering exactly the abstained
    heterozygotes that are mostly brown in the OpenSNP data."""
    gt = {**BLUE_GT, "rs12913832": "AG"}
    r = predict_irisplex(gt)
    assert r["prediction"] == "brown"
    assert r["p_brown"] > r["p_blue"]


def test_strand_agnostic_nonpalindromic():
    # rs12913832 CC == GG (both 0 brown-effect; C is complement of blue-G)
    assert effect_allele_count("rs12913832", "CC") == 0
    assert effect_allele_count("rs12913832", "GG") == 0
    assert effect_allele_count("rs12913832", "AA") == 2
    assert effect_allele_count("rs12913832", "TT") == 2   # T = complement of brown-A
    # rs1800407 effect T: AA (complement) counts as 2; CC as 0
    assert effect_allele_count("rs1800407", "AA") == 2
    assert effect_allele_count("rs1800407", "CC") == 0


def test_palindromic_rs16891982_literal_only():
    # C/G palindromic: count LITERAL C only (no complement), per the forward-strand assumption
    assert effect_allele_count("rs16891982", "CC") == 2
    assert effect_allele_count("rs16891982", "CG") == 1
    assert effect_allele_count("rs16891982", "GG") == 0


def test_missing_rs12913832_is_indeterminate():
    gt = {k: v for k, v in BLUE_GT.items() if k != "rs12913832"}
    assert predict_irisplex(gt)["status"] == "MISSING_RS12913832"


def test_incomplete_snp_set_no_imputation():
    gt = {k: v for k, v in BLUE_GT.items() if k != "rs1393350"}
    r = predict_irisplex(gt)
    assert r["status"] == "INCOMPLETE_SNP_SET"
    assert r["missing"] == "rs1393350"


def test_uncallable_genotype_counts_none():
    assert effect_allele_count("rs12913832", "--") is None
    assert effect_allele_count("rs12913832", "A") is None   # single allele
