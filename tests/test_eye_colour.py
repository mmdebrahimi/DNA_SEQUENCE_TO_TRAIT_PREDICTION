"""Pins the eye-colour decoder (the first off-pathogen cell) + its OpenSNP ingestion.

THE load-bearing test is strand-agnosticism: rs12913832 reports as A/G (dbSNP) OR C/T (23andMe complement),
and the call MUST be identical (sourced 2026-06-28: blue allele = {G,C}, brown = {A,T}; a memory-inversion
was caught here). Plus the DTC-file parser (23andMe + AncestryDNA shapes) + the free-text label binner.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.eye_colour import call_eye_colour  # noqa: E402
from scripts.eye_colour_opensnp_validate import bin_eye_colour, genotype_from_dtc_file  # noqa: E402


def test_dbsnp_strand_GG_blue_AA_brown():
    assert call_eye_colour("GG")["prediction"] == "blue"
    assert call_eye_colour("AA")["prediction"] == "brown"
    assert call_eye_colour("AG")["prediction"] == "intermediate"
    assert call_eye_colour("GA")["prediction"] == "intermediate"


def test_23andme_complement_strand_CC_blue_TT_brown_SAME_CALL():
    # the load-bearing invariant: the complement strand gives the IDENTICAL phenotype call
    assert call_eye_colour("CC")["prediction"] == "blue"      # CC == GG -> blue
    assert call_eye_colour("TT")["prediction"] == "brown"     # TT == AA -> brown
    assert call_eye_colour("CT")["prediction"] == "intermediate"


def test_slash_and_missing_genotypes():
    assert call_eye_colour("A/G")["prediction"] == "intermediate"
    assert call_eye_colour("--")["prediction"] == "INDETERMINATE"
    assert call_eye_colour("")["prediction"] == "INDETERMINATE"
    assert call_eye_colour("G")["prediction"] == "INDETERMINATE"      # single allele


def test_allele_counts_reported():
    c = call_eye_colour("GG")
    assert c["n_blue_alleles"] == 2 and c["n_brown_alleles"] == 0


def test_parse_23andme_file(tmp_path):
    f = tmp_path / "u1.txt"
    f.write_text("# rsid\tchromosome\tposition\tgenotype\n"
                 "rs1\t1\t100\tAA\n"
                 "rs12913832\t15\t28365618\tGG\n", encoding="utf-8")
    assert genotype_from_dtc_file(f) == "GG"


def test_parse_ancestrydna_file(tmp_path):
    f = tmp_path / "u2.txt"
    f.write_text("#AncestryDNA\nrsid\tchromosome\tposition\tallele1\tallele2\n"
                 "rs12913832\t15\t28365618\tC\tT\n", encoding="utf-8")
    assert genotype_from_dtc_file(f) == "CT"


def test_parse_missing_rsid_returns_none(tmp_path):
    f = tmp_path / "u3.txt"
    f.write_text("rs1\t1\t100\tAA\n", encoding="utf-8")
    assert genotype_from_dtc_file(f) is None


def test_bin_eye_colour():
    assert bin_eye_colour("Blue") == "blue"
    assert bin_eye_colour("dark brown") == "brown"
    assert bin_eye_colour("green") == "other"
    assert bin_eye_colour("blue-green") == "other"          # green wins (a real third class)
    assert bin_eye_colour("hazel") == "other"
    assert bin_eye_colour("") is None
