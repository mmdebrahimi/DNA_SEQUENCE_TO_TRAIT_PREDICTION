"""DPYD cell — CPIC activity-score catalog + caller (the fluoropyrimidine-toxicity pharmacogene).

Pins the activity-score bins (Amstutz 2018), the four actionable-haplotype coords (Ensembl-GRCh38-verified),
and the end-to-end caller on synthetic VCFs. Mirrors the CYP2C9 activity-score test shape.
"""
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import dpyd_catalog as dp  # noqa: E402
from dna_decode.pgx.runner import call_dpyd  # noqa: E402


def _vcf(*rows: str) -> Path:
    body = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP\n" + "".join(rows))
    p = Path(tempfile.mktemp(suffix=".vcf"))
    p.write_text(body, encoding="utf-8")
    return p


# ---- catalog: CPIC activity-score bins ----

def test_activity_score_sums_allele_values():
    assert dp.activity_score("*1", "*1") == 2.0
    assert dp.activity_score("*1", "*2A") == 1.0
    assert dp.activity_score("*2A", "*2A") == 0.0
    assert dp.activity_score("*1", "c.2846A>T") == 1.5
    assert dp.activity_score("c.2846A>T", "HapB3") == 1.0
    assert dp.activity_score("*2A", "unknownX") is None


def test_phenotype_bins_match_cpic():
    assert dp.diplotype_phenotype("*1", "*1") == "Normal Metabolizer"        # AS 2.0
    assert dp.diplotype_phenotype("*1", "*2A") == "Intermediate Metabolizer"  # AS 1.0
    assert dp.diplotype_phenotype("*1", "c.2846A>T") == "Intermediate Metabolizer"  # AS 1.5
    assert dp.diplotype_phenotype("*2A", "*2A") == "Poor Metabolizer"        # AS 0.0
    assert dp.diplotype_phenotype("*2A", "c.2846A>T") == "Poor Metabolizer"  # AS 0.5
    assert dp.diplotype_phenotype("*13", "unknownX") == "Indeterminate"


def test_no_function_alleles_are_zero_decreased_are_half():
    assert dp.ACTIVITY_VALUE["*2A"] == 0.0 and dp.ACTIVITY_VALUE["*13"] == 0.0
    assert dp.ACTIVITY_VALUE["c.2846A>T"] == 0.5 and dp.ACTIVITY_VALUE["HapB3"] == 0.5
    assert dp.ACTIVITY_VALUE["*1"] == 1.0


def test_catalog_four_actionable_haplotypes_all_chr1():
    stars = {d.star for d in dp.CORE_DEFINING}
    assert stars == {"*2A", "*13", "c.2846A>T", "HapB3"}
    assert all(d.chrom == "1" for d in dp.CORE_DEFINING)          # DPYD is chr1
    assert all(d.pos > 97_000_000 and d.pos < 97_600_000 for d in dp.CORE_DEFINING)  # DPYD locus band
    # every defining allele carries an activity value
    assert all(d.star in dp.ACTIVITY_VALUE for d in dp.CORE_DEFINING)


def test_defining_coords_pinned_grch38():
    """Ensembl-GRCh38-verified (2026-07-07). A silent coord drift would mis-call every real VCF."""
    by = {d.star: d for d in dp.CORE_DEFINING}
    assert (by["*2A"].pos, by["*2A"].ref, by["*2A"].alt) == (97450058, "C", "T")
    assert (by["*13"].pos, by["*13"].ref, by["*13"].alt) == (97515787, "A", "C")
    assert (by["c.2846A>T"].pos, by["c.2846A>T"].ref, by["c.2846A>T"].alt) == (97082391, "T", "A")
    assert (by["HapB3"].pos, by["HapB3"].ref, by["HapB3"].alt) == (97579893, "G", "C")


# ---- caller: end-to-end on synthetic VCFs ----

def test_call_dpyd_het_no_function_is_intermediate():
    rec = call_dpyd(_vcf("1\t97450058\trs3918290\tC\tT\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert rec["gene"] == "DPYD"
    assert rec["diplotype"] == "*1/*2A"
    assert rec["phenotype"] == "Intermediate Metabolizer"
    assert rec["activity_score"] == 1.0


def test_call_dpyd_hom_no_function_is_poor():
    rec = call_dpyd(_vcf("1\t97450058\trs3918290\tC\tT\t.\tPASS\t.\tGT\t1/1\n"), sample_id="S")
    assert rec["diplotype"] == "*2A/*2A"
    assert rec["phenotype"] == "Poor Metabolizer"
    assert rec["activity_score"] == 0.0


def test_call_dpyd_het_decreased_is_intermediate_as_1_5():
    rec = call_dpyd(_vcf("1\t97082391\trs67376798\tT\tA\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert "c.2846A>T" in rec["diplotype"]
    assert rec["activity_score"] == 1.5
    assert rec["phenotype"] == "Intermediate Metabolizer"


def test_call_dpyd_reference_is_normal():
    rec = call_dpyd(_vcf("1\t97450058\trs3918290\tC\tT\t.\tPASS\t.\tGT\t0/0\n"), sample_id="S")
    assert rec["diplotype"] == "*1/*1"
    assert rec["phenotype"] == "Normal Metabolizer"
    assert rec["activity_score"] == 2.0


def test_call_dpyd_record_is_faithful_to_cpic_not_clinical():
    rec = call_dpyd(_vcf("1\t97450058\trs3918290\tC\tT\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True
    assert "NOT a clinical tool" in rec["caveat"]
    assert rec["trait"] == "pgx_metabolizer_phenotype"
