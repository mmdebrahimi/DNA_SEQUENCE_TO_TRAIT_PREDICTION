"""CYP4F2 (warfarin triad) + ABCG2 (statin pair) single-SNP function-readout cells."""
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx.cyp4f2 import call_cyp4f2  # noqa: E402
from dna_decode.pgx.abcg2 import call_abcg2  # noqa: E402
from dna_decode.pgx import cyp4f2 as cf, abcg2 as ab  # noqa: E402


def _vcf(chrom, pos, rsid, ref, alt, gt) -> Path:
    body = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP\n"
            f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    p = Path(tempfile.mktemp(suffix=".vcf"))
    p.write_text(body, encoding="utf-8")
    return p


def test_cyp4f2_coord_pinned_grch38():
    """Ensembl-GRCh38-verified (2026-07-07): CYP4F2 *3 rs2108622 chr19:15879621 C>T (T = 433M reduced)."""
    assert (cf.RSID, cf.CHROM, cf.POS, cf.REF, cf.ALT) == ("rs2108622", "19", 15879621, "C", "T")


def test_abcg2_coord_pinned_grch38():
    """Ensembl-GRCh38-verified (2026-07-07): ABCG2 rs2231142 chr4:88131171 G>T (T = 141K poor fn)."""
    assert (ab.RSID, ab.CHROM, ab.POS, ab.REF, ab.ALT) == ("rs2231142", "4", 88131171, "G", "T")


def test_cyp4f2_genotype_to_function_and_dose():
    ref = call_cyp4f2(_vcf("19", 15879621, "rs2108622", "C", "T", "0/0"))
    het = call_cyp4f2(_vcf("19", 15879621, "rs2108622", "C", "T", "0/1"))
    hom = call_cyp4f2(_vcf("19", 15879621, "rs2108622", "C", "T", "1/1"))
    assert ref["function"] == "Normal Function" and ref["warfarin_dose_direction"] == "no_adjustment"
    assert het["function"] == "Intermediate" and het["variant_genotype"] == "Val/Met"
    assert hom["function"] == "Reduced Function" and hom["warfarin_dose_direction"] == "higher_dose"
    assert hom["star_proxy"] == "*3/*3"


def test_abcg2_genotype_to_function_and_exposure():
    ref = call_abcg2(_vcf("4", 88131171, "rs2231142", "G", "T", "0/0"))
    het = call_abcg2(_vcf("4", 88131171, "rs2231142", "G", "T", "0/1"))
    hom = call_abcg2(_vcf("4", 88131171, "rs2231142", "G", "T", "1/1"))
    assert ref["function"] == "Normal Function" and ref["rosuvastatin_exposure"] == "typical_exposure"
    assert het["function"] == "Decreased Function" and het["variant_genotype"] == "Gln/Lys"
    assert hom["function"] == "Poor Function" and hom["rosuvastatin_exposure"] == "high_exposure"


def test_assumed_reference_flag_when_site_absent():
    rec = call_cyp4f2(_vcf("1", 999, "rsOther", "A", "G", "1/1"))  # no CYP4F2 record
    assert rec["status"] == "assumed_reference"
    assert "assumed_reference_at_uncalled_site" in rec["flags"]
    assert rec["function"] == "Normal Function"


def test_not_clinical_tool_caveat():
    assert "NOT a clinical tool" in call_cyp4f2(_vcf("19", 15879621, "rs2108622", "C", "T", "0/1"))["caveat"]
    assert "NOT a clinical tool" in call_abcg2(_vcf("4", 88131171, "rs2231142", "G", "T", "0/1"))["caveat"]


def test_registered_in_pgx_genes():
    from dna_decode.pgx import PGX_GENES
    assert "cyp4f2" in PGX_GENES and "abcg2" in PGX_GENES
