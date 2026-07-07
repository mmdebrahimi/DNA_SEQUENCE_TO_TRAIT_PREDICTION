"""NUDT15 cell — CPIC activity-score catalog + caller (the 2nd thiopurine-toxicity pharmacogene)."""
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import nudt15_catalog as nu  # noqa: E402
from dna_decode.pgx.runner import call_nudt15  # noqa: E402


def _vcf(*rows: str) -> Path:
    body = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP\n" + "".join(rows))
    p = Path(tempfile.mktemp(suffix=".vcf"))
    p.write_text(body, encoding="utf-8")
    return p


def test_activity_score_and_phenotype_bins():
    assert nu.activity_score("*1", "*1") == 2.0
    assert nu.activity_score("*1", "*3") == 1.0
    assert nu.activity_score("*3", "*3") == 0.0
    assert nu.diplotype_phenotype("*1", "*1") == "Normal Metabolizer"
    assert nu.diplotype_phenotype("*1", "*3") == "Intermediate Metabolizer"
    assert nu.diplotype_phenotype("*3", "*3") == "Poor Metabolizer"
    assert nu.diplotype_phenotype("*3", "unknownX") == "Indeterminate"


def test_defining_coord_pinned_grch38():
    """Ensembl-GRCh38-verified (2026-07-07): NUDT15 *3 rs116855232 chr13:48045719 C>T."""
    d = nu.CORE_DEFINING[0]
    assert (d.star, d.rsid, d.chrom, d.pos, d.ref, d.alt) == ("*3", "rs116855232", "13", 48045719, "C", "T")
    assert nu.ACTIVITY_VALUE["*3"] == 0.0 and nu.ACTIVITY_VALUE["*1"] == 1.0


def test_call_nudt15_het_no_function_is_intermediate():
    rec = call_nudt15(_vcf("13\t48045719\trs116855232\tC\tT\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert rec["gene"] == "NUDT15"
    assert rec["diplotype"] == "*1/*3"
    assert rec["phenotype"] == "Intermediate Metabolizer"
    assert rec["activity_score"] == 1.0


def test_call_nudt15_hom_no_function_is_poor():
    rec = call_nudt15(_vcf("13\t48045719\trs116855232\tC\tT\t.\tPASS\t.\tGT\t1/1\n"), sample_id="S")
    assert rec["diplotype"] == "*3/*3"
    assert rec["phenotype"] == "Poor Metabolizer"
    assert rec["activity_score"] == 0.0


def test_call_nudt15_reference_is_normal():
    rec = call_nudt15(_vcf("13\t48045719\trs116855232\tC\tT\t.\tPASS\t.\tGT\t0/0\n"), sample_id="S")
    assert rec["diplotype"] == "*1/*1"
    assert rec["phenotype"] == "Normal Metabolizer"


def test_record_faithful_to_cpic_not_clinical():
    rec = call_nudt15(_vcf("13\t48045719\trs116855232\tC\tT\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True
    assert "NOT a clinical tool" in rec["caveat"]
    assert rec["trait"] == "pgx_metabolizer_phenotype"


def test_registered_in_pgx_genes():
    from dna_decode.pgx import PGX_GENES
    assert "nudt15" in PGX_GENES
