"""UGT1A1 cell — CPIC activity-score + the *28 TA-repeat STRUCTURAL WALL (tag-SNP proxy via rs887829)."""
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import ugt1a1_catalog as ug  # noqa: E402
from dna_decode.pgx.runner import call_ugt1a1  # noqa: E402


def _vcf(*rows: str) -> Path:
    body = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP\n" + "".join(rows))
    p = Path(tempfile.mktemp(suffix=".vcf"))
    p.write_text(body, encoding="utf-8")
    return p


def test_activity_score_and_phenotype_bins():
    assert ug.activity_score("*1", "*1") == 2.0
    assert ug.activity_score("*1", "*80") == 1.5
    assert ug.activity_score("*80", "*80") == 1.0
    assert ug.diplotype_phenotype("*1", "*1") == "Normal Metabolizer"
    assert ug.diplotype_phenotype("*1", "*80") == "Intermediate Metabolizer"
    assert ug.diplotype_phenotype("*80", "*80") == "Poor Metabolizer"     # *28/*28-equivalent (Gilbert)
    assert ug.diplotype_phenotype("*6", "*80") == "Poor Metabolizer"


def test_defining_coords_pinned_grch38():
    """Ensembl-GRCh38-verified (2026-07-07): *80 rs887829 chr2:233759924 C>T; *6 rs4148323 chr2:233760498 G>A."""
    by = {d.star: d for d in ug.CORE_DEFINING}
    assert (by["*80"].rsid, by["*80"].chrom, by["*80"].pos, by["*80"].ref, by["*80"].alt) == ("rs887829", "2", 233759924, "C", "T")
    assert (by["*6"].rsid, by["*6"].chrom, by["*6"].pos, by["*6"].ref, by["*6"].alt) == ("rs4148323", "2", 233760498, "G", "A")


def test_call_ugt1a1_het_tag_is_intermediate():
    # rs887829 (*80, *28-tag) het -> *1/*80 -> Intermediate
    rec = call_ugt1a1(_vcf("2\t233759924\trs887829\tC\tT\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert rec["gene"] == "UGT1A1"
    assert rec["diplotype"] == "*1/*80"
    assert rec["phenotype"] == "Intermediate Metabolizer"
    assert rec["activity_score"] == 1.5


def test_call_ugt1a1_hom_tag_is_poor():
    rec = call_ugt1a1(_vcf("2\t233759924\trs887829\tC\tT\t.\tPASS\t.\tGT\t1/1\n"), sample_id="S")
    assert rec["diplotype"] == "*80/*80"
    assert rec["phenotype"] == "Poor Metabolizer"


def test_structural_wall_flagged_load_bearing():
    """The *28 TA-repeat wall + LD-tag posture MUST be surfaced — never a silent direct-repeat claim."""
    rec = call_ugt1a1(_vcf("2\t233759924\trs887829\tC\tT\t.\tPASS\t.\tGT\t0/1\n"), sample_id="S")
    assert rec["star28_ta_repeat_unassessed"] is True
    assert rec["caller"]["is_ld_tag_proxy_not_direct_repeat_call"] is True
    assert "TA-repeat" in rec["caveat"] or "TA-REPEAT" in rec["caveat"]
    assert "NOT a clinical tool" in rec["caveat"]


def test_registered_in_pgx_genes():
    from dna_decode.pgx import PGX_GENES
    assert "ugt1a1" in PGX_GENES
