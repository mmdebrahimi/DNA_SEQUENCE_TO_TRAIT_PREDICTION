"""Pins the CYP2C9 (activity-score) + VKORC1 (warfarin sensitivity) PGx cells -- the warfarin pair."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import cyp2c9_catalog as c9  # noqa: E402
from dna_decode.pgx.caller import call_diplotype  # noqa: E402
from dna_decode.pgx.vkorc1 import call_vkorc1  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
# grounded GRCh38 coords
_C9 = {"*2": (94942290, "C", "T", "rs1799853"), "*3": (94981296, "A", "C", "rs1057910")}
_VKORC1 = (16, 31096368, "C", "T", "rs9923231")


def _c9_vcf(tmp_path, rows, name="c9.vcf"):
    lines = [_HEADER.rstrip("\n")]
    for star, gt in rows:
        pos, ref, alt, rsid = _C9[star]
        lines.append(f"chr10\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _call_c9(vcf, sample=None):
    return call_diplotype(vcf, sample=sample, defining=c9.CORE_DEFINING, sentinels=c9.SENTINELS,
                          reference_allele=c9.REFERENCE_ALLELE, phenotype_fn=c9.diplotype_phenotype,
                          gene=c9.GENE)


# --- catalog (CPIC activity-score) ---
def test_c9_coords_grounded():
    by = {d.star: d for d in c9.CORE_DEFINING}
    assert (by["*2"].pos, by["*2"].rsid) == (94942290, "rs1799853")
    assert (by["*3"].pos, by["*3"].rsid) == (94981296, "rs1057910")


@pytest.mark.parametrize("a1,a2,score,pheno", [
    ("*1", "*1", 2.0, "Normal Metabolizer"),
    ("*1", "*2", 1.5, "Intermediate Metabolizer"),
    ("*1", "*3", 1.0, "Intermediate Metabolizer"),
    ("*2", "*2", 1.0, "Intermediate Metabolizer"),
    ("*2", "*3", 0.5, "Poor Metabolizer"),
    ("*3", "*3", 0.0, "Poor Metabolizer"),
])
def test_c9_activity_score_phenotype(a1, a2, score, pheno):
    assert c9.activity_score(a1, a2) == score
    assert c9.diplotype_phenotype(a1, a2) == pheno
    assert c9.diplotype_phenotype(a2, a1) == pheno


def test_c9_unknown_allele_indeterminate():
    assert c9.diplotype_phenotype("*2", "*8") == "Indeterminate"


# --- caller (VCF) ---
def test_c9_star3_hom_poor(tmp_path):
    vcf = _c9_vcf(tmp_path, [("*3", "1/1")])
    r = _call_c9(vcf)
    assert r.diplotype == "*3/*3" and r.phenotype == "Poor Metabolizer"


def test_c9_star1_star2_het_intermediate(tmp_path):
    vcf = _c9_vcf(tmp_path, [("*2", "0/1")])
    r = _call_c9(vcf)
    assert r.diplotype == "*1/*2" and r.phenotype == "Intermediate Metabolizer"


def test_c9_all_ref_normal(tmp_path):
    vcf = _c9_vcf(tmp_path, [("*2", "0/0"), ("*3", "0/0")])
    r = _call_c9(vcf)
    assert r.diplotype == "*1/*1" and r.phenotype == "Normal Metabolizer"


def test_c9_runner_record(tmp_path):
    from dna_decode.pgx.runner import call_cyp2c9
    vcf = _c9_vcf(tmp_path, [("*2", "0/1"), ("*3", "0/1")])  # unphased *2 + *3
    rec = call_cyp2c9(vcf, sample_id="W1")
    assert rec["gene"] == "CYP2C9"
    assert rec["activity_score"] == 0.5            # *2/*3
    assert rec["phenotype"] == "Poor Metabolizer"
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True


# --- VKORC1 (single-SNP, minus-strand) ---
def _vk_vcf(tmp_path, gt, name="vk.vcf"):
    chrom, pos, ref, alt, rsid = _VKORC1
    p = tmp_path / name
    p.write_text(_HEADER + f"chr{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n", encoding="utf-8")
    return p


def test_vkorc1_homalt_high_sensitivity(tmp_path):
    r = call_vkorc1(_vk_vcf(tmp_path, "1|1"))
    assert r["alt_count"] == 2 and r["cdna_genotype"] == "A/A"
    assert r["warfarin_sensitivity"] == "High sensitivity" and r["dose_category"] == "low_dose"


def test_vkorc1_het_intermediate(tmp_path):
    r = call_vkorc1(_vk_vcf(tmp_path, "0|1"))
    assert r["cdna_genotype"] == "G/A" and r["warfarin_sensitivity"] == "Intermediate sensitivity"


def test_vkorc1_ref_normal(tmp_path):
    r = call_vkorc1(_vk_vcf(tmp_path, "0|0"))
    assert r["alt_count"] == 0 and r["cdna_genotype"] == "G/G"
    assert r["warfarin_sensitivity"] == "Normal sensitivity"


def test_vkorc1_absent_assumed_reference(tmp_path):
    p = tmp_path / "empty.vcf"
    p.write_text(_HEADER + "chr16\t99999999\trsX\tA\tG\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    r = call_vkorc1(p)
    assert r["status"] == "assumed_reference"
    assert "assumed_reference_at_uncalled_site" in r["flags"]


def test_cli_gene_routing(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_c9_vcf(tmp_path, [("*3", "1/1")])), "--gene", "cyp2c9", "--json-only"])
    assert rc == 3 or rc == 0   # *3/*3 PM is a valid call -> exit 0
    out = capsys.readouterr().out
    assert "CYP2C9" in out and "Poor Metabolizer" in out
    rc2 = main([str(_vk_vcf(tmp_path, "1|1")), "--gene", "vkorc1"])
    assert rc2 == 0
    assert "High sensitivity" in capsys.readouterr().out


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
