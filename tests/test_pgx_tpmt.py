"""Pins the TPMT PGx cell — the first COMPOUND-allele cell (*3A = *3B + *3C in cis).

Validated 85/85 vs the GeT-RM CDC consolidated consensus (real-data number is a committed artifact; VCF
gitignored). These tests pin the Ensembl-verified coords, the function-pair CPIC phenotype, and — the
load-bearing part — the compound resolution: two SNPs in cis -> *3A, each alone -> *3B / *3C.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import PGX_GENES  # noqa: E402
from dna_decode.pgx import tpmt_catalog as tp  # noqa: E402
from dna_decode.pgx.compound_caller import assemble_compound_diplotype  # noqa: E402
from dna_decode.pgx.runner import call_tpmt  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
# tag -> (pos, ref, alt, rsid), Ensembl-verified GRCh38 (chr6 minus strand; genomic ALT)
_C = {"c460T": (18138997, "C", "T", "rs1800460"), "c719C": (18130687, "T", "C", "rs1142345")}


def _vcf(tmp_path, rows, name="tpmt.vcf"):
    lines = [_HEADER.rstrip("\n")]
    for tag, gt in rows:
        pos, ref, alt, rsid = _C[tag]
        lines.append(f"chr6\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _call(vcf, sample=None):
    return assemble_compound_diplotype(vcf, tp.COMPONENTS, tp.COMPOUND_RULES,
                                       reference_allele=tp.REFERENCE_ALLELE,
                                       phenotype_fn=tp.diplotype_phenotype, gene=tp.GENE, sample=sample)


def test_tpmt_in_pgx_genes():
    assert "tpmt" in PGX_GENES


def test_tpmt_coords_grounded():
    by = {d.star: d for d in tp.COMPONENTS}
    assert (by["c460T"].pos, by["c460T"].ref, by["c460T"].alt, by["c460T"].rsid) == (18138997, "C", "T", "rs1800460")
    assert (by["c719C"].pos, by["c719C"].ref, by["c719C"].alt, by["c719C"].rsid) == (18130687, "T", "C", "rs1142345")
    assert all(d.chrom == "6" for d in tp.COMPONENTS)


@pytest.mark.parametrize("a1,a2,pheno", [
    ("*1", "*1", "Normal Metabolizer"),
    ("*1", "*3A", "Intermediate Metabolizer"),
    ("*1", "*3C", "Intermediate Metabolizer"),
    ("*3A", "*3A", "Poor Metabolizer"),
    ("*3B", "*3C", "Poor Metabolizer"),
])
def test_tpmt_function_pair_phenotype(a1, a2, pheno):
    assert tp.diplotype_phenotype(a1, a2) == pheno
    assert tp.diplotype_phenotype(a2, a1) == pheno


# --- the LOAD-BEARING compound tests ---
def test_tpmt_star3a_compound_cis(tmp_path):
    # both component SNPs on the SAME haplotype (phased 1|0, 1|0) -> *3A, NOT *3B+*3C
    r = _call(_vcf(tmp_path, [("c460T", "1|0"), ("c719C", "1|0")]))
    assert r.diplotype == "*1/*3A" and r.phenotype == "Intermediate Metabolizer"


def test_tpmt_star3a_hom(tmp_path):
    r = _call(_vcf(tmp_path, [("c460T", "1|1"), ("c719C", "1|1")]))
    assert r.diplotype == "*3A/*3A" and r.phenotype == "Poor Metabolizer"


def test_tpmt_star3b_and_star3c_trans(tmp_path):
    # the two SNPs on OPPOSITE haplotypes (phased 1|0, 0|1) -> *3B / *3C, NOT *3A
    r = _call(_vcf(tmp_path, [("c460T", "1|0"), ("c719C", "0|1")]))
    assert set(r.diplotype.split("/")) == {"*3B", "*3C"} and r.phenotype == "Poor Metabolizer"


def test_tpmt_star3c_single_het(tmp_path):
    r = _call(_vcf(tmp_path, [("c719C", "0|1")]))
    assert r.diplotype == "*1/*3C"


def test_tpmt_all_ref(tmp_path):
    r = _call(_vcf(tmp_path, [("c460T", "0|0"), ("c719C", "0|0")]))
    assert r.diplotype == "*1/*1" and r.phenotype == "Normal Metabolizer"


def test_tpmt_runner_record(tmp_path):
    rec = call_tpmt(_vcf(tmp_path, [("c460T", "1|0"), ("c719C", "1|0")]), sample_id="T1")
    assert rec["gene"] == "TPMT"
    assert rec["diplotype"] == "*1/*3A"
    assert rec["phenotype_abbrev"] == "IM"
    assert rec["caller"]["is_compound_allele_caller"] is True
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True


def test_tpmt_cli(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_vcf(tmp_path, [("c460T", "1|1"), ("c719C", "1|1")])), "--gene", "tpmt", "--json-only"])
    assert rc == 3 or rc == 0
    assert "TPMT" in capsys.readouterr().out


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
