"""Pins the CYP3A5 PGx cell — tacrolimus expressor/non-expressor phenotype (function-pair, chr7).

First PGx gene outside the chr10 CYP2C cluster. Validated 8/8 vs the REAL GeT-RM CDC multi-lab consensus
(UNDERPOWERED n=8) — see wiki/pgx_getrm_concordance_cyp3a5_2026-07-05. These tests pin the Ensembl-verified
GRCh38 coords (incl. the *7 INSERTION), the function-pair CPIC phenotype, and the calling logic on
synthetic VCFs. The real-data 8/8 number is a committed artifact (VCF is gitignored).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import PGX_GENES  # noqa: E402
from dna_decode.pgx import cyp3a5_catalog as c3  # noqa: E402
from dna_decode.pgx.caller import call_diplotype  # noqa: E402
from dna_decode.pgx.runner import call_cyp3a5  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
# Ensembl-verified GRCh38 forward-strand coords (CYP3A5 chr7 minus strand). *7 is a VCF left-aligned insertion.
_C3 = {
    "*3": (99672916, "T", "C", "rs776746"),
    "*6": (99665212, "C", "T", "rs10264272"),
    "*7": (99652770, "T", "TA", "rs41303343"),   # insertion
}


def _c3_vcf(tmp_path, rows, name="c3.vcf"):
    lines = [_HEADER.rstrip("\n")]
    for star, gt in rows:
        pos, ref, alt, rsid = _C3[star]
        lines.append(f"chr7\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _call_c3(vcf, sample=None):
    return call_diplotype(vcf, sample=sample, defining=c3.CORE_DEFINING, sentinels=c3.SENTINELS,
                          reference_allele=c3.REFERENCE_ALLELE, phenotype_fn=c3.diplotype_phenotype,
                          gene=c3.GENE)


# --- catalog / provenance ---
def test_c3_in_pgx_genes():
    assert "cyp3a5" in PGX_GENES


def test_c3_coords_grounded_ensembl():
    by = {d.star: d for d in c3.CORE_DEFINING}
    assert (by["*3"].pos, by["*3"].ref, by["*3"].alt, by["*3"].rsid) == (99672916, "T", "C", "rs776746")
    assert (by["*6"].pos, by["*6"].ref, by["*6"].alt, by["*6"].rsid) == (99665212, "C", "T", "rs10264272")
    assert (by["*7"].pos, by["*7"].ref, by["*7"].alt, by["*7"].rsid) == (99652770, "T", "TA", "rs41303343")
    assert all(d.chrom == "7" for d in c3.CORE_DEFINING)


@pytest.mark.parametrize("a1,a2,pheno", [
    ("*1", "*1", "Normal Metabolizer"),       # expressor
    ("*1", "*3", "Intermediate Metabolizer"),
    ("*1", "*6", "Intermediate Metabolizer"),
    ("*1", "*7", "Intermediate Metabolizer"),
    ("*3", "*3", "Poor Metabolizer"),         # non-expressor
    ("*6", "*7", "Poor Metabolizer"),
    ("*3", "*6", "Poor Metabolizer"),
    ("*7", "*7", "Poor Metabolizer"),
])
def test_c3_function_pair_phenotype(a1, a2, pheno):
    assert c3.diplotype_phenotype(a1, a2) == pheno
    assert c3.diplotype_phenotype(a2, a1) == pheno


def test_c3_unknown_allele_indeterminate():
    assert c3.diplotype_phenotype("*3", "*99") == "Indeterminate"


# --- caller (VCF) — incl. the *7 INSERTION (the novel indel case) ---
def test_c3_star3_hom_nonexpressor(tmp_path):
    r = _call_c3(_c3_vcf(tmp_path, [("*3", "1/1")]))
    assert r.diplotype == "*3/*3" and r.phenotype == "Poor Metabolizer"


def test_c3_star7_insertion_het(tmp_path):
    # the *7 insertion (REF=T ALT=TA) must be matched natively by the caller
    r = _call_c3(_c3_vcf(tmp_path, [("*7", "0/1")]))
    assert r.diplotype == "*1/*7" and r.phenotype == "Intermediate Metabolizer"


def test_c3_star7_hom_insertion(tmp_path):
    r = _call_c3(_c3_vcf(tmp_path, [("*7", "1/1")]))
    assert r.diplotype == "*7/*7" and r.phenotype == "Poor Metabolizer"


def test_c3_compound_star6_star7(tmp_path):
    r = _call_c3(_c3_vcf(tmp_path, [("*6", "0/1"), ("*7", "0/1")]))
    assert set(r.diplotype.split("/")) == {"*6", "*7"} and r.phenotype == "Poor Metabolizer"


def test_c3_all_ref_expressor(tmp_path):
    r = _call_c3(_c3_vcf(tmp_path, [("*3", "0/0"), ("*6", "0/0"), ("*7", "0/0")]))
    assert r.diplotype == "*1/*1" and r.phenotype == "Normal Metabolizer"


# --- runner record + CLI ---
def test_c3_runner_record(tmp_path):
    rec = call_cyp3a5(_c3_vcf(tmp_path, [("*3", "0/1")]), sample_id="T1")
    assert rec["gene"] == "CYP3A5"
    assert rec["diplotype"] == "*1/*3"
    assert rec["phenotype"] == "Intermediate Metabolizer"
    assert rec["phenotype_abbrev"] == "IM"
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True
    assert rec["caller"]["calling_independently_validatable"] is True


def test_c3_cli_routing(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_c3_vcf(tmp_path, [("*7", "1/1")])), "--gene", "cyp3a5", "--json-only"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "CYP3A5" in out and "Poor Metabolizer" in out


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
