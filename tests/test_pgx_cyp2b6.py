"""Pins the CYP2B6 PGx cell — *6-proxy from the 516G>T signal (efavirenz), single-SNP v0.

Validated 62/62 on clean *1/*6 truth vs the GeT-RM CDC consolidated consensus (committed artifact; VCF
gitignored). Honest scope: rs2279343 (785A>G, the 2nd *6 component) is absent from the 1000G 30x panel,
so v0 is a single-SNP *6-proxy that cannot split *6 from *9 (documented).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import PGX_GENES  # noqa: E402
from dna_decode.pgx import cyp2b6_catalog as c6  # noqa: E402
from dna_decode.pgx.caller import call_diplotype  # noqa: E402
from dna_decode.pgx.runner import call_cyp2b6  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
_POS, _REF, _ALT, _RSID = 41006936, "G", "T", "rs3745274"


def _vcf(tmp_path, gt, name="c6.vcf"):
    p = tmp_path / name
    p.write_text(_HEADER + f"chr19\t{_POS}\t{_RSID}\t{_REF}\t{_ALT}\t.\tPASS\t.\tGT\t{gt}\n", encoding="utf-8")
    return p


def _call(vcf, sample=None):
    return call_diplotype(vcf, sample=sample, defining=c6.CORE_DEFINING, sentinels=c6.SENTINELS,
                          reference_allele=c6.REFERENCE_ALLELE, phenotype_fn=c6.diplotype_phenotype, gene=c6.GENE)


def test_cyp2b6_in_pgx_genes():
    assert "cyp2b6" in PGX_GENES


def test_cyp2b6_coord_grounded():
    d = c6.CORE_DEFINING[0]
    assert (d.star, d.rsid, d.chrom, d.pos, d.ref, d.alt) == ("*6", "rs3745274", "19", 41006936, "G", "T")


def test_cyp2b6_star6_het(tmp_path):
    r = _call(_vcf(tmp_path, "0|1"))
    assert r.diplotype == "*1/*6" and r.phenotype == "Intermediate Metabolizer"


def test_cyp2b6_star6_hom(tmp_path):
    r = _call(_vcf(tmp_path, "1|1"))
    assert r.diplotype == "*6/*6" and r.phenotype == "Poor Metabolizer"


def test_cyp2b6_ref_normal(tmp_path):
    r = _call(_vcf(tmp_path, "0|0"))
    assert r.diplotype == "*1/*1" and r.phenotype == "Normal Metabolizer"


def test_cyp2b6_runner_record_flags_single_snp_proxy(tmp_path):
    rec = call_cyp2b6(_vcf(tmp_path, "0|1"), sample_id="E1")
    assert rec["gene"] == "CYP2B6"
    assert rec["diplotype"] == "*1/*6"
    assert rec["caller"]["is_single_snp_proxy"] is True   # honesty: not the full 2-SNP compound
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True
    # the *9 blind-spot is documented in the undetectable list
    assert any("785" in u or "star9" in u for u in rec["undetectable"])


def test_cyp2b6_cli(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_vcf(tmp_path, "1|1")), "--gene", "cyp2b6", "--json-only"])
    assert rc == 3 or rc == 0
    assert "CYP2B6" in capsys.readouterr().out


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
