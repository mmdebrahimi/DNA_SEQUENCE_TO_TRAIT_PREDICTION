"""Offline tests for the pigment-cell VCF input path. Synthetic VCF fixtures; no network/no D:.

Pins: rsID extraction from a VCF, GT->genotype mapping (phased/unphased, REF/ALT indexing), missing/uncallable
site omission, indel skip, and the CLI --vcf path reproducing the inline-genotype result.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pigment import predict_eye_color  # noqa: E402
from dna_decode.pigment.vcf_input import genotypes_from_vcf  # noqa: E402
from dna_decode.pigment.cli import main as pigment_main  # noqa: E402

# A minimal VCF: all 6 IrisPlex SNPs. HERC2 rs12913832 GG (blue-associated); others set to counted-allele-absent.
# counted alleles: rs12913832=A rs1800407=T rs12896399=T rs16891982=C rs1393350=A rs12203592=T.
_VCF_BLUE = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
15\t28120472\trs12913832\tA\tG\t.\tPASS\t.\tGT\t1/1
15\t27985172\trs1800407\tC\tT\t.\tPASS\t.\tGT\t0/0
14\t92307319\trs12896399\tG\tT\t.\tPASS\t.\tGT\t0/0
5\t33951588\trs16891982\tC\tG\t.\tPASS\t.\tGT\t1|1
11\t89277878\trs1393350\tG\tA\t.\tPASS\t.\tGT\t0/0
6\t396321\trs12203592\tC\tT\t.\tPASS\t.\tGT\t0/0
"""


def _write(tmp_path, text, name="x.vcf"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_extract_all_six(tmp_path):
    g = genotypes_from_vcf(_write(tmp_path, _VCF_BLUE))
    assert set(g) == {"rs12913832", "rs1800407", "rs12896399", "rs16891982", "rs1393350", "rs12203592"}
    assert g["rs12913832"] == "GG"          # 1/1 with REF=A,ALT=G -> GG (0 counted A)
    assert g["rs1800407"] == "CC"           # 0/0 with REF=C -> CC


def test_vcf_predicts_blue(tmp_path):
    g = genotypes_from_vcf(_write(tmp_path, _VCF_BLUE))
    res = predict_eye_color(g)
    assert res.call == "blue" and res.counted_alleles["rs12913832"] == 0


def test_phased_and_unphased(tmp_path):
    g = genotypes_from_vcf(_write(tmp_path, _VCF_BLUE))
    assert g["rs16891982"] == "GG"          # 1|1 phased with REF=C,ALT=G


def test_missing_and_indel_sites_omitted(tmp_path):
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "15\t1\trs12913832\tA\tG\t.\tPASS\t.\tGT\t./.\n"        # missing -> omitted
        "15\t2\trs1800407\tC\tCTT\t.\tPASS\t.\tGT\t0/1\n"       # indel -> omitted
        "6\t3\trs12203592\tC\tT\t.\tPASS\t.\tGT\t0/1\n"         # callable SNV -> kept
    )
    g = genotypes_from_vcf(_write(tmp_path, vcf))
    assert "rs12913832" not in g and "rs1800407" not in g
    assert g["rs12203592"] == "CT"


def test_cli_vcf_path(tmp_path, capsys):
    rc = pigment_main(["--vcf", _write(tmp_path, _VCF_BLUE), "--json"])
    assert rc == 0
    import json
    assert json.loads(capsys.readouterr().out)["call"] == "blue"


def test_cli_vcf_missing_file(capsys):
    rc = pigment_main(["--vcf", "does_not_exist_12345.vcf"])
    assert rc == 2 and "error" in capsys.readouterr().err


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
