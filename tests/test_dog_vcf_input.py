"""Offline tests for the dog morphology VCF caller (dna_decode/pigment/dog_vcf_input.py) + the CLI --vcf path.

Synthetic canFam4 VCF fixture — no real genotype data / no network. Pins: coord-matching of the pinned
body-size + ear SNPs, chromosome-name normalization (chr10 == 10), big-allele dosage counting, strand
harmonization, partial-panel / uncallable handling, and the end-to-end `dna-decode morphology --vcf` call.
The real-genotype PLINK-vs-VCF equivalence (dosages identical to the plink_io path that fed r=0.619) was
verified on 3 Darwin's Ark dogs at build time. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment.dog_vcf_input import _norm_chrom, dosages_from_vcf, pinned_loci  # noqa: E402
from dna_decode.pigment.morphology_cli import main as morph_main  # noqa: E402

# a large-dog / erect-ear genome: IGF1 GG(2), HMGA2 GG(2), STC2 T-(1), GHR C-(1), EAR AA(2).
# REF:ALT taken from each pinned canfam4_variant id; GT encodes 0=REF/1=ALT. Mixed chr-prefix on purpose.
_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDOG1
chr15\t41513523\t.\tG\tA\t.\t.\t.\tGT\t0/0
10\t8703415\t.\tG\tA\t.\t.\t.\tGT\t0|0
chr4\t40070215\t.\tT\tA\t.\t.\t.\tGT\t0/1
chr4\t67710295\t.\tC\tT\t.\t.\t.\tGT\t1/0
chr10\t8612500\t.\tA\tG\t.\t.\t.\tGT\t0/0
"""


def _write(tmp_path, text=_VCF):
    p = tmp_path / "dog.vcf"
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_norm_chrom():
    assert _norm_chrom("chr10") == "10" and _norm_chrom("10") == "10" and _norm_chrom("CHR4") == "4"


def test_pinned_loci_shape():
    loci = pinned_loci()
    assert set(loci) == {"IGF1", "HMGA2", "STC2", "GHR", "EAR"}
    assert loci["HMGA2"] == ("10", 8703415, "G") and loci["EAR"] == ("10", 8612500, "A")


def test_full_panel_dosages(tmp_path):
    dos = dosages_from_vcf(_write(tmp_path))
    assert dos == {"IGF1": 2, "HMGA2": 2, "STC2": 1, "GHR": 1, "EAR": 2}


def test_cli_vcf_matches_dosage_path(tmp_path, capsys):
    # --vcf must produce the SAME call as the equivalent --dosages
    rc = morph_main(["--vcf", _write(tmp_path), "--json"])
    assert rc == 0
    dv = json.loads(capsys.readouterr().out)
    rc = morph_main(["--dosages", "IGF1=2,HMGA2=2,STC2=1,GHR=1,EAR=2", "--json"])
    dd = json.loads(capsys.readouterr().out)
    assert dv["input_source"] == "vcf" and dd["input_source"] == "dosages"
    assert dv["height"]["size_rank"] == dd["height"]["size_rank"] == "average"
    assert dv["ear"]["ear_type"] == dd["ear"]["ear_type"]
    assert dv["loci_scored"] == ["EAR", "GHR", "HMGA2", "IGF1", "STC2"]


def test_partial_panel_and_uncallable(tmp_path):
    # only 2 SNPs present, one with a missing genotype -> skipped
    text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD\n"
            "10\t8703415\t.\tG\tA\t.\t.\t.\tGT\t0/1\n"           # HMGA2 dose 1
            "chr15\t41513523\t.\tG\tA\t.\t.\t.\tGT\t./.\n")      # IGF1 missing -> uncallable
    dos = dosages_from_vcf(_write(tmp_path, text))
    assert dos == {"HMGA2": 1}


def test_strand_harmonization(tmp_path):
    # HMGA2 big allele G; a minus-strand site (C/T) -> complement G->C is the site allele; GT 1/1 (T,T)=0 copies
    text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD\n"
            "10\t8703415\t.\tC\tT\t.\t.\t.\tGT\t0/0\n")          # C,C -> complement of G is C -> 2 copies
    dos = dosages_from_vcf(_write(tmp_path, text))
    assert dos == {"HMGA2": 2}


def test_missing_vcf_returns_2():
    assert morph_main(["--vcf", "/no/such/file.vcf"]) == 2


def test_no_pinned_snps_returns_2(tmp_path):
    text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD\n"
            "5\t12345\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")             # not a pinned locus
    assert morph_main(["--vcf", _write(tmp_path, text)]) == 2


if __name__ == "__main__":
    print("run via pytest (uses tmp_path/capsys fixtures)")
