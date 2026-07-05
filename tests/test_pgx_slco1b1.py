"""Pins the SLCO1B1 PGx cell — single-SNP rs4149056 (c.521T>C, *5) statin-myopathy function readout.

Honest tier = KNOWLEDGE_BASELINE (single-SNP genotype->function, like VKORC1), NOT an independent
star-diplotype number. Plus-strand: genomic T>C == cDNA 521T>C (no flip).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import PGX_GENES  # noqa: E402
from dna_decode.pgx.slco1b1 import ALT, CHROM, POS, REF, RSID, call_slco1b1  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")


def _vcf(tmp_path, gt, name="s.vcf"):
    p = tmp_path / name
    p.write_text(_HEADER + f"chr{CHROM}\t{POS}\t{RSID}\t{REF}\t{ALT}\t.\tPASS\t.\tGT\t{gt}\n", encoding="utf-8")
    return p


def test_slco1b1_in_pgx_genes():
    assert "slco1b1" in PGX_GENES


def test_slco1b1_coords_grounded():
    # rs4149056 GRCh38 chr12:21178615 T>C, plus-strand (Ensembl-verified)
    assert (CHROM, POS, REF, ALT) == ("12", 21178615, "T", "C")


def test_slco1b1_ref_normal(tmp_path):
    r = call_slco1b1(_vcf(tmp_path, "0|0"))
    assert r["alt_count"] == 0 and r["variant_genotype"] == "T/T"
    assert r["star_proxy"] == "*1/*1" and r["function"] == "Normal Function"
    assert r["myopathy_risk"] == "typical_risk"


def test_slco1b1_het_decreased(tmp_path):
    r = call_slco1b1(_vcf(tmp_path, "0|1"))
    assert r["variant_genotype"] == "T/C" and r["star_proxy"] == "*1/*5"
    assert r["function"] == "Decreased Function" and r["myopathy_risk"] == "intermediate_risk"


def test_slco1b1_homalt_poor(tmp_path):
    r = call_slco1b1(_vcf(tmp_path, "1|1"))
    assert r["variant_genotype"] == "C/C" and r["star_proxy"] == "*5/*5"
    assert r["function"] == "Poor Function" and r["myopathy_risk"] == "high_risk"


def test_slco1b1_absent_assumed_reference(tmp_path):
    p = tmp_path / "empty.vcf"
    p.write_text(_HEADER + "chr1\t999\trsX\tA\tG\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    r = call_slco1b1(p)
    assert r["status"] == "assumed_reference"
    assert "assumed_reference_at_uncalled_site" in r["flags"]
    assert r["function"] == "Normal Function"


def test_slco1b1_cli_routing(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_vcf(tmp_path, "1|1")), "--gene", "slco1b1"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "SLCO1B1" in out and "Poor Function" in out and "high_risk" in out


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
