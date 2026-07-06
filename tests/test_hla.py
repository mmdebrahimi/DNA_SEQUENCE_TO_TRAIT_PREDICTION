"""Pins the HLA drug-hypersensitivity cell (dna-hla) — tag-SNP carriage via LD proxy.

The new substrate after the PGx cells (user-ratified 2026-07-06). v0 anchor = HLA-B*57:01 / abacavir
(rs2395029, the gold-standard clinical tag). These tests pin the NCBI-verified GRCh38 tag coords, the
carriage logic, and the LOAD-BEARING honesty rail: the tag is an LD PROXY (is_ld_proxy=True), NEVER
sequence-based typing, and the sample-level concordance vs real HLA truth is the SCORED number (not the
literature LD alone).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.hla import HLA_ALLELES  # noqa: E402
from dna_decode.hla.caller import call_hla  # noqa: E402
from dna_decode.hla.catalog import CATALOG, get  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")


def _vcf(tmp_path, key, gt, name="hla.vcf"):
    a = get(key)
    p = tmp_path / name
    p.write_text(_HEADER + f"chr{a.chrom}\t{a.pos}\t{a.rsid}\t{a.ref}\t{a.tag_alt}\t.\tPASS\t.\tGT\t{gt}\n",
                 encoding="utf-8")
    return p


def test_alleles_declared():
    assert set(HLA_ALLELES) == set(CATALOG)
    assert "b5701" in HLA_ALLELES


def test_b5701_tag_coords_grounded_ncbi():
    a = get("b5701")
    assert (a.rsid, a.chrom, a.pos, a.ref, a.tag_alt) == ("rs2395029", "6", 31464003, "T", "G")
    assert a.drug == "abacavir" and a.proxy_tier == "gold_standard"


def test_b5801_a3101_provisional():
    assert get("b5801").proxy_tier == "provisional" and get("b5801").drug == "allopurinol"
    assert get("a3101").proxy_tier == "provisional" and get("a3101").drug == "carbamazepine"


@pytest.mark.parametrize("gt,copies,carrier,zyg", [
    ("0|0", 0, False, "non-carrier"),
    ("0|1", 1, True, "heterozygous carrier"),
    ("1|1", 2, True, "homozygous carrier"),
])
def test_carriage_calls(tmp_path, gt, copies, carrier, zyg):
    rec = call_hla(_vcf(tmp_path, "b5701", gt), "b5701", sample="SAMPLE")
    assert rec["tag_copies"] == copies and rec["carrier"] is carrier and rec["zygosity"] == zyg
    assert rec["caller"]["is_ld_proxy"] is True          # LOAD-BEARING: never sequence-based typing
    if carrier:
        assert "AVOID abacavir" in rec["risk_call"]


def test_absent_record_is_noncarrier_flagged(tmp_path):
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(_HEADER + "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n", encoding="utf-8")
    rec = call_hla(vcf, "b5701", sample="SAMPLE")
    assert rec["carrier"] is False and rec["status"] == "assumed_reference"
    assert "assumed_reference_at_uncalled_site" in rec["flags"]


def test_caveat_names_proxy_and_pending_truth(tmp_path):
    rec = call_hla(_vcf(tmp_path, "b5701", "0|1"), "b5701", sample="SAMPLE")
    assert "LD PROXY" in rec["caveat"] and "NOT sequence-based" in rec["caveat"]
    assert "NOT a clinical tool" in rec["caveat"]


def test_cli_runs(tmp_path, capsys):
    from dna_decode.hla.cli import main
    rc = main([str(_vcf(tmp_path, "b5701", "1|1")), "--allele", "b5701", "--sample", "SAMPLE"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "HLA-B*57:01" in out and "abacavir" in out and "AVOID abacavir" in out


def test_cli_missing_vcf_exits_2(tmp_path):
    from dna_decode.hla.cli import main
    assert main([str(tmp_path / "nope.vcf")]) == 2


# --- registry integration: the HLA cells are on the trust surface ---
def test_hla_cells_registered():
    from dna_decode.data.cell_registry import cli_routable_manifest, hla_cells
    targets = {c.target for c in hla_cells()}
    assert targets == set(HLA_ALLELES) == cli_routable_manifest()["dna-hla"]
    assert all(c.track == "hla" and c.route == "dna-hla" for c in hla_cells())


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
