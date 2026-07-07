"""Offline pins for the ClinVar/Mendelian per-VCF decoder wrapper (`scripts/clinvar_decode_vcf.py`).

The real run is on committed panels + real PGP-UK VCFs (proven 2026-07-06,
wiki/clinvar_pgp_uk_realization_2026-07-06.md). These tests pin the pure logic on a synthetic panel + VCF:
only CARRIED alts are decoded, the pathogenic + benign paths both fire, and not-in-panel is abstained.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from dna_decode.data.clinvar import ClinVarDecoder  # noqa: E402
import clinvar_decode_vcf as m  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")


def _decoder():
    return ClinVarDecoder({
        ("7", "100", "A", "T"): {"significance": "Pathogenic", "review_status": "reviewed_by_expert_panel",
                                 "gene": "CFTR", "disease": "Cystic fibrosis", "clinvar_id": "1"},
        ("11", "200", "C", "G"): {"significance": "Benign", "review_status": "criteria_provided,_single_submitter",
                                  "gene": "HBB", "disease": "beta-thal", "clinvar_id": "2"},
    })


def test_carried_alts():
    assert m._carried_alts("A", "T", "0/1") == ["T"]
    assert m._carried_alts("A", "T,G", "1/2") == ["T", "G"]
    assert m._carried_alts("A", "T", "0/0") == []      # ref/ref carries nothing
    assert m._carried_alts("A", "T", "./.") == []      # no-call


def test_decode_reports_pathogenic_and_benign_carried(tmp_path):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(_HEADER
                   + "7\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"    # carries the pathogenic CFTR alt
                   + "11\t200\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\n"   # carries the benign HBB alt (hom)
                   + "7\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\n",   # ref -> not carried (dedup-safe: different line)
                   encoding="utf-8")
    rep = m.decode_vcf(vcf, _decoder(), sample="S1")
    assert rep["n_pathogenic"] == 1 and rep["n_benign"] == 1
    assert rep["pathogenic_hits"][0]["gene"] == "CFTR" and rep["pathogenic_hits"][0]["stars"] == 3


def test_not_in_panel_is_indeterminate_not_a_hit(tmp_path):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(_HEADER + "1\t999\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    rep = m.decode_vcf(vcf, _decoder(), sample="S1")
    assert rep["n_pathogenic"] == 0 and rep["n_benign"] == 0
    assert rep["n_indeterminate_not_in_panel"] == 1   # honest abstain, never a fabricated benign


def test_ref_only_variant_not_decoded(tmp_path):
    # a panel-matching site the person is REF at must not be reported (they don't carry it)
    vcf = tmp_path / "s.vcf"
    vcf.write_text(_HEADER + "7\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\n", encoding="utf-8")
    rep = m.decode_vcf(vcf, _decoder(), sample="S1")
    assert rep["n_pathogenic"] == 0 and rep["n_indeterminate_not_in_panel"] == 0


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
