"""Offline pins for the PGP-UK targeted-liftover PGx realizer (`scripts/pgx_decode_pgp_uk.py`).

The real end-to-end run is network-dependent (proven on 3 real PGP-UK individuals 2026-07-05,
wiki/pgx_pgp_uk_realization_2026-07-05.md). These tests pin the pure logic: the GRCh37 map covers every
catalog site, and the GRCh38 mini-VCF builder handles lift / allele-mismatch / assumed-ref / uncovered-chrom
correctly.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import pgx_decode_pgp_uk as m  # noqa: E402


def test_grch37_map_covers_every_catalog_site():
    """A missing GRCh37 position would silently drop a gene on real GRCh37 VCFs."""
    missing = [rs for rs, *_ in m._all_variants() if rs not in m.GRCH37_POS or m.GRCH37_POS[rs] is None]
    assert not missing, f"GRCh37 position missing for: {missing}"


def test_all_variants_nonempty_and_spans_genes():
    chroms = {chrom for _rs, chrom, *_ in m._all_variants()}
    # CYP2C cluster (10), CYP3A5 (7), CYP2B6 (19), TPMT (6), VKORC1 (16), SLCO1B1 (12)
    assert {"6", "7", "10", "12", "16", "19"} <= chroms


def test_minivcf_lifts_concordant_snp():
    # found rs4244285 het; catalog GRCh38 10:94781859 G>A
    found = {"rs4244285": ("G", "A", "0/1")}
    vcf, audit = m.build_grch38_minivcf(found, chroms_seen={"10"}, sample_id="S")
    assert "94781859\trs4244285\tG\tA" in vcf and "GT\t0/1" in vcf
    assert any(a["rsid"] == "rs4244285" and a["state"] == "lifted" for a in audit)


def test_minivcf_absent_on_covered_chrom_is_assumed_ref():
    # chrom 10 covered but no variant found -> the CYP2C19 sites get 0/0 (assumed ref)
    vcf, audit = m.build_grch38_minivcf({}, chroms_seen={"10"}, sample_id="S")
    assert "rs4244285\tG\tA\t.\tPASS\t.\tGT\t0/0" in vcf
    assert any(a["rsid"] == "rs4244285" and a["state"] == "assumed_ref_0/0_wgs" for a in audit)


def test_minivcf_absent_on_uncovered_chrom_is_omitted_not_ref():
    # no chromosome seen -> never fabricate a ref call; the site is omitted (caller -> honest no_input)
    vcf, audit = m.build_grch38_minivcf({}, chroms_seen=set(), sample_id="S")
    assert "rs4244285" not in vcf
    assert any(a["rsid"] == "rs4244285" and a["state"] == "chrom_not_covered_omitted" for a in audit)


def test_minivcf_allele_mismatch_dropped():
    # a strand/allele mismatch must NOT be lifted (would flip the GT meaning)
    found = {"rs4244285": ("G", "C", "0/1")}   # catalog alt is A, VCF alt is C
    vcf, audit = m.build_grch38_minivcf(found, chroms_seen={"10"}, sample_id="S")
    assert "rs4244285\tG\tA" not in vcf
    assert any(a["rsid"] == "rs4244285" and a["state"] == "allele_mismatch" for a in audit)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
