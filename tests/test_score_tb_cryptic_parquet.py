"""Offline tests for the CRyPTIC-parquet TB baseline adapter (pure parser + label loader + match wiring)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.score_tb_cryptic_parquet import (  # noqa: E402
    load_labels, parse_cryptic_variant,
)
from dna_decode.data.tb_who_catalogue import Determinant  # noqa: E402
from dna_decode.organism_rules import tb_amr, tb_vcf  # noqa: E402


def test_parse_cryptic_variant_snv():
    assert parse_cryptic_variant("761155c>t") == (761155, "C", "T")
    assert parse_cryptic_variant("2155168c>g") == (2155168, "C", "G")
    assert parse_cryptic_variant(" 761139c>t ") == (761139, "C", "T")


def test_parse_cryptic_variant_rejects_nonsnv():
    assert parse_cryptic_variant("rpoB_p.Asp435Phe") is None      # protein form
    assert parse_cryptic_variant("761155_ins_g") is None          # indel
    assert parse_cryptic_variant("761155del") is None
    assert parse_cryptic_variant("") is None and parse_cryptic_variant(None) is None


def test_parsed_snv_matches_genomic_determinant():
    """A parsed CRyPTIC SNV feeds VariantCall -> score_drug genomic match against a WHO determinant."""
    pos, ref, alt = parse_cryptic_variant("761155c>t")
    calls = {pos: tb_vcf.VariantCall(pos=pos, ref=ref, alt=alt, gt="1/1")}
    det = Determinant(drug="Rifampicin", gene="rpoB", variant="rpoB_p.Ser450Leu", grade="1) Assoc w R",
                      tier="1", chrom="NC_000962.3", pos=761155, ref="C", alt="T")
    assert tb_amr.score_drug("rifampicin", calls, [det]).prediction == "R"
    # an isolate without the determinant SNV -> S (callability unassessed, no regeno)
    assert tb_amr.score_drug("rifampicin", {}, [det]).prediction == "S"


def test_codon_mnv_determinant_matched_via_per_base_snvs():
    """A 2-base codon determinant matches an isolate carrying the constituent single-base CRyPTIC rows."""
    det = Determinant(drug="Rifampicin", gene="rpoB", variant="rpoB_p.Asp435Phe", grade="1) Assoc w R",
                      tier="1", chrom="NC_000962.3", pos=761154, ref="GA", alt="TT")  # 2-base MNV
    # isolate carries both single-base components (as separate CRyPTIC per-position rows would)
    calls = {761154: tb_vcf.VariantCall(761154, "G", "T", "1/1"),
             761155: tb_vcf.VariantCall(761155, "A", "T", "1/1")}
    assert tb_amr.score_drug("rifampicin", calls, [det]).prediction == "R"
    # carrying only ONE component -> no match (S)
    assert tb_amr.score_drug("rifampicin", {761154: tb_vcf.VariantCall(761154, "G", "T", "1/1")},
                             [det]).prediction == "S"


def test_load_labels(tmp_path):
    p = tmp_path / "reuse.csv"
    p.write_text(
        "UNIQUEID,RIF_BINARY_PHENOTYPE,RIF_PHENOTYPE_QUALITY,INH_BINARY_PHENOTYPE,INH_PHENOTYPE_QUALITY\n"
        "iso1,R,HIGH,S,HIGH\n"
        "iso2,S,LOW,R,HIGH\n"            # RIF quality LOW -> dropped for RIF
        "iso3,,HIGH,R,HIGH\n"           # RIF blank -> dropped for RIF
        "iso4,R,HIGH,U,HIGH\n",         # INH 'U' -> dropped for INH
        encoding="utf-8")
    rif = load_labels(p, "RIF")
    assert rif == {"iso1": "R", "iso4": "R"}        # iso2 LOW, iso3 blank excluded
    inh = load_labels(p, "INH")
    assert inh == {"iso1": "S", "iso2": "R", "iso3": "R"}   # iso4 'U' excluded


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
