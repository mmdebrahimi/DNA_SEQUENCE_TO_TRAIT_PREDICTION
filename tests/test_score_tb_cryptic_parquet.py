"""Offline tests for the CRyPTIC-parquet TB baseline adapter (pure parser + label loader + match wiring)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.score_tb_cryptic_parquet import (  # noqa: E402
    DRUG_CODE, _resolve_drugs, _score_drug_from_calls, indel_determinant_targets, load_labels,
    parse_cryptic_variant, who_indel_to_cryptic_string,
)
from dna_decode.data.tb_who_catalogue import DRUG_CATALOGUE_NAME, Determinant  # noqa: E402
from dna_decode.organism_rules import tb_amr, tb_vcf  # noqa: E402


# --- second-line extension (2026-07-14): DRUG_CODE now spans first- + second-line + new drugs ------
def test_drug_code_covers_second_line():
    # the second-line + new-drug reuse-table codes are wired (the only thing that gated second-line)
    expected = {"moxifloxacin": "MXF", "levofloxacin": "LEV", "amikacin": "AMI", "kanamycin": "KAN",
                "ethionamide": "ETH", "bedaquiline": "BDQ", "clofazimine": "CFZ", "delamanid": "DLM",
                "linezolid": "LZD", "ethambutol": "EMB"}
    for drug, code in expected.items():
        assert DRUG_CODE.get(drug) == code, f"{drug} should map to reuse code {code}"
    # first-line still intact
    assert DRUG_CODE["rifampicin"] == "RIF" and DRUG_CODE["isoniazid"] == "INH"


def test_drug_code_keys_are_catalogue_names():
    # every DRUG_CODE key MUST be a DRUG_CATALOGUE_NAME key, else load_determinants(drug) KeyErrors
    for drug in DRUG_CODE:
        assert drug in DRUG_CATALOGUE_NAME, f"{drug} is in DRUG_CODE but not DRUG_CATALOGUE_NAME"


def test_reuse_codes_are_uppercase_2to3_letter():
    # reuse-table phenotype-column prefixes are 2-3 uppercase letters (RIF/INH/MXF/BDQ/...)
    for code in DRUG_CODE.values():
        assert code.isupper() and 2 <= len(code) <= 3, f"unexpected reuse code shape: {code!r}"


# --- single-parquet-pass multi-drug mode (2026-07-14): score many drugs in ONE stream ---------------
def test_resolve_drugs():
    assert _resolve_drugs("all") == list(DRUG_CODE)
    # 'all-remaining' excludes the 4 already-scored (RIF/INH/MXF/BDQ)
    rem = _resolve_drugs("all-remaining")
    assert set(rem) == set(DRUG_CODE) - {"rifampicin", "isoniazid", "moxifloxacin", "bedaquiline"}
    assert "moxifloxacin" not in rem and "levofloxacin" in rem
    assert _resolve_drugs("levofloxacin,amikacin") == ["levofloxacin", "amikacin"]
    assert _resolve_drugs(" linezolid , ethambutol ") == ["linezolid", "ethambutol"]
    import pytest
    with pytest.raises(SystemExit):
        _resolve_drugs("not_a_drug")


def test_multi_matches_single_scoring_from_shared_calls():
    """The load-bearing multi-drug claim: scoring drug A from a SHARED calls+indel_hits set (carrying
    OTHER drugs' positions + indel strings) is BYTE-IDENTICAL to scoring A from a clean per-drug set —
    foreign positions are filtered, foreign indel strings are intersected out."""
    det_snv = Determinant("Rifampicin", "rpoB", "rpoB_p.Ser450Leu", "1) Assoc w R", "1",
                          "NC_000962.3", 100, "C", "T")
    det_indel = Determinant("Rifampicin", "rpoB", "rpoB_p.Xdup", "1) Assoc w R", "1",
                            "NC_000962.3", 200, "A", "ATTC")   # -> "200_ins_ttc"
    dets = [det_snv, det_indel]
    indel_dets = [det_indel]
    targets = indel_determinant_targets(dets)                  # {"200_ins_ttc": det_indel}
    assert "200_ins_ttc" in targets
    labels = {"iso1": "R", "iso2": "S", "iso3": "R"}
    # one shared lineage SNP (pos 500) so all 3 isolates form ONE non-singleton cluster (else score_cohort
    # returns the BLOCKED-no-lineage dict). pos 500 is a barcode position -> survives the wanted-filter.
    from dna_decode.data.tb_lineage_barcode import BarcodeSNP
    barcode = [BarcodeSNP(pos=500, lineage="lineage4", allele="G")]
    bc = lambda: tb_vcf.VariantCall(500, "A", "G", "1/1")      # noqa: E731

    # SHARED: iso1 has A's SNV(100) + a FOREIGN pos 300; iso3 relies on A's indel + a FOREIGN indel string
    calls_shared = {
        "iso1": {100: tb_vcf.VariantCall(100, "C", "T", "1/1"), 300: tb_vcf.VariantCall(300, "C", "G", "1/1"),
                 500: bc()},
        "iso2": {300: tb_vcf.VariantCall(300, "C", "G", "1/1"), 500: bc()},
        "iso3": {500: bc()},
    }
    indel_hits_shared = {"iso3": {"200_ins_ttc", "400_ins_g"}, "iso1": {"400_ins_g"}}

    # CLEAN: only A's positions + the barcode pos / only A's indel string
    calls_clean = {"iso1": {100: tb_vcf.VariantCall(100, "C", "T", "1/1"), 500: bc()},
                   "iso2": {500: bc()}, "iso3": {500: bc()}}
    indel_hits_clean = {"iso3": {"200_ins_ttc"}}

    kw = dict(drug="rifampicin", dets=dets, indel_dets=indel_dets, targets=targets, n_complex=0,
              barcode=barcode, labels=labels, max_isolates=0, match_indels=True)
    res_shared = _score_drug_from_calls(calls=calls_shared, indel_hits=indel_hits_shared, **kw)
    res_clean = _score_drug_from_calls(calls=calls_clean, indel_hits=indel_hits_clean, **kw)

    assert res_shared == res_clean                             # full equivalence
    # sanity: iso1 R via SNV, iso3 R via indel, iso2 S -> perfect on this toy set
    assert res_clean["raw"]["sens"] == 1.0 and res_clean["raw"]["spec"] == 1.0
    # foreign indel "400_ins_g" did NOT flip iso1 (it is not one of A's targets)
    assert res_shared["indel_matching"]["isolates_with_indel_hit"] == 1     # only iso3
    # all 3 carry the barcode pos 500 -> all in-scope; the FOREIGN pos 300 was filtered (proven by ==)
    assert res_shared["n_with_inscope_calls"] == 3


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


def test_who_indel_to_cryptic_string_verified_examples():
    """EXACT conversion verified against real CRyPTIC VARIANTS rows (scripts/_tb_indel_probe.py)."""
    # rpoB Phe433dup insertion: 761102 A>ATTC -> 761102_ins_ttc (NOT 761103 — ins anchors at pos+len(ref)-1)
    assert who_indel_to_cryptic_string(761102, "A", "ATTC") == "761102_ins_ttc"
    # rpoB Thr444dup: 761135 G>GACC -> 761135_ins_acc
    assert who_indel_to_cryptic_string(761135, "G", "GACC") == "761135_ins_acc"
    # rpoB deletion: 761100 CAATTCATGG>C -> 761101_del_aattcatgg (del anchors at pos+len(alt))
    assert who_indel_to_cryptic_string(761100, "CAATTCATGG", "C") == "761101_del_aattcatgg"
    # katG single-base deletion: 2153894 GC>G -> 2153895_del_c
    assert who_indel_to_cryptic_string(2153894, "GC", "G") == "2153895_del_c"


def test_who_indel_to_cryptic_string_complex_returns_none():
    # alt does not share a clean prefix/suffix with ref (true delins) -> unmapped (needs ref left-alignment)
    assert who_indel_to_cryptic_string(100, "A", "CCAGCATT") is None     # ins not anchored on ref base
    assert who_indel_to_cryptic_string(100, "ACGT", "TGCA") is None      # equal-len would be MNV not indel here
    # an SNV (equal length) is not an indel target either
    assert who_indel_to_cryptic_string(761155, "C", "T") is None


def test_indel_determinant_targets_skips_snv_and_complex():
    dets = [
        Determinant("Rifampicin", "rpoB", "rpoB_p.Ser450Leu", "1) Assoc w R", "1", "NC_000962.3",
                    761155, "C", "T"),                                   # SNV -> skipped
        Determinant("Rifampicin", "rpoB", "rpoB_p.Phe433dup", "1) Assoc w R", "1", "NC_000962.3",
                    761102, "A", "ATTC"),                                # ins -> 761102_ins_ttc
        Determinant("Isoniazid", "katG", "katG_LoF", "1) Assoc w R", "1", "NC_000962.3",
                    2153894, "GC", "G"),                                 # del -> 2153895_del_c
        Determinant("Isoniazid", "katG", "katG_complex", "1) Assoc w R", "1", "NC_000962.3",
                    100, "A", "CCAGCATT"),                               # complex -> skipped
    ]
    targets = indel_determinant_targets(dets)
    assert set(targets) == {"761102_ins_ttc", "2153895_del_c"}
    assert targets["761102_ins_ttc"].variant == "rpoB_p.Phe433dup"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
