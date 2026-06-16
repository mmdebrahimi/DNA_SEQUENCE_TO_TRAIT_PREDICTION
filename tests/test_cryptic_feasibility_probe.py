"""Pin the CRyPTIC feasibility probe's pure logic (census, RRDR VCF parse, VCF URL construction)."""
from scripts.cryptic_feasibility_probe import (
    RRDR_HI, RRDR_LO, _vcf_url, census, clonality_proxy, rrdr_variant_present,
)


def _row(rif_call, rif_qual="HIGH", inh_call="S", inh_qual="HIGH", uid="site.05.subj.1"):
    r = {c: "NA" for d in ("AMI", "BDQ", "CFZ", "DLM", "EMB", "ETH", "INH", "KAN", "LEV", "LZD", "MXF", "RIF", "RFB")
         for c in (f"{d}_BINARY_PHENOTYPE", f"{d}_PHENOTYPE_QUALITY")}
    r["RIF_BINARY_PHENOTYPE"], r["RIF_PHENOTYPE_QUALITY"] = rif_call, rif_qual
    r["INH_BINARY_PHENOTYPE"], r["INH_PHENOTYPE_QUALITY"] = inh_call, inh_qual
    r["UNIQUEID"] = uid
    return r


def test_census_counts_high_quality_only():
    rows = [_row("R", "HIGH"), _row("R", "LOW"), _row("S", "HIGH"), _row("S", "HIGH")]
    cen = census(rows)
    assert cen["rifampicin"]["R_all"] == 2 and cen["rifampicin"]["R_high"] == 1
    assert cen["rifampicin"]["S_high"] == 2


def test_census_powering_verdict():
    rows = [_row("R", "HIGH") for _ in range(20)] + [_row("S", "HIGH") for _ in range(20)]
    assert census(rows)["rifampicin"]["verdict"] == "POWERED"
    assert census(rows[:21])["rifampicin"]["verdict"] == "UNDERPOWERED"  # 20R/1S -> S short


def _vcf(records):
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    body = "".join(f"NC_000962.3\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n" for pos, ref, alt in records)
    return (head + body).encode("utf-8")


def test_rrdr_variant_present_in_window():
    assert rrdr_variant_present(_vcf([(RRDR_LO + 50, "C", "T")])) is True


def test_rrdr_variant_absent_outside_window():
    assert rrdr_variant_present(_vcf([(RRDR_LO - 100, "C", "T"), (RRDR_HI + 100, "G", "A")])) is False


def test_rrdr_ignores_non_alt():
    assert rrdr_variant_present(_vcf([(RRDR_LO + 10, "C", ".")])) is False


def test_vcf_url_strips_relative_reproducibility_prefix():
    url = _vcf_url("../reproducibility/00/01/08/61/10861/site.02.masked.vcf.gz")
    assert url == ("https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/"
                   "00/01/08/61/10861/site.02.masked.vcf.gz")


def test_clonality_proxy_counts_sites():
    rows = [_row("R", uid="site.05.subj.1"), _row("R", uid="site.05.subj.2"), _row("S", uid="site.10.subj.1")]
    c = clonality_proxy(rows)
    assert c["distinct_sites"] == 2 and c["top_sites"]["05"] == 2
