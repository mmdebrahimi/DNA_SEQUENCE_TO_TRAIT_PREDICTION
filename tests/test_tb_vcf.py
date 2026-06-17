"""Step 1 — TB VCF genotype-aware parser + callability (brainstorm C1 + C3)."""
from __future__ import annotations

from dna_decode.organism_rules import tb_vcf

_FMT = "GT:DP:DPF:COV:FRS:GT_CONF:GT_CONF_PERCENTILE"
_HDR = (
    "##fileformat=VCFv4.2\n"
    '##FILTER=<ID=MIN_GCP,Description="Minimum GT_CONF_PERCENTILE of 5.0">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def _rec(pos, ref, alt, filt, gt, dp=200, frs="1.0", gcp="78.0"):
    samp = f"{gt}:{dp}:1.1:0,{dp}:{frs}:1500:{gcp}"
    return f"{tb_vcf.CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\t{_FMT}\t{samp}\n"


def _masked(*recs):
    return _HDR + "".join(recs)


# --- C1: the call rule ---------------------------------------------------------------------------

def test_pass_nonref_call_parses_with_provenance():
    calls = tb_vcf.parse_masked_calls(_masked(_rec(761155, "C", "T", "PASS", "1/1", frs="0.99")))
    assert 761155 in calls
    c = calls[761155]
    assert (c.ref, c.alt, c.gt) == ("C", "T", "1/1")
    assert c.frs == 0.99 and c.dp == 200 and c.gt_conf_percentile == 78.0
    # provenance only — there is no GCP attribute / FORMAT field
    assert not hasattr(c, "gcp")


def test_ref_call_0_0_with_alt_is_excluded():
    calls = tb_vcf.parse_masked_calls(_masked(_rec(100, "A", "G", "PASS", "0/0")))
    assert calls == {}


def test_nocall_excluded():
    calls = tb_vcf.parse_masked_calls(_masked(_rec(100, "A", "G", "PASS", "./.")))
    assert calls == {}


def test_filter_failed_excluded():
    calls = tb_vcf.parse_masked_calls(_masked(_rec(761155, "C", "T", "MIN_GCP", "1/1")))
    assert calls == {}


def test_multiallelic_resolved_by_gt_index():
    # katG S315T canonical: 2155168 C>G is ALT index 2 of A,G,T -> GT 2/2.
    calls = tb_vcf.parse_masked_calls(_masked(_rec(2155168, "C", "A,G,T", "PASS", "2/2")))
    assert calls[2155168].alt == "G"


def test_no_gcp_format_key_is_read():
    # Regression guard: the FORMAT carries no GCP key; parsing must not assume one.
    assert "GCP" not in _FMT.split(":")
    calls = tb_vcf.parse_masked_calls(_masked(_rec(761155, "C", "T", "PASS", "1/1")))
    assert 761155 in calls  # parsed fine without any GCP lookup


# --- C3: callability from the regeno VCF ---------------------------------------------------------

def test_callable_positions_ref_alt_nocall_absent():
    regeno = _masked(
        _rec(761155, "C", "T", "PASS", "1/1"),   # alt call -> callable
        _rec(2155168, "C", ".", "PASS", "0/0"),  # explicit ref -> callable
        _rec(500, "A", "G", "PASS", "./."),       # no-call -> uncallable
    )
    flags = tb_vcf.callable_positions(regeno, [761155, 2155168, 500, 999])
    assert flags[761155] is True
    assert flags[2155168] is True
    assert flags[500] is False     # ./.
    assert flags[999] is False     # absent from regeno


def test_window_callable_all_vs_one_nocall():
    ok = _masked(_rec(10, "A", ".", "PASS", "0/0"), _rec(11, "C", ".", "PASS", "0/0"),
                 _rec(12, "G", ".", "PASS", "0/0"))
    assert tb_vcf.is_window_callable(ok, 10, 12) is True
    bad = _masked(_rec(10, "A", ".", "PASS", "0/0"), _rec(11, "C", "T", "PASS", "./."),
                  _rec(12, "G", ".", "PASS", "0/0"))
    assert tb_vcf.is_window_callable(bad, 10, 12) is False


def test_window_invalid_raises():
    import pytest
    with pytest.raises(ValueError):
        tb_vcf.is_window_callable(_HDR, 12, 10)


def test_vcf_paths_for_reads_both_columns():
    row = {"VCF": " ../x.masked.vcf.gz ", "REGENOTYPED_VCF": "../x.regeno.vcf.gz"}
    assert tb_vcf.vcf_paths_for(row) == ("../x.masked.vcf.gz", "../x.regeno.vcf.gz")
