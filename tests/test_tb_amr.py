"""Step 3 — TB decoder cell (RIF + INH) v1a plumbing + the coordinate-alignment fixture."""
from __future__ import annotations

import pytest

from dna_decode.data import tb_who_catalogue as cat
from dna_decode.data.tb_who_catalogue import Determinant
from dna_decode.organism_rules import tb_amr, tb_vcf

_FMT = "GT:DP:DPF:COV:FRS:GT_CONF:GT_CONF_PERCENTILE"
_HDR = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"


def _rec(pos, ref, alt, gt, filt="PASS"):
    return f"{tb_vcf.CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\t{_FMT}\t{gt}:200:1.1:0,200:1.0:1500:78\n"


def _vcf(*recs):
    return _HDR + "".join(recs)


# synthetic grade-1/2 determinants (the two canonical sentinels)
RPOB_S450L = Determinant("Rifampicin", "rpoB", "rpoB_p.Ser450Leu", "1) Assoc w R", "1",
                         tb_vcf.CHROM, 761155, "C", "T")
KATG_S315T = Determinant("Isoniazid", "katG", "katG_p.Ser315Thr", "1) Assoc w R", "1",
                         tb_vcf.CHROM, 2155168, "C", "G")
INHA_15 = Determinant("Isoniazid", "inhA", "inhA_c.-15C>T", "2) Assoc w R - Interim", "2",
                      tb_vcf.CHROM, 1673425, "C", "T")


def test_rif_determinant_match_is_R():
    calls = tb_vcf.parse_masked_calls(_vcf(_rec(761155, "C", "T", "1/1")))
    out = tb_amr.score_rif(calls, [RPOB_S450L])
    assert out.prediction == tb_amr.R
    assert out.matched[0].variant == "rpoB_p.Ser450Leu"
    assert out.rule_status == "KNOWLEDGE_BASELINE" and out.input_type == "vcf_h37rv"


def test_inh_determinant_match_is_R_and_scope_multilocus():
    calls = tb_vcf.parse_masked_calls(_vcf(_rec(2155168, "C", "G", "1/1")))
    out = tb_amr.score_inh(calls, [KATG_S315T, INHA_15])
    assert out.prediction == tb_amr.R
    assert out.coverage_scope == ("inhA", "katG")  # ratified A — multi-locus scope


def test_determinant_matched_inside_mnv():
    # isolate carries S450L (761155 C>T) inside a 3-base MNV record 761154 CCG>CTG (the dominant FN cause)
    calls = tb_vcf.parse_masked_calls(_vcf(_rec(761154, "CCG", "CTG", "1/1")))
    out = tb_amr.score_rif(calls, [RPOB_S450L])
    assert out.prediction == tb_amr.R
    assert out.matched[0].variant == "rpoB_p.Ser450Leu"


def test_no_determinant_all_callable_is_S():
    calls = {}  # no variant
    regeno = _vcf(_rec(761155, "C", ".", "0/0"))  # determinant position callable (explicit ref)
    out = tb_amr.score_rif(calls, [RPOB_S450L], regeno_text=regeno)
    assert out.prediction == tb_amr.S
    assert out.n_uncallable_positions == 0 and out.callability_assessed is True


def test_no_determinant_uncallable_window_is_ABSTAIN():
    calls = {}
    regeno = _vcf(_rec(761155, "C", "T", "./."))  # determinant position no-call -> uncallable
    out = tb_amr.score_rif(calls, [RPOB_S450L], regeno_text=regeno)
    assert out.prediction == tb_amr.ABSTAIN
    assert out.n_uncallable_positions == 1  # C3: never susceptible-by-absence


def test_no_regeno_is_S_with_callability_unassessed():
    out = tb_amr.score_rif({}, [RPOB_S450L])
    assert out.prediction == tb_amr.S and out.callability_assessed is False


def test_wrong_alt_at_determinant_position_does_not_match():
    # a different ALT at 761155 (C>A, not the S450L C>T) must NOT fire the determinant
    calls = tb_vcf.parse_masked_calls(_vcf(_rec(761155, "C", "A", "1/1")))
    out = tb_amr.score_rif(calls, [RPOB_S450L], regeno_text=_vcf(_rec(761155, "C", "A", "1/1")))
    assert out.prediction != tb_amr.R


# ---- coordinate-alignment fixture (real catalogue; skips if absent) ------------------------------

_REAL = pytest.mark.skipif(not cat.catalogue_available(),
                           reason="pinned WHO catalogue not present (gitignored)")


@_REAL
def test_alignment_fixture_rpoB_S450L_resolves_to_R():
    dets = cat.load_determinants("rifampicin")
    calls = tb_vcf.parse_masked_calls(_vcf(_rec(761155, "C", "T", "1/1")))
    out = tb_amr.score_rif(calls, dets)
    assert out.prediction == tb_amr.R
    assert any(d.variant == "rpoB_p.Ser450Leu" for d in out.matched)


@_REAL
def test_alignment_fixture_katG_S315T_resolves_to_R():
    dets = cat.load_determinants("isoniazid")
    calls = tb_vcf.parse_masked_calls(_vcf(_rec(2155168, "C", "G", "1/1")))
    out = tb_amr.score_inh(calls, dets)
    assert out.prediction == tb_amr.R
    assert any(d.variant == "katG_p.Ser315Thr" for d in out.matched)
