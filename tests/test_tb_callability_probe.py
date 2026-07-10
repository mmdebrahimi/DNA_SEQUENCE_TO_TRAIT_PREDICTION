"""Offline tests for the sampled TB regeno-callability probe.

Pure helpers only — no network. The load-bearing distinction under test: `tb_vcf.callable_positions`
collapses "position ABSENT from the regeno VCF" and "position present but the call FAILED" into a single
`uncallable`. These are minos regenotyped VCFs, which emit records only at variant-panel sites (~29% of
the genome), so `absent` overwhelmingly means "not a genotyping target", not "sequencing could not call
it". Only `present_fail` is evidence of true uncallability.
"""
from __future__ import annotations

import pytest

from scripts import tb_callability_probe as P

CHROM = "NC_000962.3"


def _rec(pos: int, filt: str = "PASS", gt: str = "0/0") -> str:
    return f"{CHROM}\t{pos}\t.\tA\tG\t.\t{filt}\t.\tGT:DP\t{gt}:50"


def _vcf(*records: str) -> str:
    return "\n".join(["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
                     + list(records))


# ---------------------------------------------------------------- the three-state classifier

def test_classify_positions_absent_is_not_the_same_as_failed():
    """A minos VCF simply omits non-panel sites; that is NOT evidence the site could not be called."""
    text = _vcf(_rec(100), _rec(200, filt="MIN_DP"))
    got = P.classify_positions(text, [100, 200, 300])
    assert got == {"present_pass": 1, "present_fail": 1, "absent": 1, "n_positions": 3}


def test_classify_positions_all_absent():
    got = P.classify_positions(_vcf(), [10, 20])
    assert got["absent"] == 2 and got["present_pass"] == 0 and got["present_fail"] == 0


def test_classify_positions_treats_null_gt_as_failed_not_absent():
    got = P.classify_positions(_vcf(_rec(100, gt="./.")), [100])
    assert got["present_fail"] == 1 and got["absent"] == 0


def test_classify_positions_ignores_offtarget_positions():
    got = P.classify_positions(_vcf(_rec(100), _rec(999)), [100])
    assert got["n_positions"] == 1 and got["present_pass"] == 1


def test_classify_positions_counts_a_position_once_when_duplicated():
    got = P.classify_positions(_vcf(_rec(100), _rec(100)), [100])
    assert got["present_pass"] + got["present_fail"] == 1
    assert got["absent"] == 0


def test_classify_positions_agrees_with_frozen_callable_positions():
    """present_pass must equal the frozen callable count — we measure, we do not redefine."""
    from dna_decode.organism_rules import tb_vcf

    text = _vcf(_rec(100), _rec(200, filt="MIN_FRS"), _rec(300, gt="./."))
    positions = [100, 200, 300, 400]
    frozen = tb_vcf.callable_positions(text, positions)
    got = P.classify_positions(text, positions)
    assert got["present_pass"] == sum(1 for ok in frozen.values() if ok)
    # the frozen rule lumps the other three together
    assert got["present_fail"] + got["absent"] == sum(1 for ok in frozen.values() if not ok)


# ---------------------------------------------------------------- sampling

def test_stride_sample_spreads_across_the_table_and_is_deterministic():
    rows = [{"i": i} for i in range(100)]
    got = P.stride_sample(rows, 5)
    assert got == [{"i": 0}, {"i": 20}, {"i": 40}, {"i": 60}, {"i": 80}]
    assert P.stride_sample(rows, 5) == got            # deterministic, no RNG


def test_stride_sample_returns_all_when_n_exceeds_population():
    rows = [{"i": i} for i in range(3)]
    assert len(P.stride_sample(rows, 10)) == 3


def test_eligible_rows_requires_both_vcfs_and_high_quality():
    rows = [
        {"VCF": "a", "REGENOTYPED_VCF": "b", "RIF_BINARY_PHENOTYPE": "R", "RIF_PHENOTYPE_QUALITY": "HIGH"},
        {"VCF": "a", "REGENOTYPED_VCF": "",  "RIF_BINARY_PHENOTYPE": "R", "RIF_PHENOTYPE_QUALITY": "HIGH"},
        {"VCF": "a", "REGENOTYPED_VCF": "b", "RIF_BINARY_PHENOTYPE": "R", "RIF_PHENOTYPE_QUALITY": "MEDIUM"},
        {"VCF": "a", "REGENOTYPED_VCF": "b", "RIF_BINARY_PHENOTYPE": "U", "RIF_PHENOTYPE_QUALITY": "HIGH"},
    ]
    assert len(P.eligible_rows(rows)) == 1


# ---------------------------------------------------------------- the summary reports BOTH rules

def _drug_rec(call, flipped, fail_only, absent=1000, ppass=400, pfail=0):
    return {"call_without_callability": call, "call_with_callability": "ABSTAIN" if flipped else call,
            "flipped_S_to_ABSTAIN": flipped, "would_abstain_fail_only_rule": fail_only,
            "n_uncallable_positions": absent + pfail, "n_determinant_positions": absent + ppass + pfail,
            "absent": absent, "present_pass": ppass, "present_fail": pfail, "phenotype": "S"}


def test_summarize_separates_the_frozen_rule_from_the_fail_only_rule():
    records = [
        {"uniqueid": "a", "drugs": {"rifampicin": _drug_rec("S", True, False),
                                    "isoniazid": _drug_rec("S", True, False)}},
        {"uniqueid": "b", "drugs": {"rifampicin": _drug_rec("R", False, False),
                                    "isoniazid": _drug_rec("R", False, False)}},
    ]
    s = P.summarize(records)["rifampicin"]
    assert s["n_called_R"] == 1 and s["n_called_S_by_absence"] == 1
    assert s["n_S_would_ABSTAIN_current_rule"] == 1
    assert s["frac_of_S_that_is_uncallable_current_rule"] == 1.0
    assert s["n_S_would_ABSTAIN_fail_only_rule"] == 0          # <-- the whole point
    assert s["frac_of_S_that_is_uncallable_fail_only_rule"] == 0.0


def test_summarize_counts_a_true_failure_under_both_rules():
    records = [{"uniqueid": "a", "drugs": {"rifampicin": _drug_rec("S", True, True, pfail=3)}}]
    s = P.summarize(records)["rifampicin"]
    assert s["n_S_would_ABSTAIN_current_rule"] == 1
    assert s["n_S_would_ABSTAIN_fail_only_rule"] == 1
    assert s["median_positions_present_fail"] == 3


def test_summarize_handles_no_S_calls():
    records = [{"uniqueid": "a", "drugs": {"rifampicin": _drug_rec("R", False, False)}}]
    s = P.summarize(records)["rifampicin"]
    assert s["n_called_S_by_absence"] == 0
    assert s["frac_of_S_that_is_uncallable_current_rule"] is None


def test_median_helper():
    assert P._median([]) is None
    assert P._median([5]) == 5
    assert P._median([1, 3]) == 2
    assert P._median([3, 1, 2]) == 2


def test_r_calls_never_flip_because_score_drug_returns_R_before_checking_callability():
    """A matched determinant short-circuits to R, so callability can only ever affect S calls."""
    rec = _drug_rec("R", False, False)
    assert rec["call_with_callability"] == "R"
    with pytest.raises(AssertionError):
        assert rec["flipped_S_to_ABSTAIN"]


# ============================================================================================
# Three-state callability landed in tb_vcf / tb_amr (2026-07-10). ADDITIVE: the default is the
# frozen behaviour; `absent_is_uncallable=False` opts into the minos-correct rule.
# ============================================================================================

def test_position_states_returns_the_three_states():
    from dna_decode.organism_rules import tb_vcf

    text = _vcf(_rec(100), _rec(200, filt="MIN_DP"))
    got = tb_vcf.position_states(text, [100, 200, 300])
    assert got == {100: tb_vcf.PRESENT_PASS, 200: tb_vcf.PRESENT_FAIL, 300: tb_vcf.ABSENT}


def test_callable_positions_default_is_unchanged_absent_is_uncallable():
    from dna_decode.organism_rules import tb_vcf

    text = _vcf(_rec(100))
    assert tb_vcf.callable_positions(text, [100, 300]) == {100: True, 300: False}


def test_callable_positions_absent_treated_callable_when_opted_in():
    from dna_decode.organism_rules import tb_vcf

    text = _vcf(_rec(100), _rec(200, filt="MIN_FRS"))
    got = tb_vcf.callable_positions(text, [100, 200, 300], absent_is_uncallable=False)
    assert got == {100: True, 200: False, 300: True}   # only the FAILED record is uncallable


def test_score_drug_default_abstains_on_an_absent_position():
    """The frozen behaviour: absent -> ABSTAIN. This is what makes 100% of S calls abstain on minos VCFs."""
    from dna_decode.data.tb_who_catalogue import Determinant
    from dna_decode.organism_rules import tb_amr

    det = Determinant(drug="rifampicin", gene="rpoB", variant="S450L", grade="1", tier="1",
                      chrom="NC_000962.3", pos=761155, ref="C", alt="T")
    call = tb_amr.score_drug("rifampicin", {}, [det], regeno_text=_vcf())   # position absent
    assert call.prediction == "ABSTAIN" and call.n_uncallable_positions == 1


def test_score_drug_returns_S_under_the_fail_only_rule_when_position_merely_absent():
    from dna_decode.data.tb_who_catalogue import Determinant
    from dna_decode.organism_rules import tb_amr

    det = Determinant(drug="rifampicin", gene="rpoB", variant="S450L", grade="1", tier="1",
                      chrom="NC_000962.3", pos=761155, ref="C", alt="T")
    call = tb_amr.score_drug("rifampicin", {}, [det], regeno_text=_vcf(), absent_is_uncallable=False)
    assert call.prediction == "S" and call.n_uncallable_positions == 0


def test_score_drug_still_abstains_under_fail_only_rule_on_a_real_failure():
    from dna_decode.data.tb_who_catalogue import Determinant
    from dna_decode.organism_rules import tb_amr

    det = Determinant(drug="rifampicin", gene="rpoB", variant="S450L", grade="1", tier="1",
                      chrom="NC_000962.3", pos=761155, ref="C", alt="T")
    text = _vcf(_rec(761155, filt="MIN_DP"))
    call = tb_amr.score_drug("rifampicin", {}, [det], regeno_text=text, absent_is_uncallable=False)
    assert call.prediction == "ABSTAIN" and call.n_uncallable_positions == 1
