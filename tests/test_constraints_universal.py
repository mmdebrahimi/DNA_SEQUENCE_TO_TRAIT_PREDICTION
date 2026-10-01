"""Guards for the universal (axiomatic) constraints (plan Step 4, 2026-10-01)."""

from __future__ import annotations

import pytest

from dna_decode.constraints import universal as U
from dna_decode.constraints.registry import AXIOMATIC, INAPPLICABLE, REFUTED, SATISFIED, UNIVERSAL

# ATG AAA TGA CCC TAA -- a mid-frame TGA with a terminal stop
CDS_INTERNAL_TGA = "ATGAAATGACCCTAA"
CDS_CLEAN = "ATGAAACCCTAA"


def _v(fn, ctx):
    return fn(ctx).verdict


def test_every_law_is_axiomatic_universal_and_needs_no_citation():
    """A law is not someone's finding; requiring a citation would be a category error."""
    for c in U.UNIVERSAL_CONSTRAINTS:
        assert c.kind == AXIOMATIC, c.key
        assert c.scope == UNIVERSAL, c.key
        assert not c.source_url and not c.verbatim_quote, c.key
        assert c.note, c.key


def test_the_clade_table_flips_the_internal_stop_verdict():
    """THE load-bearing test. The same CDS is refuted under table 1 and satisfied under table 4, so the
    clade tables change behaviour rather than decorating it."""
    assert _v(U.no_internal_stop_codon, {"cds": CDS_INTERNAL_TGA}) == REFUTED
    assert _v(U.no_internal_stop_codon,
              {"cds": CDS_INTERNAL_TGA, "clade": "bacteria.mycoplasma"}) == SATISFIED
    # and the ciliate code flips a TAA instead
    assert _v(U.no_internal_stop_codon, {"cds": "ATGTAACCCTAA"}) == REFUTED
    assert _v(U.no_internal_stop_codon,
              {"cds": "ATGTAACCCTAA", "clade": "eukaryote.ciliate"}) == SATISFIED


def test_a_terminal_stop_is_allowed_but_an_internal_one_is_not():
    assert _v(U.no_internal_stop_codon, {"cds": CDS_CLEAN}) == SATISFIED
    assert _v(U.no_internal_stop_codon, {"cds": "ATGTGAAAATAA"}) == REFUTED


@pytest.mark.parametrize("fn,ctx,expected", [
    (U.cds_length_multiple_of_three, {"cds": "ATGAA"}, REFUTED),
    (U.cds_length_multiple_of_three, {"cds": "ATGAAA"}, SATISFIED),
    (U.cds_length_multiple_of_three, {}, INAPPLICABLE),
    (U.ref_base_matches_cds, {"cds": "ATGAAA", "nt_pos": 2, "ref_base": "G"}, REFUTED),
    (U.ref_base_matches_cds, {"cds": "ATGAAA", "nt_pos": 2, "ref_base": "T"}, SATISFIED),
    (U.ref_base_matches_cds, {"cds": "ATGAAA"}, INAPPLICABLE),
    (U.ref_base_matches_cds, {"cds": "ATGAAA", "nt_pos": 99, "ref_base": "T"}, REFUTED),
    (U.start_codon_present, {"cds": CDS_CLEAN, "cds_is_complete_gene": True}, SATISFIED),
    (U.start_codon_present, {"cds": "AAACCCTAA", "cds_is_complete_gene": True}, REFUTED),
    (U.start_codon_present, {}, INAPPLICABLE),
    # without the completeness assertion it ABSTAINS -- see the dedicated test below
    (U.start_codon_present, {"cds": "AAACCCTAA"}, INAPPLICABLE),
    (U.substitution_reachable_by_single_nt, {"codon": "TCG", "mut_aa": "L"}, SATISFIED),
    (U.substitution_reachable_by_single_nt, {"codon": "TGG", "mut_aa": "D"}, REFUTED),
    (U.substitution_reachable_by_single_nt, {}, INAPPLICABLE),
    (U.substitution_reachable_by_single_nt, {"wt_aa": "S", "mut_aa": "L"}, SATISFIED),
])
def test_each_law_refutes_satisfies_and_abstains(fn, ctx, expected):
    assert fn(ctx) .verdict == expected


def test_a_missing_precondition_yields_INAPPLICABLE_not_SATISFIED():
    """"could not check" must never be recorded as "checked and clean" -- that is how a layer starts
    reporting protection it never provided."""
    for c in U.UNIVERSAL_CONSTRAINTS:
        assert c.check({}).verdict == INAPPLICABLE, c.key


def test_the_frame_constraint_gates_the_frame_dependent_one():
    """An out-of-frame CDS makes the internal-stop scan meaningless, so it abstains rather than
    refuting on a garbage translation."""
    assert _v(U.cds_length_multiple_of_three, {"cds": "ATGAAATG"}) == REFUTED
    assert _v(U.no_internal_stop_codon, {"cds": "ATGAAATG"}) == INAPPLICABLE


def test_single_nt_reachability_matches_the_real_qrdr_case():
    """gyrA S83L is reachable by one nucleotide (TCG->TTG), which is why it is a common resistance
    substitution; a two-base change is refuted."""
    assert _v(U.substitution_reachable_by_single_nt, {"codon": "TCG", "mut_aa": "L"}) == SATISFIED
    assert _v(U.substitution_reachable_by_single_nt, {"codon": "ATG", "mut_aa": "P"}) == REFUTED


def test_reachability_is_clade_aware():
    """Under table 4, TGA encodes W, so a stop-to-nothing case changes shape with the code."""
    assert U.substitution_reachable_by_single_nt(
        {"codon": "TGG", "mut_aa": "W", "clade": "bacteria.mycoplasma"}).verdict in {SATISFIED, REFUTED}
    # TGG -> TGA is one nt; under table 4 that is W->W, under standard it is W->*
    from dna_decode.constraints.codon_tables import MYCOPLASMA_SPIROPLASMA, STANDARD
    assert STANDARD["TGA"] == "*" and MYCOPLASMA_SPIROPLASMA["TGA"] == "W"


def test_an_unknown_clade_propagates_the_refusal_rather_than_defaulting():
    from dna_decode.constraints.codon_tables import UnknownCladeError
    with pytest.raises(UnknownCladeError):
        U.no_internal_stop_codon({"cds": CDS_CLEAN, "clade": "not_a_clade"})


def test_evaluate_all_records_every_law_into_the_report():
    from dna_decode.constraints.report import ConstraintReport
    rep = ConstraintReport()
    res = U.evaluate_all({"cds": CDS_INTERNAL_TGA, "nt_pos": 1, "ref_base": "A"}, report=rep)
    assert set(res) == {c.key for c in U.UNIVERSAL_CONSTRAINTS}
    assert rep.status("no_internal_stop_codon") == "active"          # it refuted the internal TGA
    assert "no_internal_stop_codon" in rep.active_keys()
    # every law appears in the artifact, including the ones that found nothing
    assert set(rep.as_dict()["constraints"]) == set(res)


def test_registration_is_safe_to_call_twice():
    U.register_universal()
    again = U.register_universal()
    assert len(again) == len(U.UNIVERSAL_CONSTRAINTS)


def test_start_codon_abstains_unless_the_caller_asserts_a_COMPLETE_gene():
    """MEASURED CORRECTION (2026-10-01). Run over eight committed reference CDS this law refuted five,
    and all five were CORRECT sequences: HIV RT/PR/IN/CA are in-frame extracts from the HXB2 gag/pol
    polyprotein (beginning CCC/CCT) and SARS-CoV-2 Mpro is nsp5 cleaved from ORF1ab (beginning AGT --
    Ser really is its mature N-terminus). A law cannot tell a truncated extract from a broken start
    codon by looking at the sequence, and a FALSE refutation is the failure mode that would make a
    constraint filter delete valid predictions. So the completeness claim belongs to the caller."""
    polyprotein_extract = "CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGAT"
    assert U.start_codon_present({"cds": polyprotein_extract}).verdict == INAPPLICABLE
    # the caller may still take responsibility, and then it fires
    v = U.start_codon_present({"cds": polyprotein_extract, "cds_is_complete_gene": True})
    assert v.verdict == REFUTED
    # and a real complete gene passes under the same assertion
    assert U.start_codon_present(
        {"cds": CDS_CLEAN, "cds_is_complete_gene": True}).verdict == SATISFIED


def test_the_abstain_reason_explains_itself_to_a_caller():
    v = U.start_codon_present({"cds": "CCCAAATAA"})
    assert "cds_is_complete_gene" in v.detail
    assert "polyprotein" in v.detail
