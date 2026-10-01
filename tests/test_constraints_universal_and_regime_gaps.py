"""Gap coverage for `constraints/universal.py` and `constraints/specializations.py`
(test-epilogue, 2026-10-01).

The plan's files cover each law's happy path, the clade flip, and the start-codon abstention. What they
leave is the ABSTENTION SURFACE -- the branches that decide `inapplicable` rather than `refuted`. That is
the half that matters most for a filter: a law that refutes when it should abstain deletes valid
predictions, which the probe already caught once with `start_codon_present` refuting five correct
polyprotein extracts. Every remaining `_na(...)` branch is pinned here.

Two further findings, both recorded rather than patched:

* `_table(ctx)` falls back to the STANDARD code whenever `clade` is falsy. Everywhere else this module
  refuses a default; here it quietly supplies one, so an absent clade key reads as table 1.
* `register_universal()` / `register_specializations()` wrap `register` in a bare `except Exception` and
  return the constraint anyway. That is intended as idempotence for a re-run, but it also swallows a
  genuine `UnsourcedConstraintError` -- the provenance gate's own registration helper fails OPEN, and
  the caller cannot tell "registered" from "refused".
"""

from __future__ import annotations

import pytest

from dna_decode.constraints import registry as R
from dna_decode.constraints import specializations as S
from dna_decode.constraints import universal as U
from dna_decode.constraints.registry import (
    AXIOMATIC,
    INAPPLICABLE,
    REFUTED,
    SATISFIED,
    SOURCED,
    UNIVERSAL,
    Constraint,
    ConstraintVerdict,
)

# ===================================================================== abstention surface

@pytest.mark.parametrize("ctx,why", [
    ({}, "nothing at all"),
    ({"wt_aa": "S"}, "wt_aa without mut_aa"),
    ({"mut_aa": "L"}, "mut_aa without wt_aa"),
    ({"codon": "TCG"}, "codon without mut_aa"),
    ({"codon": "NNN", "mut_aa": "L"}, "codon not in the table"),
    ({"codon": "TC", "mut_aa": "L"}, "codon is not a triplet"),
    ({"codon": "TCGA", "mut_aa": "L"}, "codon is four bases"),
    ({"wt_aa": "X", "mut_aa": "L"}, "no codon encodes wt_aa under this code"),
    ({"wt_aa": "B", "mut_aa": "L"}, "ambiguity letter is not an amino acid here"),
])
def test_reachability_ABSTAINS_on_every_incomplete_input_rather_than_refuting(ctx, why):
    """A missing precondition must never read as 'needs more than one nucleotide change'. That verdict
    would delete a prediction the law never actually examined."""
    v = U.substitution_reachable_by_single_nt(ctx)
    assert v.verdict == INAPPLICABLE, f"{why}: got {v.verdict} ({v.detail})"
    assert v.detail, "an abstention must say what was missing"


def test_reachability_still_REFUTES_a_genuine_two_nucleotide_hop():
    """Non-vacuity for the block above: the law is not simply abstaining on everything."""
    assert U.substitution_reachable_by_single_nt(
        {"codon": "ATG", "mut_aa": "P"}).verdict == REFUTED
    assert U.substitution_reachable_by_single_nt(
        {"wt_aa": "W", "mut_aa": "D"}).verdict == REFUTED
    assert U.substitution_reachable_by_single_nt(
        {"wt_aa": "S", "mut_aa": "L"}).verdict == SATISFIED


def test_reachability_treats_a_STOP_as_an_ordinary_target_amino_acid():
    """`mut_aa='*'` is how a nonsense prediction arrives; it must be scored, not abstained on."""
    assert U.substitution_reachable_by_single_nt({"codon": "TCA", "mut_aa": "*"}).verdict == SATISFIED
    assert U.substitution_reachable_by_single_nt({"codon": "ATG", "mut_aa": "*"}).verdict == REFUTED


def test_the_codon_path_WINS_over_wt_aa_and_does_not_cross_check_them():
    """Known limit: when both are supplied the codon decides and no consistency check runs, so a
    caller passing a mismatched pair gets the codon's answer silently."""
    ctx = {"codon": "TCG", "wt_aa": "W", "mut_aa": "L"}       # TCG is Ser, wt_aa says Trp
    assert U.substitution_reachable_by_single_nt(ctx).verdict == SATISFIED
    assert "TCG" in U.substitution_reachable_by_single_nt(ctx).detail


@pytest.mark.parametrize("ctx", [
    {},
    {"cds": "ATGTCG"},
    {"cds": "ATGTCG", "nt_pos": 1},
    {"cds": "ATGTCG", "ref_base": "A"},
    {"cds": "", "nt_pos": 1, "ref_base": "A"},
    {"cds": "ATGTCG", "nt_pos": 1, "ref_base": ""},
])
def test_ref_base_check_ABSTAINS_without_all_three_inputs(ctx):
    assert U.ref_base_matches_cds(ctx).verdict == INAPPLICABLE


@pytest.mark.parametrize("pos", [0, -1, 7, 1000])
def test_an_out_of_range_coordinate_is_REFUTED_not_abstained(pos):
    """The opposite call from a missing input: the caller gave a coordinate and it is wrong. That is
    exactly the frame/offset error this law exists to catch, so it must refute."""
    v = U.ref_base_matches_cds({"cds": "ATGTCG", "nt_pos": pos, "ref_base": "A"})
    assert v.verdict == REFUTED and "outside the cds" in v.detail


def test_the_ref_base_check_is_case_insensitive_on_both_sides():
    for cds, ref in (("atgtcg", "a"), ("ATGTCG", "a"), ("atgtcg", "A")):
        assert U.ref_base_matches_cds({"cds": cds, "nt_pos": 1, "ref_base": ref}).verdict == SATISFIED
    assert U.ref_base_matches_cds(
        {"cds": "ATGTCG", "nt_pos": 1, "ref_base": "C"}).verdict == REFUTED


@pytest.mark.parametrize("ctx", [
    {},
    {"cds": ""},
    {"cds": "AT", "cds_is_complete_gene": True},
    {"cds": "ATGTCG"},                                   # completeness not asserted
    {"cds": "ATGTCG", "cds_is_complete_gene": False},
])
def test_start_codon_ABSTAINS_on_a_short_or_uncommitted_sequence(ctx):
    assert U.start_codon_present(ctx).verdict == INAPPLICABLE


def test_start_codon_refutes_ONLY_when_completeness_is_asserted_and_the_codon_is_wrong():
    assert U.start_codon_present(
        {"cds": "ATGTCG", "cds_is_complete_gene": True}).verdict == SATISFIED
    v = U.start_codon_present({"cds": "CCCTCG", "cds_is_complete_gene": True})
    assert v.verdict == REFUTED and "CCC" in v.detail


def test_the_frame_law_abstains_on_an_empty_cds_rather_than_calling_zero_in_frame():
    """`len('') % 3 == 0`, so a naive implementation would report an empty sequence as in frame."""
    assert U.cds_length_multiple_of_three({}).verdict == INAPPLICABLE
    assert U.cds_length_multiple_of_three({"cds": ""}).verdict == INAPPLICABLE
    assert U.cds_length_multiple_of_three({"cds": "ATGG"}).verdict == REFUTED
    assert U.cds_length_multiple_of_three({"cds": "ATG"}).verdict == SATISFIED


def test_the_internal_stop_law_abstains_on_an_OUT_OF_FRAME_cds_rather_than_refuting():
    """Ordering matters: reading a frameshifted sequence in triplets manufactures stops that are not
    there, so this law defers to the frame law instead of reporting them."""
    v = U.no_internal_stop_codon({"cds": "ATGTGAA"})
    assert v.verdict == INAPPLICABLE and "frame" in v.detail


def test_a_LONE_stop_codon_is_accepted_because_the_terminal_codon_is_exempt():
    """`codons[:-1]` means a single-codon CDS has no internal positions at all. Pinned so the exemption
    is a stated scope rather than an accident."""
    assert U.no_internal_stop_codon({"cds": "TGA"}).verdict == SATISFIED
    assert U.no_internal_stop_codon({"cds": "ATGTGA"}).verdict == SATISFIED
    assert U.no_internal_stop_codon({"cds": "TGAATG"}).verdict == REFUTED


def test_an_unresolvable_codon_is_IGNORED_by_the_internal_stop_law():
    """`tab.get(c) == '*'` means an `N`-containing codon is neither a stop nor an error. A real-world
    consensus carries them, so this is the right call -- but it means the law cannot see a stop hidden
    behind an ambiguity code."""
    assert U.no_internal_stop_codon({"cds": "ATGNNNTCG"}).verdict == SATISFIED


# ===================================================================== the quiet standard-code default

def test_an_ABSENT_clade_quietly_defaults_to_the_STANDARD_code():
    """Recorded as a deliberate asymmetry: `table_for` REFUSES an unknown clade, but `_table` supplies
    table 1 when no clade is given at all. 'Caller said nothing' and 'caller said something wrong' are
    treated differently, which is defensible -- and worth being explicit about, because a caller who
    forgets the key gets a silent table-1 answer on a Mycoplasma sequence."""
    from dna_decode.constraints.codon_tables import STANDARD
    for ctx in ({}, {"clade": None}, {"clade": ""}):
        assert U._table(ctx) is STANDARD

    mycoplasma_cds = "ATGTGATGGTAA"                      # internal TGA: stop under 1, Trp under 4
    assert U.no_internal_stop_codon({"cds": mycoplasma_cds}).verdict == REFUTED
    assert U.no_internal_stop_codon(
        {"cds": mycoplasma_cds, "clade": "bacteria.mycoplasma"}).verdict == SATISFIED


def test_an_UNKNOWN_clade_still_RAISES_through_every_clade_aware_law():
    """The refusal must not be caught and converted to an abstention anywhere in the law layer."""
    from dna_decode.constraints.codon_tables import UnknownCladeError
    for fn, ctx in (
        (U.no_internal_stop_codon, {"cds": "ATGTCG", "clade": "martian"}),
        (U.substitution_reachable_by_single_nt, {"codon": "TCG", "mut_aa": "L", "clade": "martian"}),
        (U.start_codon_present, {"cds": "ATGTCG", "cds_is_complete_gene": True, "clade": "martian"}),
    ):
        with pytest.raises(UnknownCladeError):
            fn(ctx)


def test_evaluate_all_works_with_NO_report_and_returns_one_verdict_per_law():
    """`report=None` is the default and is the shape a library caller uses."""
    results = U.evaluate_all({"cds": "ATGTCGTGG"})
    assert list(results) == [c.key for c in U.UNIVERSAL_CONSTRAINTS]
    assert all(isinstance(v, ConstraintVerdict) for v in results.values())
    assert results["cds_length_multiple_of_three"].verdict == SATISFIED


def test_the_law_ORDER_puts_the_frame_check_before_everything_frame_dependent():
    """The module documents the ordering as part of the contract; a reshuffle would make the internal
    stop law abstain on inputs the frame law had not yet validated."""
    keys = [c.key for c in U.UNIVERSAL_CONSTRAINTS]
    assert keys.index("cds_length_multiple_of_three") < keys.index("no_internal_stop_codon")
    assert keys.index("cds_length_multiple_of_three") < keys.index("start_codon_present")


# ===================================================================== fail-open registration helpers

class _IsolatedRegistry:
    def __enter__(self):
        R.clear_registry()
        return R

    def __exit__(self, *exc):
        R.clear_registry()
        return False


def test_register_universal_does_NOT_swallow_a_provenance_refusal(monkeypatch):
    """FIXED 2026-10-01. The helper caught bare `Exception` for idempotence, which also swallowed
    `UnsourcedConstraintError` -- so the provenance gate's own helper failed OPEN and handed back an
    unregistered constraint as though it had landed. It now catches only `ConstraintConflictError`
    (the genuine already-registered case), so a provenance refusal PROPAGATES."""
    unsourced = Constraint("deliberately_unsourced", UNIVERSAL, SOURCED,
                           lambda ctx: ConstraintVerdict(SATISFIED))
    with _IsolatedRegistry():
        with pytest.raises(R.UnsourcedConstraintError):
            R.register(unsourced)                       # direct call refuses, as designed

        monkeypatch.setattr(U, "UNIVERSAL_CONSTRAINTS", (unsourced,))
        with pytest.raises(R.UnsourcedConstraintError):
            U.register_universal()                      # and the helper no longer hides it
        assert R.all_constraints() == ()


def test_register_universal_is_still_idempotent_for_the_real_laws(monkeypatch):
    """Non-vacuity for the narrowing: the idempotence the bare except was there for must still work."""
    with _IsolatedRegistry():
        first = U.register_universal()
        again = U.register_universal()
        assert len(first) == len(again) == len(U.UNIVERSAL_CONSTRAINTS)
def test_register_universal_is_genuinely_idempotent_which_is_what_the_broad_except_is_FOR():
    with _IsolatedRegistry():
        first = U.register_universal()
        second = U.register_universal()
        assert len(first) == len(second) == len(U.UNIVERSAL_CONSTRAINTS)
        assert len(R.all_constraints()) == len(U.UNIVERSAL_CONSTRAINTS)
        assert all(a is b for a, b in zip(first, second, strict=True)), (
            "a re-run returned different objects; the existing registration was not reused")


def test_register_specializations_is_idempotent_and_coexists_with_the_laws():
    with _IsolatedRegistry():
        U.register_universal()
        S.register_specializations()
        S.register_specializations()
        n_laws = len(U.UNIVERSAL_CONSTRAINTS)
        n_spec = len(S.genetic_code_specializations())
        assert len(R.all_constraints()) == n_laws + n_spec
        assert R.detect_conflicts() == []


def test_the_shipped_set_resolves_the_genetic_code_key_per_clade_without_conflict():
    """End-to-end over the real registry: every clade with a specialization resolves to its own table
    id, and a clade without one inherits nothing for that key."""
    with _IsolatedRegistry():
        U.register_universal()
        S.register_specializations()
        for clade in S._CODE_CLADES:
            got = R.constraints_for(taxon=clade)[S.GENETIC_CODE_KEY]
            assert got.payload["clade"] == clade
            assert got.kind == SOURCED and got.source_url and got.verbatim_quote
        assert S.GENETIC_CODE_KEY not in R.constraints_for(taxon="fungi")
        assert S.GENETIC_CODE_KEY not in R.constraints_for()


# ===================================================================== regime axis accessors

def test_regime_scope_REFUSES_an_unknown_key_rather_than_returning_an_empty_record():
    """An empty scoping record would read as 'this regime has no recorded verdict', which is a very
    different claim from 'there is no such regime'."""
    with pytest.raises(KeyError, match="unknown regime"):
        S.regime_scope("not_a_regime")
    with pytest.raises(KeyError):
        S.regime_scope("")


def test_regime_table_covers_EVERY_regime_and_none_are_invented():
    from dna_decode.eval.regime import REGIMES

    table = S.regime_table()
    assert set(table) == {r.key for r in REGIMES}
    assert len(table) == len(REGIMES) >= 2
    for key, row in table.items():
        assert row == S.regime_scope(key)
        assert set(row) == {"key", "population", "endpoint", "method", "verdict", "artifact",
                            "is_closed", "advisory"}
        assert row["key"] == key and row["verdict"]


def test_exactly_the_recorded_CLOSED_NEGATIVE_regime_is_flagged_closed():
    """`is_closed` is the field a caller routes on, so an over-broad flag would redirect live work."""
    closed = sorted(k for k, v in S.regime_table().items() if v["is_closed"])
    assert closed == ["natural_organism_zeroshot"], closed
    assert all(S.regime_table()[k]["verdict"] == "CLOSED_NEGATIVE" for k in closed)


def test_every_regime_row_states_that_it_cannot_refute_a_sequence():
    """The advisory string is the guard against a reader treating a regime as a constraint."""
    for row in S.regime_table().values():
        assert "never refutes" in row["advisory"]


def test_no_regime_is_registered_as_a_CONSTRAINT_on_the_regime_axis():
    """The plan said 'regime specializations as constraints'; the implementation deliberately did not,
    because an always-inapplicable constraint is exactly the inert filter `report.py` exists to catch.
    Pinned so a future change back to constraints is a visible decision."""
    with _IsolatedRegistry():
        U.register_universal()
        S.register_specializations()
        assert [c.scope for c in R.all_constraints() if c.axis == "regime"] == []


# ===================================================================== deferred clade properties

def test_each_deferred_clade_property_names_the_evidence_it_would_need():
    """The gate refusing its own author's under-sourced entries is the behaviour it exists for; the
    record of WHAT would unblock each is the part that keeps it actionable."""
    assert S.DEFERRED_CLADE_PROPERTIES, "an empty deferral list makes the refusal story unverifiable"
    for name, row in S.DEFERRED_CLADE_PROPERTIES.items():
        assert set(row) == {"buckets", "why_it_matters", "evidence_needed"}, name
        assert row["buckets"].strip(), name
        # the two narrative fields carry the actionable content; `buckets` is legitimately terse
        # ("bacteria, archaea"), so it is checked for presence only.
        assert len(row["why_it_matters"]) > 20 and len(row["evidence_needed"]) > 20, name
        assert "citable" in row["evidence_needed"], name
        assert name not in {c.key for c in S.genetic_code_specializations()}


def test_a_deferred_property_cannot_be_registered_as_AXIOMATIC_to_dodge_the_provenance_gate():
    """The gate only bites on `kind=sourced`, so the live loophole is mislabelling a clade FACT as a
    LAW. Nothing can enforce that mechanically -- this pins the distinction in the test suite so the
    dodge is at least named."""
    with _IsolatedRegistry():
        dodge = Constraint("splicing_applies", "taxon:eukaryote", AXIOMATIC,
                           lambda ctx: ConstraintVerdict(SATISFIED))
        R.register(dodge)                       # accepted: the registry cannot see the mislabel
        assert R.all_constraints()[0].kind == AXIOMATIC
        assert "splicing_applies" in S.DEFERRED_CLADE_PROPERTIES, (
            "splicing is a clade FACT and belongs on the sourced path; if it is ever registered for "
            "real it must carry source_url + verbatim_quote, not be relabelled axiomatic")
