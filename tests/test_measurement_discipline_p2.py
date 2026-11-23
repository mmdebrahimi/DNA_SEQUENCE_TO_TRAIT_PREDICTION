"""Pin the two disciplines from P2's final measurement-discipline document that this arm lacked.

Source: C:\\Users\\Laptop Portable\\pythonprojects\\dna-generator-research\\P2-tokenization-axis\\
results\\lessons-learned-measurement-discipline.md (2026-11-22, final for that axis).

Discipline 5 -- never claim "free" from a non-significant difference -- is tested against P2's OWN
withdrawn claim, which is the strongest fixture available: the function must refuse the exact
reasoning that was published and then retracted.

Discipline 1 -- the generated set is the unit of replication -- is tested by reconstructing the
sqrt(2) sem understatement numerically rather than asserting it, because the whole reason the defect
survived review is that a flat list of scoring rows looks like a flat list of replicates.

Discipline 3 -- sidedness consistency -- is pinned on this arm's existing sign-flip test. It was
already consistent; the pin exists so it cannot stop being.

No torch, no substrate, no network.
"""
from __future__ import annotations

import math
import statistics

import pytest

from dna_decode.eval.equivalence import (
    P2_REFERENCE_MARGINS,
    EquivalenceRefusal,
    equivalence,
)
from dna_decode.eval.signflip import (
    SignFlipRefusal,
    collapse_to_rows,
    family_wise_signflip,
)

# P2's withdrawn equivalence claim, verbatim from the document's table.
P2_REPATCH_POINT = 0.0098
P2_REPATCH_CI = (-0.0120, 0.0316)


# ------------------------------------------------------- discipline 5: equivalence vs a margin

def test_it_refuses_P2s_own_withdrawn_free_claim():
    """"annotation-anchored re-patching is free" rested on +0.0098 with a CI reaching +0.0316,
    against an arm spread of 0.0226. The interval contains zero AND permits a cost larger than the
    whole spread it was measured against. Those two facts together are the retraction."""
    r = equivalence(P2_REPATCH_CI, P2_REFERENCE_MARGINS["MARGIN_LOOSE"],
                    source="P2 observed arm spread (E. coli); fixture only")
    assert r.contains_zero is True, "the CI does contain zero -- that is why it read as 'no cost'"
    assert r.equivalent is False, "but it is NOT equivalent, which is what 'free' would require"
    assert r.verdict == "UNDERPOWERED_NOT_EQUIVALENT"
    assert r.worst_case_abs == pytest.approx(0.0316)
    assert r.worst_case_abs > P2_REFERENCE_MARGINS["MARGIN_LOOSE"], (
        "the permitted worst case exceeds the smallest difference anyone would act on")


def test_the_stricter_margin_refuses_it_even_harder():
    r = equivalence(P2_REPATCH_CI, P2_REFERENCE_MARGINS["MARGIN_STRICT"],
                    source="P2 leak-family per-candidate MDD (E. coli); fixture only")
    assert r.equivalent is False and r.verdict == "UNDERPOWERED_NOT_EQUIVALENT"


def test_contains_zero_and_equivalent_are_separate_fields():
    """The entire failure mode is reading the first as the second."""
    r = equivalence((-0.001, 0.001), 0.01, source="measured spread of the thing being compared")
    assert r.contains_zero and r.equivalent and r.verdict == "EQUIVALENT_WITHIN_MARGIN"
    d = r.as_dict()
    assert "contains_zero" in d and "equivalent" in d
    assert "SEPARATE" in d["reading"]


def test_an_interval_clear_of_the_margin_is_different_not_merely_unresolved():
    r = equivalence((0.05, 0.09), 0.01, source="measured spread")
    assert not r.contains_zero and not r.equivalent
    assert r.verdict == "DIFFERENT_BEYOND_MARGIN"


@pytest.mark.parametrize("margin", [0.0, -0.01])
def test_a_zero_or_negative_margin_is_refused_because_it_is_a_test_against_zero(margin):
    with pytest.raises(EquivalenceRefusal) as e:
        equivalence((-0.001, 0.001), margin, source="anything")
    assert "ZERO" in str(e.value)


@pytest.mark.parametrize("source", ["", "   ", None, 7])
def test_an_unsourced_margin_is_refused(source):
    """A margin is a claim about consequences, so an unsourced one cannot be checked later."""
    with pytest.raises(EquivalenceRefusal):
        equivalence((-0.001, 0.001), 0.01, source=source)


def test_an_inverted_or_non_finite_interval_is_refused():
    with pytest.raises(EquivalenceRefusal):
        equivalence((0.05, -0.05), 0.1, source="measured spread")
    with pytest.raises(EquivalenceRefusal):
        equivalence((float("nan"), 0.05), 0.1, source="measured spread")


def test_the_reference_margins_are_marked_DO_NOT_TRANSFER():
    """The function transfers across substrates; the numbers do not. The human bar carries an
    explicit prohibition on cross-substrate comparison, so these E. coli figures are provenance and
    a test fixture, never a human-arm margin."""
    assert "DO_NOT_TRANSFER" in P2_REFERENCE_MARGINS
    assert "E. coli" in P2_REFERENCE_MARGINS["substrate"]
    assert P2_REFERENCE_MARGINS["superseded"]["value"] == 0.0130, (
        "the superseded pooled-40 figure must stay recorded, or a later reader re-derives it")
    assert P2_REFERENCE_MARGINS["MARGIN_STRICT"] == 0.0075
    assert P2_REFERENCE_MARGINS["family_wise_mdd"] == 0.0101


# ------------------------------------------------------- discipline 1: the unit of replication

def _nested_records(n_gen=20, n_split=2, n_cand=3, effect=0.02, seed=1):
    """20 generated sets, each re-scored under n_split splits. The split seeds re-score ONE
    artifact, so within a gen_seed the value is identical -- the cleanest possible version of the
    defect, where every extra row carries exactly zero new information."""
    import random
    g = random.Random(seed)
    out = []
    for gs in range(n_gen):
        per_cand = [effect + g.gauss(0, 0.01) for _ in range(n_cand)]
        for ss in range(n_split):
            for ci, v in enumerate(per_cand):
                out.append({"gen_seed": gs, "split_seed": ss, "candidate": ci, "value": v})
    return out


def test_collapsing_recovers_the_true_replicate_count():
    recs = _nested_records()
    rows, cands, audit = collapse_to_rows(recs, replicate_key="gen_seed",
                                          candidate_key="candidate", value_key="value")
    assert len(rows) == 20, "20 generated sets are 20 replicates, not 40"
    assert cands == [0, 1, 2]
    assert audit.n_records_in == 120 and audit.n_replicates == 20
    assert audit.collapse_factor == pytest.approx(2.0)
    assert audit.nothing_collapsed is False
    assert "collapsed to 20 replicate" in audit.reads


def test_pooling_the_inner_loop_understates_the_sem_by_root_two():
    """Reconstructed numerically, not asserted. With the split seeds carrying zero extra
    information, the pooled sem is smaller by exactly sqrt(2) -- and smaller is the direction that
    manufactures significance."""
    recs = _nested_records(n_cand=1)
    rows, _, _ = collapse_to_rows(recs, replicate_key="gen_seed",
                                  candidate_key="candidate", value_key="value")
    collapsed = [r[0] for r in rows]
    pooled = [r["value"] for r in recs]

    def sem(xs):
        return statistics.stdev(xs) / math.sqrt(len(xs))

    assert len(pooled) == 2 * len(collapsed)
    ratio = sem(collapsed) / sem(pooled)
    assert ratio == pytest.approx(math.sqrt(2), rel=0.02), (
        f"expected the pooled sem to be sqrt(2) too small; ratio was {ratio}")
    assert sem(pooled) < sem(collapsed), "the defect's direction is toward significance"


def test_the_collapsed_rows_feed_the_family_wise_test_directly():
    """The two disciplines compose: collapse first, then correct for the family."""
    recs = _nested_records()
    rows, _, _ = collapse_to_rows(recs, replicate_key="gen_seed",
                                  candidate_key="candidate", value_key="value")
    r = family_wise_signflip(rows, statistic="mean")
    assert r.n_rows == 20 and r.n_candidates == 3
    assert 0.0 < r.family_wise_p <= 1.0


def test_a_replicate_missing_one_candidate_is_refused():
    """A row missing a candidate is not a whole-row sign-flip unit."""
    recs = _nested_records(n_gen=4, n_split=1, n_cand=2)
    recs = [r for r in recs if not (r["gen_seed"] == 2 and r["candidate"] == 1)]
    with pytest.raises(SignFlipRefusal) as e:
        collapse_to_rows(recs, replicate_key="gen_seed", candidate_key="candidate",
                         value_key="value")
    assert "cell(s) absent" in str(e.value)


def test_nothing_collapsed_is_reported_loudly_because_the_key_may_be_wrong():
    """A per-row id as the replicate key makes every row its own replicate -- the same pooling
    defect wearing the opposite disguise, and silent unless reported."""
    recs = _nested_records(n_gen=5, n_split=2, n_cand=1)
    for i, r in enumerate(recs):
        r["row_id"] = i
    _, _, audit = collapse_to_rows(recs, replicate_key="row_id", candidate_key="candidate",
                                  value_key="value")
    assert audit.n_replicates == 10 and audit.nothing_collapsed is True
    assert "NOTHING COLLAPSED" in audit.reads
    assert "too fine" in audit.reads


def test_a_non_finite_value_is_refused_before_it_reaches_the_mean():
    recs = _nested_records(n_gen=3, n_split=1, n_cand=1)
    recs[0]["value"] = float("nan")
    with pytest.raises(SignFlipRefusal):
        collapse_to_rows(recs, replicate_key="gen_seed", candidate_key="candidate",
                         value_key="value")


# ------------------------------------------------------- discipline 3: sidedness consistency

def test_the_family_wise_and_per_candidate_bars_are_BOTH_two_sided():
    """P2 mixed a two-sided per-candidate bar with a one-sided family-wise bar, which narrows the
    band in whichever direction flatters the claim. This arm's test was already consistent; the pin
    exists so it cannot stop being.

    Checked behaviourally: negating every value must leave both p values unchanged. A one-sided bar
    on either side would move under negation.
    """
    import random
    g = random.Random(4)
    rows = [[g.gauss(0.03, 0.02) for _ in range(4)] for _ in range(12)]
    flipped = [[-v for v in r] for r in rows]

    a = family_wise_signflip(rows, statistic="mean")
    b = family_wise_signflip(flipped, statistic="mean")
    assert a.family_wise_p == pytest.approx(b.family_wise_p), (
        "the family-wise p moved under negation, so that bar is one-sided while the per-candidate "
        "bar is two-sided -- the exact inconsistency P2 reported")
    for pa, pb in zip(a.per_candidate_p, b.per_candidate_p):
        assert pa == pytest.approx(pb)


def test_sidedness_holds_for_the_t_statistic_too():
    import random
    g = random.Random(9)
    rows = [[g.gauss(0.03, 0.02) for _ in range(3)] for _ in range(10)]
    a = family_wise_signflip(rows, statistic="max_abs_t")
    b = family_wise_signflip([[-v for v in r] for r in rows], statistic="max_abs_t")
    assert a.family_wise_p == pytest.approx(b.family_wise_p)
