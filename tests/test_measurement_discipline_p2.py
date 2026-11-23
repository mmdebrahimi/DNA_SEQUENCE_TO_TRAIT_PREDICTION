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


# ------------------------------------------------------- discipline 0: the MDD before the null

def test_the_mdd_is_reported_and_the_sems_are_the_real_ones():
    """Discipline 0, added by P2 ahead of the other ten: read what the design could detect BEFORE
    reading the p. The sems are checked against an independent computation, not against the
    implementation that produced them."""
    import random
    g = random.Random(3)
    rows = [[g.gauss(0.0, 0.05) for _ in range(4)] for _ in range(12)]
    r = family_wise_signflip(rows, statistic="max_abs_t")
    assert r.null_family_p95 > 0
    assert len(r.per_candidate_sem) == 4 and len(r.family_wise_mdd) == 4
    for j, sem in enumerate(r.per_candidate_sem):
        col = [row[j] for row in rows]
        expected = statistics.stdev(col) / math.sqrt(len(col))
        assert sem == pytest.approx(expected), f"candidate {j} sem is not sd/sqrt(n)"


def test_the_t_statistic_mdd_is_p95_times_each_candidates_OWN_sem():
    """Per discipline 4, candidates do not share a sem, so there are k minimum detectable
    differences and not one."""
    import random
    g = random.Random(8)
    # Deliberately unequal spreads so a shared-sem bug cannot pass.
    rows = [[g.gauss(0, 0.01), g.gauss(0, 0.10), g.gauss(0, 0.50)] for _ in range(14)]
    r = family_wise_signflip(rows, statistic="max_abs_t")
    for sem, mdd in zip(r.per_candidate_sem, r.family_wise_mdd):
        assert mdd == pytest.approx(r.null_family_p95 * sem)
    assert len(set(round(v, 6) for v in r.family_wise_mdd)) == 3, (
        "three candidates with different spreads must yield three different MDDs")


def test_the_mean_statistic_mdd_is_NOT_multiplied_by_a_sem():
    """The units guard. For `mean` the null p95 is ALREADY an effect size, so multiplying it by a
    sem would give a number in squared units that still reads like an effect size -- the kind of
    error that produces a plausible wrong bar rather than a crash."""
    import random
    g = random.Random(2)
    rows = [[g.gauss(0, 0.05) for _ in range(3)] for _ in range(11)]
    r = family_wise_signflip(rows, statistic="mean")
    for mdd in r.family_wise_mdd:
        assert mdd == pytest.approx(r.null_family_p95)
        assert mdd != pytest.approx(r.null_family_p95 * statistics.median(r.per_candidate_sem))


def test_a_design_limited_null_says_so_in_words():
    """P2's own headline failed this: an MDD of 0.0317 over a spread of 0.0226 means '0 of 10
    resolve' was fixed by n=20 before any data was seen."""
    import random
    g = random.Random(11)
    n = 20
    sd = 0.0106 * math.sqrt(n)
    rows = [[g.gauss(0.0, sd) for _ in range(10)] for _ in range(n)]
    r = family_wise_signflip(rows, statistic="max_abs_t", observed_spread=0.0226)
    joined = " ".join(r.notes)
    assert "DESIGN-LIMITED" in joined
    assert "statement about the DESIGN" in joined
    assert "replicate rows" in joined, "it must say how many rows would resolve the spread"
    assert statistics.median(r.family_wise_mdd) > 0.0226


def test_a_design_that_can_resolve_its_spread_is_told_apart_from_one_that_cannot():
    import random
    g = random.Random(5)
    rows = [[g.gauss(0.0, 0.004) for _ in range(3)] for _ in range(30)]
    r = family_wise_signflip(rows, statistic="max_abs_t", observed_spread=0.05)
    joined = " ".join(r.notes)
    assert "DESIGN-LIMITED" not in joined
    assert "informative about the candidates" in joined


def test_omitting_the_spread_emits_no_design_verdict():
    """Backwards compatible: the comparison needs a spread the CALLER supplies, because what counts
    as a difference worth resolving is not derivable from the data being tested."""
    import random
    g = random.Random(1)
    rows = [[g.gauss(0, 0.02) for _ in range(3)] for _ in range(10)]
    r = family_wise_signflip(rows, statistic="mean")
    joined = " ".join(r.notes)
    assert "DESIGN-LIMITED" not in joined and "informative about the candidates" not in joined
    assert r.null_family_p95 > 0, "the MDD is still computed; only the VERDICT needs a spread"


def test_the_p95_is_a_value_the_null_can_actually_produce():
    """Nearest-rank, not interpolated. Under exact enumeration the null is a finite set of
    attainable values, so an interpolated percentile would report a bar no sign vector produces."""
    rows = [[0.03, 0.04], [0.02, 0.05], [0.04, 0.03], [0.05, 0.02]]
    r = family_wise_signflip(rows, statistic="mean")
    assert r.exact
    attainable = set()
    for bits in range(2 ** len(rows)):
        signs = [1.0 if (bits >> i) & 1 else -1.0 for i in range(len(rows))]
        flipped = [[s * v for v in row] for s, row in zip(signs, rows)]
        means = [statistics.fmean([row[j] for row in flipped]) for j in range(2)]
        attainable.add(round(max(abs(m) for m in means), 9))
    assert round(r.null_family_p95, 9) in attainable


def test_the_record_states_the_power_convention_and_the_reading_order():
    import random
    g = random.Random(6)
    rows = [[g.gauss(0, 0.02) for _ in range(3)] for _ in range(10)]
    d = family_wise_signflip(rows, statistic="mean").as_dict()
    blk = d["DISCIPLINE_0_what_this_design_could_detect"]
    assert "50%" in blk["power_convention"], (
        "the MDD is a 50%-power figure and must say so; an 80% target is materially larger")
    # The key carries the emphasis; the VALUE has to carry the instruction. Asserting "BEFORE"
    # against the value was checking capitalisation rather than meaning.
    assert "read_this_BEFORE_the_p" in blk
    guidance = blk["read_this_BEFORE_the_p"]
    assert "statement about the DESIGN" in guidance
    assert "before reading family_wise_p" in guidance
    assert blk["median_family_wise_mdd"] is not None


# ------------------------------------------------- the capacity-scaling confound (P2, 2026-11-22)

def _protocol():
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "wiki" / \
        "human_replication_and_multiplicity_protocol_2026-11-22.json"
    return json.loads(p.read_text(encoding="utf-8"))["P2_FINAL_LESSONS_2026_11_22"]


def _oracle_fold9():
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "wiki" / \
        "glm_human_activity_oracle_fold9_2026-11-22.json"
    return json.loads(p.read_text(encoding="utf-8"))


def test_the_recorded_band_matches_the_oracle_artifact_it_is_derived_from():
    """The band is kmer6 minus kmer4 on held-out fold 9. Pinned against the artifact so the two
    cannot drift -- a transcribed number would go stale silently the next time the oracle is
    re-measured."""
    blk = _protocol()["CAPACITY_SCALING_CONFOUND"]
    rec = blk["THE_SHARPER_FORM_ON_THIS_ARM_THE_CHANNEL_IS_COMPOSITION_NOT_SUBSTRINGS"][
        "THE_BAND_IS_MEASURED_NOT_SPECULATED"]["per_cell"]
    o = _oracle_fold9()
    for cell in ("HepG2", "K562"):
        arms = o["cells"][cell]["geometries"]["element_200nt_adaptors_stripped"]["arms"]
        assert rec[cell]["kmer4"] == pytest.approx(arms["kmer4"])
        assert rec[cell]["kmer6"] == pytest.approx(arms["kmer6"])
        assert rec[cell]["band_above_4mers"] == pytest.approx(arms["kmer6"] - arms["kmer4"], abs=1e-4)
        assert rec[cell]["band_as_fraction_of_oracle_total"] == pytest.approx(
            (arms["kmer6"] - arms["kmer4"]) / arms["kmer6"], abs=1e-3)


def test_the_band_is_large_enough_in_K562_to_matter():
    """If this ever falls to nothing the confound stops mattering and the reporting rule can be
    relaxed. While it is a fifth of the oracle's signal, it cannot be."""
    o = _oracle_fold9()
    arms = o["cells"]["K562"]["geometries"]["element_200nt_adaptors_stripped"]["arms"]
    frac = (arms["kmer6"] - arms["kmer4"]) / arms["kmer6"]
    assert frac > 0.15, (
        f"the 5-6-mer band is {frac:.1%} of the K562 oracle's signal; the reporting rule below "
        f"exists because a generator could occupy it without exceeding composition")


def test_the_reporting_rule_says_beyond_4MER_and_never_beyond_composition():
    """THE load-bearing pin. The frozen verdict is named DESIGNS_TO_SPEC_BEYOND_4MERS and is
    accurate; the exposure is entirely in prose that rounds it up to 'beyond composition'."""
    blk = _protocol()["CAPACITY_SCALING_CONFOUND"]
    rule = blk["the_reporting_rule_adopted_here"]
    assert "4-MER" in rule
    assert "not be phrased as beyond composition" in rule
    assert "flatters the larger models" in rule, (
        "a size-curve reading must state the confound's DIRECTION, since it runs opposite to the "
        "usual one")


def test_the_protocol_records_that_copying_does_NOT_bound_this_channel():
    """The tempting and wrong inference: copying at chance does not close a composition channel.
    A model can match 6-mer frequencies while never quoting more than a few bases verbatim."""
    blk = _protocol()["CAPACITY_SCALING_CONFOUND"]
    sharp = blk["THE_SHARPER_FORM_ON_THIS_ARM_THE_CHANNEL_IS_COMPOSITION_NOT_SUBSTRINGS"]
    assert "COPYING CONTROL DOES NOT BOUND THIS" in sharp["the_consequence_that_matters"]
    assert "must NOT be read as closing this channel" in sharp["the_consequence_that_matters"]


def test_the_trigger_condition_is_recorded_as_MET_not_hypothetical():
    blk = _protocol()["CAPACITY_SCALING_CONFOUND"]
    assert "THE_TRIGGER_IS_MET_HERE" in blk
    assert "four-size curve" in blk["THE_TRIGGER_IS_MET_HERE"]["how"]
    # The emphasis lives in the KEY; the substance has to live in the VALUE. Asserting the key's
    # own wording against its value is the second instance of that mistake in this session, so
    # this pins what the sentence has to ESTABLISH instead.
    linked = blk["THE_TRIGGER_IS_MET_HERE"]["and_it_is_the_EXPOSED_family"]
    assert "RULE_2" in linked, "it must name the rule that marks this family as exposed"
    assert "sign-flip" in linked
    assert "both biases point the same way" in linked, (
        "the two biases compounding is the reason this matters, not merely coinciding")


def test_the_order_5_remedy_is_recorded_as_a_user_authority_option_not_taken():
    """Amending the comparator named by a frozen bar is not this session's call; adding a secondary
    null is additive. The distinction has to survive in the record."""
    blk = _protocol()["CAPACITY_SCALING_CONFOUND"]
    opt = blk["the_empirical_option_is_a_USER_AUTHORITY_call"]
    assert "order-5" in opt["what_would_close_it"]
    assert "ALONGSIDE" in opt["what_would_close_it"]
    assert "amend a frozen bar" in opt["why_it_is_not_taken_unilaterally"]


# ----------------------------------- the feature-family rule (P2 correction, 2026-11-22)

def test_the_deployed_scorer_is_SINGLE_FAMILY_so_the_band_is_one_dimensional():
    """The corrected rule is 'match the null to the instrument on EVERY feature family'. This arm is
    exempt on every axis but order ONLY because the deployed scorer is one mode. If that ever becomes
    an ensemble, a second band opens and this test is what notices.

    Position-invariance itself is a property of counting k-mers over the whole sequence, which the
    hand-derived fixture in tests/test_human_oracle_precondition.py already pins; this test pins that
    no OTHER family is concatenated into the scored feature vector.
    """
    import importlib.util
    import sys as _sys
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_null_under_test", repo / "scripts" / "glm_human_stratified_null.py")
    mod = importlib.util.module_from_spec(spec)
    _sys.modules["_null_under_test"] = mod
    spec.loader.exec_module(mod)

    assert mod.ORACLE_MODE == "kmer6", "the deployed scorer mode changed; re-derive the band"
    assert isinstance(mod.ORACLE_MODE, str), (
        "a list or tuple here would make the scorer multi-family and reopen the rule on every "
        "family it then spans")

    from scripts.glm_human_activity_oracle import build
    # Exactly 4**6 columns: no positional block appended to the scored vector.
    assert len(build(["ACGT" * 50], "kmer6", 10)[0]) == 4 ** 6
    assert len(build(["ACGT" * 50], "kmer4", 10)[0]) == 4 ** 4
    # And a compound mode is refused rather than silently concatenated.
    for compound in ("kmer6+positional", "both", "kmer6,positional"):
        with pytest.raises(ValueError):
            build(["ACGT" * 50], compound, 10)


def test_the_protocol_records_my_scope_error_rather_than_only_the_corrected_rule():
    """A correction that keeps only the fixed rule loses why the rule needed fixing. I checked one
    feature family and declared a peer's whole design clean."""
    cap = _protocol()["CAPACITY_SCALING_CONFOUND"]["FEATURE_FAMILY_RULE_CORRECTED"]
    assert "EXACTLY ZERO" in cap["MY_CLAIM_WAS_WRONG_AND_THE_ERROR_WAS_SCOPE"]
    assert "EVERY FEATURE FAMILY" in cap["THE_CORRECTED_RULE"]
    assert cap["WHY_THIS_ARM_IS_EXEMPT_AND_IT_IS_VERIFIED_NOT_ASSUMED"][
        "checked_because_I_had_just_made_this_exact_error_about_P2"] is True


def test_the_record_keeps_the_cost_of_the_property_that_grants_the_exemption():
    """Position-invariance spares this arm the positional band AND blinds the judge to positional
    biology. Reporting either without the other is the over-read."""
    cap = _protocol()["CAPACITY_SCALING_CONFOUND"]["FEATURE_FAMILY_RULE_CORRECTED"]
    limit = cap["THE_SAME_PROPERTY_THAT_CLOSES_THE_BAND_IS_A_LIMIT_ON_THE_SCIENCE"]
    assert "cannot reward position-dependent" in limit
    assert "neither should be reported without the other" in limit


def test_the_control_must_match_the_channel_is_recorded_as_the_generalisation():
    cap = _protocol()["CAPACITY_SCALING_CONFOUND"]["FEATURE_FAMILY_RULE_CORRECTED"]
    gen = cap["P2s_generalisation_of_my_copying_correction_which_is_better_than_mine_was"]
    assert "MATCH THE CHANNEL" in gen
    assert "SUBSTRINGS only" in gen
