"""Offline tests for the conv ⊕ max_tx rank-average hybrid arm.

The arm's own hypothesis is made falsifiable here, on synthetic data, BEFORE it meets the real table:
an orthogonal second modality must lift, and a near-chance one must not. If those two did not hold on
data we construct, a real-data result either way would be uninterpretable.
"""
from __future__ import annotations

import gzip
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import glm_tile_hybrid_arm as arm  # noqa: E402
from glm_tile_hybrid_arm import (  # noqa: E402
    PREREGISTERED,
    as_table,
    combine,
    orthogonality,
    render_md,
    shuffled,
    verdict,
)
from glm_tile_predictions import TilePredictionError, table_payload  # noqa: E402

from dna_decode.forward.variant_effect import rank_average_hybrid  # noqa: E402
from dna_decode.glm.genomewide import spearman  # noqa: E402


# --------------------------------------------------------------------------------------------------
# the combine plumbing
# --------------------------------------------------------------------------------------------------
def test_as_table_keys_are_positions_as_strings():
    assert as_table([0.5, 0.1, 0.9]) == {"0": 0.5, "1": 0.1, "2": 0.9}


def test_combine_returns_original_index_order():
    """`rank_average_hybrid` returns a dict over sorted STRING keys ('10' < '2'), so the caller must
    re-index by position or every vector past index 9 silently reorders."""
    n = 12
    preds = {"a": [float(i) for i in range(n)], "b": [float(i) for i in range(n)]}
    out = combine(preds, ("a", "b"))
    assert len(out) == n
    assert out == sorted(out), "two identical ascending scorers must combine to an ascending vector"


def test_combine_is_monotone_when_both_modalities_agree():
    preds = {"a": [1.0, 2.0, 3.0, 4.0], "b": [10.0, 20.0, 30.0, 40.0]}
    assert combine(preds, ("a", "b")) == sorted(combine(preds, ("a", "b")))


def test_rank_average_hybrid_contract_is_respected():
    """Its own refusals: fewer than 2 tables, and an empty key intersection."""
    with pytest.raises(ValueError, match="needs >=2"):
        rank_average_hybrid([{"0": 1.0}])
    with pytest.raises(ValueError, match="no variant"):
        rank_average_hybrid([{"0": 1.0}, {"1": 1.0}])


def test_combine_survives_the_two_digit_key_boundary():
    """Element-wise at n=11, the first length where a string sort disagrees with a numeric one.

    `rank_average_hybrid` keys on sorted STRINGS, so '10' sorts between '1' and '2'. At n<=10 the two
    orders coincide and the bug is invisible; n=11 is the smallest length that exposes it. The
    expected vector is written out rather than checked for monotonicity, so a reindex that merely
    happened to stay sorted could not pass.
    """
    n = 11
    v = [float(i) for i in range(n)]
    assert combine({"a": v, "b": v}, ("a", "b")) == [i / (n - 1) for i in range(n)]


def test_combine_preserves_the_input_ORDER_not_merely_sortedness():
    """The decisive form: a SHUFFLED vector combined with itself must rank-correlate 1.0 with itself.

    The existing ascending-input test passes even if positions are permuted into sorted order, because
    its correct answer IS sorted. With a shuffled input any mis-indexing shows up as rho < 1.
    """
    v = [float(i) for i in range(25)]
    random.Random(4).shuffle(v)
    out = combine({"a": v, "b": v}, ("a", "b"))
    assert spearman(out, v) == pytest.approx(1.0, abs=1e-9)


def test_combine_of_three_modalities_also_holds_the_order():
    """The `hybrid_conv_maxtx_gc` combo takes three tables, not two."""
    v = [float(i) for i in range(30)]
    random.Random(5).shuffle(v)
    out = combine({"a": v, "b": v, "c": v}, ("a", "b", "c"))
    assert len(out) == 30
    assert spearman(out, v) == pytest.approx(1.0, abs=1e-9)


def test_combine_on_a_length_mismatch_raises_one_way_and_TRUNCATES_the_other():
    """Recorded as it IS, not as it should be — and why that is safe here.

    Positional keys intersect to the shorter vector. When the FIRST name is the longer one, the
    re-index walks past the intersection and raises. When it is the shorter one, `n` is the short
    length and the result is silently truncated. The protection against the second case is upstream,
    not here: `table_payload` and `load_table` both refuse a ragged table, so a mismatch cannot reach
    `combine` from a real seed file. Pinned so the asymmetry is visible rather than assumed away.
    """
    with pytest.raises(KeyError):
        combine({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0]}, ("a", "b"))
    assert len(combine({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0]}, ("b", "a"))) == 2


def test_shuffled_preserves_the_multiset_but_not_the_order():
    v = [float(i) for i in range(50)]
    s = shuffled(v, seed=0)
    assert sorted(s) == sorted(v)
    assert s != v


def test_shuffled_is_deterministic_per_seed():
    v = [float(i) for i in range(30)]
    assert shuffled(v, 3) == shuffled(v, 3)
    assert shuffled(v, 3) != shuffled(v, 4)


# --------------------------------------------------------------------------------------------------
# the arm's HYPOTHESIS, both directions, on constructed data
# --------------------------------------------------------------------------------------------------
def _truth(n: int = 400) -> list[float]:
    return [float(i) for i in range(n)]


def test_an_ORTHOGONAL_second_modality_LIFTS():
    """The positive direction. Two noisy views of the truth with INDEPENDENT noise must combine
    better than either alone -- this is the published mechanism the arm is testing for."""
    rng = random.Random(0)
    truth = _truth()
    a = [t + rng.gauss(0, 120) for t in truth]
    b = [t + rng.gauss(0, 120) for t in truth]
    sa, sb = spearman(a, truth), spearman(b, truth)
    sh = spearman(combine({"a": a, "b": b}, ("a", "b")), truth)
    assert sh > max(sa, sb), f"hybrid {sh:.4f} did not beat max(a={sa:.4f}, b={sb:.4f})"


def test_a_NEAR_CHANCE_second_modality_REDUCES_the_signal():
    """The negative direction, and the measured failure mode (ESM2+ProSST 0.204 vs ESM2 0.352).

    Without this test a real-data negative could not be distinguished from a broken combine.
    """
    rng = random.Random(1)
    truth = _truth()
    good = [t + rng.gauss(0, 60) for t in truth]
    junk = [rng.gauss(0, 1) for _ in truth]
    sg = spearman(good, truth)
    sh = spearman(combine({"a": good, "b": junk}, ("a", "b")), truth)
    assert sh < sg, f"hybrid {sh:.4f} should be WORSE than the good modality alone {sg:.4f}"


def test_orthogonality_separates_redundant_from_independent():
    rng = random.Random(2)
    truth = _truth()
    a = [t + rng.gauss(0, 50) for t in truth]
    assert orthogonality(a, a) == pytest.approx(1.0, abs=1e-9)
    b = [rng.gauss(0, 1) for _ in truth]
    assert abs(orthogonality(a, b)) < 0.25


# --------------------------------------------------------------------------------------------------
# the frozen bar
# --------------------------------------------------------------------------------------------------
def _p(median: float, wins: int, n: int = 10) -> dict:
    return {"median": median, "wins": wins, "n": n, "ci95": [median - 0.01, median + 0.01]}


def test_verdict_hybrid_adds():
    arms = {"_paired_hybrid_vs_comparator": _p(0.03, 9),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "HYBRID_ADDS"


def test_verdict_weak_when_below_margin():
    arms = {"_paired_hybrid_vs_comparator": _p(0.005, 9),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "WEAK"


def test_verdict_weak_when_below_seed_count():
    arms = {"_paired_hybrid_vs_comparator": _p(0.03, 5),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "WEAK"


def test_verdict_no_gain_on_a_negative_margin():
    arms = {"_paired_hybrid_vs_comparator": _p(-0.004, 2),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "NO_GAIN"


def test_verdict_REFUSES_when_the_control_also_lifts():
    """A positive that the control reproduces is an artifact of rank-averaging, not evidence.

    Ordered BEFORE the HYBRID_ADDS branch on purpose: otherwise a run where both lift would print a
    headline and bury its own control.
    """
    arms = {"_paired_hybrid_vs_comparator": _p(0.03, 9),
            "_paired_control_vs_comparator": _p(0.03, 9)}
    v, detail = verdict(arms)
    assert v == "INDETERMINATE_CONTROL_ALSO_LIFTS"
    assert "artifact" in detail


def test_verdict_REFUSES_on_a_drifted_table():
    arms = {"_reconcile_problems": ["gc seed 0 = 0.2000, not 0.3 +-0.01"],
            "_paired_hybrid_vs_comparator": _p(0.05, 10),
            "_paired_control_vs_comparator": _p(0.0, 3)}
    assert verdict(arms)[0] == "INDETERMINATE_RECONCILE_FAILED"


def test_verdict_REFUSES_on_a_missing_scorer():
    assert verdict({"_missing_scorer": ["conv"]})[0] == "INDETERMINATE_MISSING_SCORER"


def test_refusals_outrank_a_headline():
    """A drifted table AND a huge apparent gain must still refuse."""
    arms = {"_missing_scorer": ["conv"], "_reconcile_problems": ["x"],
            "_paired_hybrid_vs_comparator": _p(0.99, 10),
            "_paired_control_vs_comparator": _p(0.0, 0)}
    assert verdict(arms)[0].startswith("INDETERMINATE")


def test_frozen_bar_mirrors_the_decomposition_arm():
    """Same margin and seed count, so 'does more information help?' is comparable across the two."""
    assert PREREGISTERED["margin"] == 0.02
    assert PREREGISTERED["min_seeds_winning"] == 8
    assert PREREGISTERED["comparator"] == "max_tx"


def test_the_precondition_is_reported_not_used_as_a_gate():
    """Predicting a combine will fail is not a substitute for measuring it (R6)."""
    assert PREREGISTERED["precondition"]["measured_first"] is True
    assert "not_a_gate" in PREREGISTERED["precondition"]


def test_the_stacking_bar_is_ADDITIVE_and_the_answered_question_s_bar_is_untouched():
    """The whole point of a new key: a future edit must not be able to satisfy the stacking pin by
    LOOSENING the parent. Both halves are asserted in ONE test so they cannot drift apart."""
    s = PREREGISTERED["stacking"]
    assert s["comparator"] == "hybrid_conv_maxtx", "stacking compares to the 2-channel, not the scalar"
    assert s["margin"] == 0.02
    assert s["min_seeds_winning"] == 8
    assert set(s["verdicts"]) >= {"STACKS", "WEAK_STACK", "NO_STACK"}
    # the answered question's bar, re-asserted here rather than only in the sibling test
    assert PREREGISTERED["margin"] == 0.02
    assert PREREGISTERED["min_seeds_winning"] == 8
    assert PREREGISTERED["comparator"] == "max_tx"


def test_the_stacking_prediction_and_its_negative_scope_are_recorded_BEFORE_any_run():
    """An unrecorded prediction is the thing this bar exists to prevent: it is what lets a negative
    read as a surprise and a positive read as retro-rationalised."""
    s = PREREGISTERED["stacking"]
    assert s["predicted_outcome"].strip(), "no prediction recorded"
    assert s["negative_scope"].strip(), "no negative scope recorded"
    # the prediction must name the measurement it rests on, not merely assert a direction
    assert "0.7812" in s["predicted_outcome"]
    # the scope must say what a negative does NOT close
    assert "not" in s["negative_scope"].lower()


def test_a_MISSING_stacking_measurement_is_not_the_same_verdict_as_a_negative_one():
    """The distinction Step 2's guards exist to preserve: an absent statistic and a measured negative
    must not collapse to one verdict."""
    v = PREREGISTERED["stacking"]["verdicts"]
    assert "NO_STACK" in v and "INDETERMINATE_STACKING_NOT_MEASURED" in v
    assert v["NO_STACK"] != v["INDETERMINATE_STACKING_NOT_MEASURED"]


# --------------------------------------------------------------------------------------------------
# the COMBOS-derived surface -- generality is asserted against COMBOS itself, never against a literal
# --------------------------------------------------------------------------------------------------
def test_COMBOS_has_four_entries_and_the_two_ORIGINAL_ones_are_UNCHANGED():
    """Asserted LITERALLY on purpose: the published 0.4254 and 0.4217 are only reproducible while
    these two keys and tuples are exactly what produced them."""
    assert len(arm.COMBOS) == 4
    assert arm.COMBOS["hybrid_conv_maxtx"] == ("conv", "max_tx")
    assert arm.COMBOS["hybrid_conv_maxtx_gc"] == ("conv", "max_tx", "gc")
    assert "hybrid_conv_maxtx_decomp" in arm.COMBOS
    assert "hybrid_conv_maxtx_gc_decomp" in arm.COMBOS


def test_the_decomposition_combos_are_supersets_so_Step_2_s_rule_GIVES_them_a_stacking_statistic():
    """Checked rather than assumed: a new combo that was not a strict superset would silently get no
    `_paired_*_vs_hybrid2` and so be judged by nothing."""
    sup = set(arm.superset_combos())
    assert {"hybrid_conv_maxtx_decomp", "hybrid_conv_maxtx_gc_decomp"} <= sup
    assert "decomposition" in arm.all_channels()
    # and the matrix row set widened with them
    assert ("max_tx", "decomposition") in arm.channel_pairs() or            ("decomposition", "max_tx") in arm.channel_pairs()


def test_the_stack_base_and_required_channels_are_DERIVED_from_COMBOS():
    """If these were restated literals, a COMBOS edit would silently desynchronise them."""
    assert arm.stack_base() == arm.COMBOS[PREREGISTERED["stacking"]["comparator"]]
    assert arm.required_channels() == arm.stack_base()
    # the union is strictly wider than the required set as soon as any combo adds a channel
    assert set(arm.all_channels()) >= set(arm.required_channels())


def test_every_channel_pair_any_combo_uses_appears_exactly_once_unordered():
    ch = arm.all_channels()
    pairs = arm.channel_pairs()
    assert len(pairs) == len(ch) * (len(ch) - 1) // 2
    assert len({frozenset(pr) for pr in pairs}) == len(pairs), "a pair is listed twice"
    for a_, b_ in pairs:
        assert a_ != b_


def test_superset_combos_are_STRICT_supersets_of_the_stacking_comparator():
    base = set(arm.stack_base())
    sup = arm.superset_combos()
    assert PREREGISTERED["stacking"]["comparator"] not in sup, "the comparator cannot stack on itself"
    for k in sup:
        assert set(arm.COMBOS[k]) > base, f"{k} is not a strict superset of {sorted(base)}"


def test_the_two_published_paired_key_names_are_preserved_verbatim():
    """Published artifacts and pre-existing tests read these two names. Renaming them for symmetry
    would break backward compatibility for a cosmetic gain."""
    assert arm.paired_key("hybrid_conv_maxtx") == "_paired_hybrid_vs_comparator"
    assert arm.paired_key("hybrid_conv_maxtx_gc") == "_paired_hybrid3_vs_comparator"
    # and a combo with no legacy name gets a generated one that cannot collide with them
    gen = arm.paired_key("some_future_combo")
    assert gen == "_paired_some_future_combo_vs_comparator"
    assert gen not in set(arm.LEGACY_PAIRED_KEY.values())


def test_the_stacking_pair_key_is_a_DIFFERENT_namespace_from_the_comparator_key():
    """`vs_hybrid2` and `vs_comparator` answer different questions; one key for both would make the
    stacking number silently overwrite the published one."""
    for k in arm.COMBOS:
        assert arm.stacking_pair_key(k) != arm.paired_key(k)
    assert set(arm.LEGACY_PAIRED_KEY.values()).isdisjoint(
        {arm.stacking_pair_key(k) for k in arm.COMBOS})


def test_skipped_combos_names_the_MISSING_channels_and_is_empty_on_a_complete_table():
    full = {"predictions": {c: [1.0, 2.0] for c in arm.all_channels()}}
    assert arm.skipped_combos([full]) == {}
    # drop ONE optional channel -- derived from COMBOS, not named literally
    optional = [c for c in arm.all_channels() if c not in arm.required_channels()]
    if optional:
        drop = optional[0]
        partial = {"predictions": {c: [1.0, 2.0] for c in arm.all_channels() if c != drop}}
        sk = arm.skipped_combos([partial])
        assert sk, f"dropping {drop} skipped nothing"
        for k, gaps in sk.items():
            assert drop in gaps
            assert drop in arm.COMBOS[k]
        # a combo that does NOT name the dropped channel must not be skipped
        for k, names in arm.COMBOS.items():
            if drop not in names:
                assert k not in sk


def test_the_SHORT_ARM_guard_refuses_rather_than_reporting_a_short_n():
    """`if vals` previously accepted 3 of 10 and recorded `"n": 3` with no flag."""
    arms = {"_short_arm": {"hybrid_conv_maxtx": {"scored": 3, "expected": 10}},
            "_paired_hybrid_vs_comparator": _p(0.99, 10),
            "_paired_control_vs_comparator": _p(0.0, 0)}
    v, detail = verdict(arms)
    assert v == "INDETERMINATE_SHORT_ARM"
    assert "3 of 10" in detail
    assert v in PREREGISTERED["verdicts"], "a verdict the frozen map does not describe"


def test_the_short_arm_refusal_OUTRANKS_a_headline_sized_gain():
    """Ordering, asserted separately from existence: a 0.99 gain must not win over the refusal."""
    arms = {"_short_arm": {"x": {"scored": 1, "expected": 10}},
            "_paired_hybrid_vs_comparator": _p(0.99, 10),
            "_paired_control_vs_comparator": _p(0.0, 0)}
    assert verdict(arms)[0].startswith("INDETERMINATE")


# --------------------------------------------------------------------------------------------------
# stacking_verdict -- the second question, kept out of the top-level field
# --------------------------------------------------------------------------------------------------
def test_stacking_verdict_REFUSES_when_no_superset_statistic_exists():
    """A missing measurement must not read as NO_STACK -- the seam where the two would collapse."""
    v, detail, per = arm.stacking_verdict({"_paired_hybrid_vs_comparator": _p(0.03, 9)})
    assert v == "INDETERMINATE_STACKING_NOT_MEASURED"
    assert per == {}
    assert "not measured" in detail


def test_stacking_verdict_applies_the_STACKING_bar_not_the_parent_bar():
    s = PREREGISTERED["stacking"]
    assert arm.stacking_verdict_for(_p(s["margin"], s["min_seeds_winning"]))[0] == "STACKS"
    assert arm.stacking_verdict_for(_p(s["margin"] / 2, s["min_seeds_winning"]))[0] == "WEAK_STACK"
    # CHANGED 2026-10-26, and this one assertion IS the change. It previously expected WEAK_STACK for
    # a gain that clears the effect-size bar on only `min_seeds_winning - 1` seeds -- which gave that
    # case the same verdict as a gain reproducing on every seed. A result can be large and
    # unreproducible; the verdict now says so instead of flattening the two together.
    assert (arm.stacking_verdict_for(_p(s["margin"], s["min_seeds_winning"] - 1))[0]
            == "INDETERMINATE_GAIN_NOT_REPRODUCIBLE")
    assert arm.stacking_verdict_for(_p(-0.01, 2))[0] == "NO_STACK"


def test_the_verdict_SEPARATES_a_reproducible_gain_from_positive_NOISE():
    """THE defect this change closes, pinned as the first thing a reader meets.

    Before: `WEAK_STACK` fired on `median > 0` alone, so positive noise on 4 of 10 seeds with a CI
    spanning zero and a gain reproducing on 10 of 10 with a CI excluding zero returned the SAME
    string. That is the identical collapse `INDETERMINATE_STACKING_NOT_MEASURED` exists to prevent
    one key away -- two materially different findings sharing one machine-readable verdict.
    """
    noise = {"median": 0.0001, "wins": 4, "n": 10, "ci95": [-0.0200, 0.0202]}
    real = {"median": 0.0145, "wins": 10, "n": 10, "ci95": [0.0108, 0.0172]}
    v_noise = arm.stacking_verdict_for(noise)[0]
    v_real = arm.stacking_verdict_for(real)[0]
    assert v_noise != v_real, "a reproducible gain and positive noise still share one verdict"
    assert v_real == "WEAK_STACK"
    assert v_noise == "INDETERMINATE_GAIN_NOT_REPRODUCIBLE"
    # and the refusal must SAY it is neither, so nobody reads it as a measured negative
    detail = arm.stacking_verdict_for(noise)[1]
    assert "did NOT reproduce" in detail
    assert "neither a measured stack nor a measured negative" in detail


def test_the_NARROWED_bar_moved_no_measured_verdict():
    """The neutrality claim the narrowed bar rests on, checked rather than argued.

    `PREREGISTERED["stacking"]["verdicts"]` names this test by name: editing a pre-registered bar
    after the result is known is only defensible if it is OUTCOME-NEUTRAL on the data already in
    hand. The three statistics below are the MEASURED ones; under the narrowed bar each must keep the
    verdict it was published with. If one moved, the correct response was to revert the narrowing.
    """
    measured = {
        "hybrid_conv_maxtx_decomp": ({"median": 0.0145, "wins": 10, "n": 10,
                                      "ci95": [0.0108, 0.0172]}, "WEAK_STACK"),
        "hybrid_conv_maxtx_gc_decomp": ({"median": 0.0103, "wins": 9, "n": 10,
                                         "ci95": [0.0039, 0.0156]}, "WEAK_STACK"),
        "hybrid_conv_maxtx_gc": ({"median": -0.0039, "wins": 4, "n": 10,
                                  "ci95": [-0.0095, 0.0026]}, "NO_STACK"),
    }
    for combo, (stat, want) in measured.items():
        got = arm.stacking_verdict_for(stat)[0]
        assert got == want, f"{combo} moved {want} -> {got}; the narrowing was NOT outcome-neutral"

    # the aggregate over exactly those three is unchanged too
    arms = {arm.stacking_pair_key(c): s for c, (s, _w) in measured.items()}
    assert arm.stacking_verdict(arms)[0] == "WEAK_STACK"


def test_the_narrowed_condition_READS_the_frozen_seed_count_and_is_not_a_new_literal():
    """No new threshold was introduced -- the narrowing reuses `min_seeds_winning`, already frozen at
    8 and already load-bearing for `STACKS`. Proven BEHAVIOURALLY rather than by parsing the source
    for digits (the f-strings carry digits of their own): feed a bar with a DIFFERENT seed count and
    the boundary must move with it.
    """
    s = PREREGISTERED["stacking"]
    stat = {"median": s["margin"] / 2, "wins": 6, "n": 10, "ci95": [0.0, 0.02]}
    # under the real bar (8 seeds) 6 wins is not reproducible
    assert arm.stacking_verdict_for(stat)[0] == "INDETERMINATE_GAIN_NOT_REPRODUCIBLE"
    # under a bar asking for 5, the same statistic IS reproducible -- so the condition reads the key
    relaxed = {"stacking": dict(s, min_seeds_winning=5)}
    assert arm.stacking_verdict_for(stat, relaxed)[0] == "WEAK_STACK"
    # and the frozen values themselves are untouched
    assert s["margin"] == 0.02
    assert s["min_seeds_winning"] == 8


def test_the_three_OTHER_verdict_strings_are_byte_unchanged():
    """Only `WEAK_STACK` was narrowed. Pinning the other three verbatim stops a future edit quietly
    widening a second branch under cover of this one."""
    v = PREREGISTERED["stacking"]["verdicts"]
    assert v["STACKS"] == ("paired median >= +0.02 over hybrid_conv_maxtx AND >= 8/10 seeds win")
    assert v["NO_STACK"] == ("no positive paired gain -- the third channel carries nothing the "
                            "measured 2-channel hybrid lacks")
    assert v["INDETERMINATE_STACKING_NOT_MEASURED"] == (
        "no superset combo produced a paired statistic -- a MISSING measurement, which must not be "
        "read as a negative one")
    # the narrowed one must now carry BOTH conditions, not just the gain
    assert "min_seeds_winning" in v["WEAK_STACK"]
    assert "REPRODUCIBLE" in v["WEAK_STACK"]


def test_an_ALL_unreproducible_positive_arm_does_NOT_aggregate_to_NO_STACK():
    """The precedence bug the new key would otherwise create one level up.

    `stacking_verdict` returns the first verdict present in a best-first tuple. If the new key is
    missing from that tuple, an arm whose every combo is positive-but-unreproducible falls through to
    the `NO_STACK` return -- reporting a measured NEGATIVE for a set of positive gains, which is the
    original defect re-created at the aggregate.
    """
    sup = arm.superset_combos()
    assert sup
    arms = {arm.stacking_pair_key(c): {"median": 0.001, "wins": 3, "n": 10, "ci95": [-0.01, 0.02]}
            for c in sup}
    v, detail, per = arm.stacking_verdict(arms)
    assert v == "INDETERMINATE_GAIN_NOT_REPRODUCIBLE"
    assert v != "NO_STACK"
    assert all(r["verdict"] == "INDETERMINATE_GAIN_NOT_REPRODUCIBLE" for r in per.values())
    assert detail.strip()


def test_the_aggregate_detail_names_EVERY_combo_including_the_ones_that_LOST():
    """The aggregate is best-of by design; the detail string must not present best-of as all-of.

    Measured on the shipped artifact before the fix: `stacking_per_combo` recorded a `NO_STACK` combo
    that `stacking_verdict_detail` never mentioned, while `render_md`'s per-combo table did list it --
    so the human surface was honest and the JSON field a programmatic reader consumes was not. This
    assertion fails against the pre-change function, which joined over the winners only.
    """
    sup = arm.superset_combos()
    assert len(sup) >= 2, "need at least one winner and one loser to make this non-vacuous"
    winner, *losers = sup
    arms = {arm.stacking_pair_key(c): {"median": -0.01, "wins": 2, "n": 10, "ci95": [-0.02, 0.0]}
            for c in losers}
    arms[arm.stacking_pair_key(winner)] = {"median": 0.01, "wins": 10, "n": 10,
                                           "ci95": [0.005, 0.015]}
    v, detail, per = arm.stacking_verdict(arms)
    assert v == "WEAK_STACK", "the aggregate verdict itself must stay best-of"
    for c in sup:
        assert c in detail, f"{c} is in stacking_per_combo but absent from the detail string"
        assert f"{c} [{per[c]['verdict']}]" in detail, f"{c} is named without its own verdict"
    # the loser's verdict must be visible in the string, not only in the per-combo block
    assert "NO_STACK" in detail
    assert any(per[c]["verdict"] != v for c in sup), "fixture does not actually mix verdicts"


def test_stacking_verdict_aggregates_over_every_superset_combo_that_RAN():
    sup = arm.superset_combos()
    assert sup, "no superset combo to aggregate over"
    arms = {arm.stacking_pair_key(c): _p(-0.01, 2) for c in sup}
    v, _d, per = arm.stacking_verdict(arms)
    assert v == "NO_STACK"
    assert set(per) == set(sup)
    # one winner is enough for the aggregate, and the per-combo block still shows the losers
    arms[arm.stacking_pair_key(sup[0])] = _p(0.05, 10)
    v2, _d2, per2 = arm.stacking_verdict(arms)
    assert v2 == "STACKS"
    assert per2[sup[0]]["verdict"] == "STACKS"
    assert all(v2k in PREREGISTERED["stacking"]["verdicts"]
               for v2k in {r["verdict"] for r in per2.values()})


def test_the_top_level_verdict_does_NOT_change_meaning_when_stacking_is_present():
    """Risk flag 5: the published `verdict` field stays the 2-channel question's answer."""
    arms = {"_paired_hybrid_vs_comparator": _p(0.04, 9),
            "_paired_control_vs_comparator": _p(0.0, 1),
            arm.stacking_pair_key(arm.superset_combos()[0]): _p(-0.5, 0)}
    assert verdict(arms)[0] == "HYBRID_ADDS", "a stacking negative changed the answered question"


# --------------------------------------------------------------------------------------------------
# render_md -- the human-facing surface. A caveat carried only in the JSON is not a disclosure.
# --------------------------------------------------------------------------------------------------
def _art(**over) -> dict:
    art = {
        "date": "2026-10-14",
        "verdict": "WEAK",
        "verdict_detail": "hybrid +0.0054 over max_tx, 7/10 seeds",
        "orthogonality": {"median": 0.6421, "min": 0.5900, "max": 0.7100,
                          "per_seed": [0.59, 0.64, 0.71]},
        "arms": {
            "max_tx": {"median": 0.3968, "n": 2, "per_seed": [0.39, 0.40]},
            "conv": {"median": 0.3234, "n": 2, "per_seed": [0.32, 0.33]},
            "hybrid_conv_maxtx": {"median": 0.4022, "n": 2, "per_seed": [0.40, 0.41]},
            "_paired_hybrid_vs_comparator": _p(0.0054, 7),
            "_paired_control_vs_comparator": _p(0.0001, 4),
        },
        "honest_limits": ["it RANKS and does not dose", "one rung, one combine rule"],
    }
    art.update(over)
    return art


def test_render_md_carries_the_verdict_and_its_detail():
    md = render_md(_art())
    assert "**VERDICT: `WEAK`**" in md
    assert "hybrid +0.0054 over max_tx" in md


def test_render_md_reports_the_orthogonality_precondition_with_its_SPREAD():
    """The median alone would hide a seed where the two scorers were redundant."""
    md = render_md(_art())
    assert "0.6421" in md and "0.5900" in md and "0.7100" in md


def test_render_md_renders_one_row_per_PUBLIC_arm_and_never_a_private_key():
    """Non-vacuous: the `_paired_*` dicts carry no `per_seed`, so a leaked private key would raise
    rather than merely render an odd row. That this call SUCCEEDS is the assertion."""
    md = render_md(_art())
    assert "| `max_tx` |" in md and "| `conv` |" in md and "| `hybrid_conv_maxtx` |" in md
    assert "_paired_hybrid_vs_comparator" not in md.split("## Paired margins")[0]


def test_render_md_omits_a_paired_line_whose_arm_did_not_run():
    """`hybrid+gc` exists only when `gc` is in the table; its absence must not print an empty row."""
    md = render_md(_art())
    assert "hybrid vs max_tx" in md
    assert "CONTROL (shuffled conv) vs max_tx" in md
    assert "hybrid+gc vs max_tx" not in md


def test_render_md_includes_the_three_way_line_when_that_arm_ran():
    a = _art()
    a["arms"]["_paired_hybrid3_vs_comparator"] = _p(0.0070, 8)
    assert "hybrid+gc vs max_tx" in render_md(a)


def test_render_md_always_renders_the_control_line():
    """A positive the control reproduces is an artifact; the control must be visible beside it."""
    assert "CONTROL (shuffled conv)" in render_md(_art(verdict="HYBRID_ADDS"))


def test_render_md_carries_every_honest_limit():
    md = render_md(_art())
    for lim in _art()["honest_limits"]:
        assert lim in md


# --------------------------------------------------------------------------------------------------
# main() -- the refusal and write paths, on synthetic seed tables. No D:, no model, no network.
# --------------------------------------------------------------------------------------------------
def _noisy(n: int, sd: float, seed: int) -> list[float]:
    r = random.Random(seed)
    return [float(i) + r.gauss(0, sd) for i in range(n)]


def _write_seed(d: Path, seed: int, preds: dict[str, list[float]],
                sp: dict[str, float], *, smoke: bool = False, n: int = 40) -> None:
    """One synthetic seed table on disk, built through the script's own validating payload builder."""
    d.mkdir(parents=True, exist_ok=True)
    payload = table_payload(seed=seed, held=["peakA"], tile_indices=list(range(n)),
                            measured=[float(i) for i in range(n)], predictions=preds,
                            spearmans=sp, smoke=smoke)
    with gzip.open(d / f"seed-{seed:02d}.json.gz", "wt", encoding="utf-8") as fh:
        json.dump(payload, fh)


def _full_seed(d: Path, seed: int, *, smoke: bool = False) -> None:
    """A table carrying all three scorers, with seed 0's recomputed Spearmans ON the published
    anchors so the reconcile gate is clean."""
    _write_seed(d, seed,
                {"gc": _noisy(40, 8, 10 + seed), "max_tx": _noisy(40, 6, 20 + seed),
                 "conv": _noisy(40, 7, 30 + seed)},
                {"gc": 0.3000 if seed == 0 else 0.3100,
                 "max_tx": 0.3968 if seed == 0 else 0.4000, "conv": 0.3234},
                smoke=smoke)


def test_main_REFUSES_an_empty_input_directory(tmp_path, monkeypatch):
    """Nothing to combine is not a NO_GAIN. It must name the producer to run first."""
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path)])
    with pytest.raises(SystemExit, match="no seed tables"):
        arm.main()


def test_main_REFUSES_a_table_without_the_scorers_it_combines(tmp_path, monkeypatch, capsys):
    """A gc-only table cannot answer this arm's question; main must exit non-zero and write nothing."""
    _write_seed(tmp_path, 0, {"gc": _noisy(40, 5, 1)}, {"gc": 0.3000})
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1"])
    assert arm.main() == 1
    out = capsys.readouterr().out
    assert "INDETERMINATE_MISSING_SCORER" in out
    assert not list(tmp_path.glob("*.json")), "a refusal must not leave an artifact behind"


def test_main_REFUSES_a_SMOKE_table_at_the_consumer_boundary(tmp_path, monkeypatch):
    """The laundering guard where it matters: a smoke table skipped the reconcile gate, so the entry
    point itself must refuse rather than relying on a caller to have checked."""
    _full_seed(tmp_path, 0, smoke=True)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1"])
    with pytest.raises(TilePredictionError, match="SMOKE"):
        arm.main()


def test_main_writes_both_the_json_and_the_markdown_sidecar(tmp_path, monkeypatch):
    """End to end on two synthetic seeds. The MD sidecar is the only surface a human reads, so a run
    that wrote JSON alone would be a silent half-report."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    assert arm.main() == 0
    assert out.is_file() and out.with_suffix(".md").is_file()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["schema"] == "glm-tile-hybrid-arm-v1"
    assert art["n_seeds"] == 2
    assert art["preregistered"]["margin"] == PREREGISTERED["margin"]
    assert art["arm_name"] == "promcalc-plus-conv-rank-average-hybrid-150bp"
    assert art["verdict"] in PREREGISTERED["verdicts"]
    assert "**VERDICT:" in out.with_suffix(".md").read_text(encoding="utf-8")


def test_main_cannot_reach_HYBRID_ADDS_on_fewer_seeds_than_the_bar_needs(tmp_path, monkeypatch):
    """The bar wants 8 winning seeds of 10. A 2-seed run must therefore not print the headline, even
    when both seeds win — the same under-scoped-run trap the selection script carries a branch for."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    arm.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["verdict"] != "HYBRID_ADDS"
    assert art["arms"]["_paired_hybrid_vs_comparator"]["n"] == 2


def test_main_SKIPS_an_optional_combo_and_RECORDS_why_rather_than_refusing_the_whole_arm(
        tmp_path, monkeypatch):
    """The Step-2 deviation, pinned.

    `_full_seed` writes gc/max_tx/conv and NO `decomposition`. Refusing the whole arm on that -- the
    plan's literal instruction, a union over COMBOS -- would make the ANSWERED question (hybrid vs
    max_tx, published 0.0339) unreproducible from any table predating the decomposition combos. So
    the optional combos are skipped UNIFORMLY and recorded, and the primary question still answers.

    What makes this safe rather than lax is the pair of guards: the skip is all-or-nothing per combo,
    and `_short_arm` refuses if any arm ever scores a seed subset anyway.
    """
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    assert arm.main() == 0, "a missing OPTIONAL channel must not refuse the whole arm"
    art = json.loads(out.read_text(encoding="utf-8"))
    a = art["arms"]

    # the primary question still answered, and never on a short arm
    assert "_paired_hybrid_vs_comparator" in a
    assert a["hybrid_conv_maxtx"]["n"] == 2
    assert "_short_arm" not in a

    # the combos needing the absent channel are skipped, recorded, and name it -- derived from COMBOS
    absent = "decomposition"
    expected = {k for k, names in arm.COMBOS.items() if absent in names}
    assert set(a["_skipped_combos"]) == expected
    for k in expected:
        assert absent in a["_skipped_combos"][k]
        assert k not in a, "a skipped combo must not also report a median"
        assert arm.paired_key(k) not in a
        assert arm.stacking_pair_key(k) not in a

    # and a skipped combo is a MISSING measurement, never a negative one
    assert art["stacking_verdict"] in PREREGISTERED["stacking"]["verdicts"]
    assert set(art["stacking_per_combo"]).isdisjoint(expected)
    # the skip reaches the human-facing surface, not only the JSON
    md = out.with_suffix(".md").read_text(encoding="utf-8")
    assert "SKIPPED" in md and absent in md


def test_main_REFUSES_when_a_REQUIRED_channel_is_absent_even_though_optional_ones_are_tolerated(
        tmp_path, monkeypatch, capsys):
    """The other half of the same trade: tolerance is scoped to OPTIONAL channels only. Derived from
    `required_channels()`, so it cannot pass by hardcoding `conv`."""
    req = arm.required_channels()[0]
    preds = {c: _noisy(40, 5, 7) for c in arm.all_channels() if c != req}
    _write_seed(tmp_path, 0, preds, {c: 0.3 for c in preds})
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1"])
    assert arm.main() == 1
    out = capsys.readouterr().out
    assert "INDETERMINATE_MISSING_SCORER" in out
    assert req in out
    assert not list(tmp_path.glob("*.json")), "a refusal must not leave an artifact behind"


def test_main_REFUSES_a_SHORT_ARM_end_to_end_under_MUTATION(tmp_path, monkeypatch, capsys):
    """Risk flag 4: on real data every table carries every channel, so this guard CANNOT fire and is
    `[inferred]` rather than `[grounded]` in production. Mutation is therefore the only way to
    exercise it: `skipped_combos` is forced to UNDER-REPORT, which is precisely the bug the guard
    exists to catch. One table carries `decomposition`, the other does not, so the decomposition
    combos score 1 of 2 seeds -- the old `if vals` would have recorded `"n": 1` beside 2-seed
    medians with no flag at all.
    """
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    # seed 0's recomputed Spearmans must sit ON the published anchors or the RECONCILE gate fires
    # first and this test measures the wrong refusal -- which is exactly what it did on first run.
    _write_seed(tables, 0,
                {c: _noisy(40, 6, 10 + i) for i, c in enumerate(arm.all_channels())},
                {"gc": 0.3000, "max_tx": 0.3968, "conv": 0.3234, "decomposition": 0.4021})
    _full_seed(tables, 1)  # no `decomposition`
    monkeypatch.setattr(arm, "skipped_combos", lambda _t: {})
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    assert arm.main() == 1
    printed = capsys.readouterr().out
    assert "INDETERMINATE_SHORT_ARM" in printed
    assert "1 of 2 seeds" in printed
    # the refusal names WHICH arms were short, and only those
    short_combos = {k for k, names in arm.COMBOS.items() if "decomposition" in names}
    for k in short_combos:
        assert k in printed
    assert not out.is_file(), "a refusal must not leave an artifact behind"


def test_main_writes_an_orthogonality_matrix_for_every_PRESENT_channel_pair(tmp_path, monkeypatch):
    """Step 4. The row set is derived from COMBOS, so a new channel cannot be silently unmeasured."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    arm.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    mx = art["orthogonality_matrix"]

    present = {c for c in arm.all_channels() if c != "decomposition"}
    want = {f"{a_}__{b_}" for a_, b_ in arm.channel_pairs()
            if a_ in present and b_ in present}
    assert set(mx) == want, "matrix rows do not match the channel pairs actually present"

    # the matrix must AGREE with the single number it generalises -- not merely coexist with it
    assert mx["conv__max_tx"]["median"] == art["orthogonality"]["median"]
    assert mx["conv__max_tx"]["per_seed"] == art["orthogonality"]["per_seed"]

    # ordered most-redundant first, so the reader meets the obstacle before the pair that worked
    meds = [m["median"] for m in mx.values()]
    assert meds == sorted(meds, reverse=True)

    for m in mx.values():
        assert len(m["per_seed"]) == 2
        assert m["min"] <= m["median"] <= m["max"]

    assert "Pairwise orthogonality" in out.with_suffix(".md").read_text(encoding="utf-8")


def test_the_published_orthogonality_KEY_keeps_its_exact_shape(tmp_path, monkeypatch):
    """The matrix is ADDITIVE. Published artifacts and `test_main_records_the_orthogonality_...`
    read `orthogonality`; replacing it with the matrix would break both."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    arm.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert set(art["orthogonality"]) == {"per_seed", "median", "min", "max"}
    assert "orthogonality_matrix" in art and art["orthogonality_matrix"] is not art["orthogonality"]


def test_render_md_tolerates_an_artifact_with_NEITHER_new_section(tmp_path):
    """Backward compatibility of the renderer itself: the published 2026-10-14 artifact carries no
    matrix and no stacking verdict, and must still render."""
    md = render_md(_art())
    assert "**VERDICT:" in md
    assert "Pairwise orthogonality" not in md
    assert "STACKING VERDICT" not in md


def test_the_SHIPPED_artifact_carries_both_new_honest_limits_and_both_REACH_the_markdown():
    """The two limits added 2026-10-26 are the scoping of a PUBLISHED claim and the record of a
    finding deliberately not fixed. A limit that lives only in the JSON is not a disclosure, so both
    must reach `render_md`'s rendered section as well.

    Reads the real artifact; skips when absent, because it is regenerated rather than committed-only.
    """
    art_path = Path(__file__).resolve().parents[1] / "wiki" / "glm_tile_hybrid_arm_2026-10-26.json"
    if not art_path.exists():
        pytest.skip("the shipped artifact is not on disk; re-run scripts/glm_tile_hybrid_arm.py")
    art = json.loads(art_path.read_text(encoding="utf-8"))
    limits = art["honest_limits"]
    assert len(limits) == 8, f"expected 8 honest limits, found {len(limits)}"

    scope, ci = limits[6], limits[7]
    # item 7 SCOPES the direction finding: the stacking pair is a reduction/superset pair, not two
    # modalities, so prediction-vector Spearman is refuted across THAT pair and not in general.
    for tok in ("max_tx", "decomposition", "REDUCTION", "modality"):
        assert tok in scope, f"the scoping limit does not name {tok!r}"
    # item 8 RECORDS the CI estimand finding and says plainly that it is not fixed and why
    for tok in ("MEAN", "MEDIAN", "t(9)", "promcalc_decomposition_arm"):
        assert tok in ci, f"the CI limit does not name {tok!r}"

    md = render_md(art)
    assert "## Honest limits" in md
    assert "REDUCTION" in md, "the scoping limit never reaches a human"
    assert "t(9)" in md, "the CI limit never reaches a human"


def test_main_records_the_orthogonality_and_the_control_arm(tmp_path, monkeypatch):
    """Both are preconditions for reading the headline, so neither may be absent from the artifact."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    arm.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert set(art["orthogonality"]) >= {"per_seed", "median", "min", "max"}
    assert len(art["orthogonality"]["per_seed"]) == 2
    assert "control_shuffled_conv" in art["arms"]
    assert "_paired_control_vs_comparator" in art["arms"]


# --------------------------------------------------------------------------------------------------
# paired_label -- the human-facing name the two NEW combos render under
# --------------------------------------------------------------------------------------------------
def test_paired_label_names_the_EXTRA_channels_and_never_collides_with_a_legacy_label():
    """The two legacy labels are published wording; the generated ones must stay distinguishable from
    them AND from each other, or two paired lines in the markdown read as the same arm."""
    assert arm.paired_label("hybrid_conv_maxtx") == "hybrid"
    assert arm.paired_label("hybrid_conv_maxtx_gc") == "hybrid+gc"
    # the generated branch names exactly the channels the combo adds over the 2-channel base
    for combo in arm.superset_combos():
        extra = [c for c in arm.COMBOS[combo] if c not in set(arm.stack_base())]
        lab = arm.paired_label(combo)
        for c in extra:
            assert c in lab, f"{combo}'s label {lab!r} does not name its added channel {c}"
        for c in arm.stack_base():
            assert c not in lab, f"{combo}'s label {lab!r} restates a base channel"
    labels = [arm.paired_label(c) for c in arm.COMBOS]
    assert len(set(labels)) == len(labels), f"two combos render under one label: {labels}"


# --------------------------------------------------------------------------------------------------
# the COMPLETE-table path -- every other main() test writes a table WITHOUT `decomposition`, so the
# stacking statistic is skipped and the stacking verdict refuses. This is the only place the second
# question is actually MEASURED end to end.
# --------------------------------------------------------------------------------------------------
def _quad_seed(d: Path, seed: int) -> None:
    """A table carrying ALL FOUR channels, seed 0's recomputed Spearmans on the published anchors."""
    _write_seed(d, seed,
                {c: _noisy(40, 6 + i, 10 * i + seed) for i, c in enumerate(arm.all_channels())},
                {"gc": 0.3000 if seed == 0 else 0.3100,
                 "max_tx": 0.3968 if seed == 0 else 0.4000,
                 "conv": 0.3234, "decomposition": 0.4021})


def _complete_run(tmp_path, monkeypatch, n_seeds: int = 2) -> tuple[dict, str]:
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in range(n_seeds):
        _quad_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", str(n_seeds),
                                      "--out", str(out)])
    assert arm.main() == 0
    return (json.loads(out.read_text(encoding="utf-8")),
            out.with_suffix(".md").read_text(encoding="utf-8"))


def test_main_MEASURES_the_stacking_question_when_EVERY_channel_is_present(tmp_path, monkeypatch):
    """The complement of `test_main_SKIPS_an_optional_combo...`: with nothing absent there is nothing
    to skip, so every superset combo must get its OWN statistic against the 2-channel hybrid and the
    stacking verdict must be a MEASURED one rather than the refusal."""
    art, _md = _complete_run(tmp_path, monkeypatch)
    a = art["arms"]
    assert "_skipped_combos" not in a, "nothing is absent, so nothing may be recorded as skipped"
    assert "_short_arm" not in a

    sup = arm.superset_combos()
    assert sup
    for combo in sup:
        assert a[combo]["n"] == art["n_seeds"], f"{combo} scored a seed subset"
        # BOTH statistics exist and are separate numbers: one against the scalar, one against the
        # 2-channel. A single key for both would make the stacking number overwrite the published one.
        assert arm.paired_key(combo) in a
        key = arm.stacking_pair_key(combo)
        assert key in a, f"{combo} has no statistic against the stacking comparator"
        assert a[key]["n"] == art["n_seeds"]
        assert a[key] != a[arm.paired_key(combo)]

    assert art["stacking_verdict"] != "INDETERMINATE_STACKING_NOT_MEASURED"
    assert art["stacking_verdict"] in PREREGISTERED["stacking"]["verdicts"]
    assert art["stacking_verdict_detail"].strip()
    assert set(art["stacking_per_combo"]) == set(sup)
    for row in art["stacking_per_combo"].values():
        assert row["verdict"] in PREREGISTERED["stacking"]["verdicts"]
        assert row["detail"].strip()


def test_the_MEASURED_stacking_verdict_REACHES_the_markdown_with_its_per_combo_block(
        tmp_path, monkeypatch):
    """A second question carried only in the JSON is not a disclosure -- and the top-level `verdict`
    line says nothing about stacking, which is what makes this section load-bearing. The per-combo
    table is the part the honest-limits entry tells the reader to prefer over the aggregate."""
    art, md = _complete_run(tmp_path, monkeypatch)
    assert f"**STACKING VERDICT: `{art['stacking_verdict']}`**" in md
    assert PREREGISTERED["stacking"]["comparator"] in md
    # the per-combo block, one row per combo that ran, with its own verdict beside it
    for combo, row in art["stacking_per_combo"].items():
        assert f"| `{combo}` | `{row['verdict']}` |" in md
    # the prediction and the scope were frozen before the run, so both must travel with the result
    assert PREREGISTERED["stacking"]["predicted_outcome"] in md
    assert PREREGISTERED["stacking"]["negative_scope"] in md
    # and the parent question's answer is still its own, separately labelled line
    assert f"**VERDICT: `{art['verdict']}`**" in md
    assert art["verdict"] in PREREGISTERED["verdicts"]


def test_the_honest_limits_record_the_MULTIPLE_COMPARISON_exposure_and_the_RANK_ONLY_scope(
        tmp_path, monkeypatch):
    """Two costs the stacking question introduces, asserted on what the run SHIPS rather than on a
    literal in the module: k combos against one bar inflates the chance one clears it by noise, and
    `orthogonality` is a rank relationship, not a claim that the channels share information."""
    art, md = _complete_run(tmp_path, monkeypatch)
    lims = " ".join(art["honest_limits"]).lower()
    assert "multiple-comparison" in lims
    assert "per-combo" in lims, "the exposure is named but the reader is not told what to read"
    assert "rank" in lims and "prediction" in lims
    for lim in art["honest_limits"]:
        assert lim in md, "an honest limit that reaches the JSON but not the markdown"


# --------------------------------------------------------------------------------------------------
# the NEW verdict key -- its aggregate PRECEDENCE, its FALLBACK-path detail, and its reachability
# through the full main() path. Every other test above reaches it by calling `stacking_verdict_for`
# directly, which cannot show that it survives the artifact write and the renderer.
# --------------------------------------------------------------------------------------------------
_NEW = "INDETERMINATE_GAIN_NOT_REPRODUCIBLE"


def test_the_aggregate_PRECEDENCE_is_a_strict_ladder_with_the_new_key_BELOW_weak_stack():
    """Inserting the new key into the precedence tuple fixes one defect and risks its mirror.

    Placed too HIGH it demotes a REPRODUCIBLE gain to indeterminate; left out entirely (the bug
    `test_an_ALL_unreproducible_positive_arm...` pins) it promotes positive noise to a measured
    NEGATIVE. Only the ORDER rules both out, and nothing above asserts on the order -- the existing
    mixed-verdict fixtures pair `WEAK_STACK` with `NO_STACK` only, which the ladder's other three
    rungs never touch. Every ordered pair is checked, so a single swap anywhere fails here.
    """
    ladder = ["STACKS", "WEAK_STACK", _NEW, "NO_STACK"]
    stat = {
        "STACKS": _p(0.05, 10),
        "WEAK_STACK": _p(0.01, 10),
        _NEW: _p(0.01, PREREGISTERED["stacking"]["min_seeds_winning"] - 1),
        "NO_STACK": _p(-0.01, 2),
    }
    # the fixture must really produce the verdict it is named for, or the ladder is checked on
    # mislabelled inputs and would pass under a reordering
    for want, s in stat.items():
        assert arm.stacking_verdict_for(s)[0] == want, f"fixture for {want} does not produce it"

    sup = arm.superset_combos()
    assert len(sup) >= 2, "need two superset combos to express a two-verdict arm"
    for i, higher in enumerate(ladder):
        for lower in ladder[i + 1:]:
            arms = {arm.stacking_pair_key(sup[0]): stat[higher],
                    arm.stacking_pair_key(sup[1]): stat[lower]}
            v, _d, per = arm.stacking_verdict(arms)
            assert v == higher, f"{higher} beside {lower} aggregated to {v}, not {higher}"
            # and the loser is never dropped from the per-combo block by the precedence choice
            assert per[sup[1]]["verdict"] == lower


def test_the_NO_STACK_FALLBACK_detail_also_names_EVERY_combo_and_in_the_SAME_format():
    """`_all_named()` exists so the hit path and the fallback cannot drift apart again, and only the
    HIT path is asserted above. An all-negative arm takes the fallback `return`, so if that one line
    ever names fewer combos -- or names them differently -- the JSON reader sees a different shape
    depending on which branch fired, which is the drift the helper was introduced to prevent.
    """
    sup = arm.superset_combos()
    assert len(sup) >= 2, "a single-combo arm cannot show that EVERY combo is named"
    lose = {"median": -0.0039, "wins": 4, "n": 10, "ci95": [-0.0095, 0.0026]}
    arms = {arm.stacking_pair_key(c): dict(lose) for c in sup}
    v, detail, per = arm.stacking_verdict(arms)
    assert v == "NO_STACK", "fixture does not take the fallback return"
    for c in sup:
        assert f"{c} [{per[c]['verdict']}]: {per[c]['detail']}" in detail, (
            f"{c} is in stacking_per_combo but not named verbatim on the fallback path")

    # SAME FORMAT as the hit path: promote one combo and every OTHER combo's rendering must be
    # byte-identical to what the fallback produced for it.
    arms[arm.stacking_pair_key(sup[0])] = _p(0.05, 10)
    v_hit, detail_hit, per_hit = arm.stacking_verdict(arms)
    assert v_hit == "STACKS", "fixture does not take the hit return"
    for c in sup[1:]:
        frag = f"{c} [{per[c]['verdict']}]: {per[c]['detail']}"
        assert frag in detail and frag in detail_hit, f"{c} renders differently per branch"
    assert detail.index(sup[1]) >= 0 and sorted(per) == sorted(per_hit)


def test_main_REACHES_the_new_key_end_to_end_and_it_LANDS_in_both_surfaces(tmp_path, monkeypatch):
    """The new key through the FULL path: main() -> artifact JSON -> `render_md`.

    `_complete_run` runs 2 seeds against a bar asking for 8, so a positive stacking gain cannot be
    reproducible -- which makes this fixture the one that exercises the new branch end to end. Every
    existing main() assertion on `stacking_verdict` is `in PREREGISTERED[...]["verdicts"]`, satisfied
    by any of the five keys, and the markdown assertion interpolates the artifact's own value; so
    nothing above would notice the branch collapsing back into `WEAK_STACK` or `NO_STACK`.
    """
    art, md = _complete_run(tmp_path, monkeypatch)
    k = PREREGISTERED["stacking"]["min_seeds_winning"]
    assert art["n_seeds"] < k, "fixture would have to be re-chosen; it no longer under-runs the bar"

    unreproducible_positives = [
        c for c in arm.superset_combos()
        if (st := art["arms"].get(arm.stacking_pair_key(c))) and st["median"] > 0 and st["wins"] < k
    ]
    assert unreproducible_positives, "fixture produced no positive-but-unreproducible gain"
    for c in unreproducible_positives:
        assert art["stacking_per_combo"][c]["verdict"] == _NEW, (
            f"{c} has a positive gain on {art['arms'][arm.stacking_pair_key(c)]['wins']} of "
            f"{art['n_seeds']} seeds and did not get the unreproducible verdict")

    # the aggregate: no combo cleared the margin AND the seed bar, and none is a measured negative,
    # so the arm-level answer must be the new key -- never `NO_STACK`, which would report a measured
    # negative for a set of positive gains.
    if len(unreproducible_positives) == len(art["stacking_per_combo"]):
        assert art["stacking_verdict"] == _NEW
    assert art["stacking_verdict"] != "NO_STACK"

    # and BOTH surfaces carry it. The key is asserted as a LITERAL here, not interpolated from the
    # artifact, which is what makes this a pin rather than a tautology.
    assert f"**STACKING VERDICT: `{_NEW}`**" in md
    for c in unreproducible_positives:
        assert f"| `{c}` | `{_NEW}` |" in md, f"{c}'s row in the per-combo table lost the verdict"
    assert "did NOT reproduce" in md
    assert "neither a measured stack nor a measured negative" in md
