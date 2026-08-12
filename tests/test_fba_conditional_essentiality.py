"""Conditional gene essentiality -- pure logic + the committed gold standard (no cobra, no network)."""
from __future__ import annotations

import pytest

from dna_decode.fba.conditional_essentiality import (
    CONDITIONS,
    GeneRecord,
    conditionally_essential_genes,
    confusion_from_calls,
    load_labels,
    mcc,
    switch_accuracy,
)


def _rec(gid, exp, fba=None, flag=False):
    exp_d = dict(zip(sorted(CONDITIONS), exp, strict=True))
    fba_d = dict(zip(sorted(CONDITIONS), fba or exp, strict=True))
    return GeneRecord(gid, gid, exp_d, fba_d, flag)


# ---- the two-sided subset ----

def test_conditionally_essential_requires_BOTH_an_essential_and_a_dispensable_condition():
    """The whole point of the metric. Always-essential and never-essential genes are NOT two-sided and
    must be excluded, or a model that calls everything dispensable scores well by default."""
    always = _rec("b1", [True, True, True, True])
    never = _rec("b2", [False, False, False, False])
    switch = _rec("b3", [True, False, False, False])

    got = {r.gene_id for r in conditionally_essential_genes([always, never, switch])}
    assert got == {"b3"}


def test_conditional_subset_is_recomputed_not_inherited_from_the_supplement_flag():
    """The supplement ships its own YES/NO flag. Recomputing from the experimental columns means a
    disagreement is VISIBLE rather than silently inherited."""
    lying = _rec("b9", [False, False, False, False], flag=True)   # flagged YES, but never essential
    assert conditionally_essential_genes([lying]) == []


# ---- the switch metric ----

def test_exact_set_match_requires_the_whole_pattern_not_just_some_cells():
    r = _rec("b1", [True, False, False, False])          # essential in exactly one medium
    keys = sorted(CONDITIONS)

    right = {c: {"b1": (c == keys[0])} for c in keys}
    assert switch_accuracy([r], right)["exact_set_match"] == 1

    off_by_one = {c: {"b1": c in (keys[0], keys[1])} for c in keys}
    got = switch_accuracy([r], off_by_one)
    assert got["exact_set_match"] == 0                    # strict: pattern is wrong
    assert got["per_condition_agreement"] == pytest.approx(0.75)   # lenient: 3 of 4 cells right


def test_a_model_that_calls_a_gene_essential_everywhere_scores_zero_on_the_switch():
    """A single-condition metric can look fine while the model has no conditional resolution at all.
    This is what the conditional metric exists to expose -- and it is what the paper's own iJO1366 does
    (exact-set 4/68)."""
    r = _rec("b1", [True, False, False, False])
    all_essential = {c: {"b1": True} for c in CONDITIONS}
    assert switch_accuracy([r], all_essential)["exact_set_match"] == 0

    none_essential = {c: {"b1": False} for c in CONDITIONS}
    assert switch_accuracy([r], none_essential)["exact_set_match"] == 0


def test_switch_metrics_are_none_rather_than_a_divide_by_zero_on_an_empty_subset():
    got = switch_accuracy([_rec("b1", [False, False, False, False])], {})
    assert got["n_conditionally_essential"] == 0
    assert got["exact_set_match_rate"] is None


def test_a_missing_condition_in_the_prediction_counts_as_not_essential_not_a_crash():
    """A model that could not be scored in one medium (e.g. no growth) must degrade to a wrong call
    there, never take the whole run down."""
    r = _rec("b1", [True, False, False, False])
    assert switch_accuracy([r], {})["exact_set_match"] == 0


# ---- confusion / MCC ----

def test_confusion_counts_over_the_shared_key_set_only():
    cm = confusion_from_calls({"a": True, "b": False, "c": True}, {"a": True, "b": True})
    assert (cm["tp"], cm["fp"], cm["fn"], cm["tn"], cm["n"]) == (1, 1, 0, 0, 2)


def test_mcc_is_zero_not_a_crash_when_the_prediction_is_single_class():
    """A model calling nothing essential has a zero margin; MCC must report 0.0, not raise."""
    assert mcc({"tp": 0, "fp": 0, "fn": 10, "tn": 90, "n": 100}) == 0.0


def test_mcc_is_one_for_a_perfect_call_and_negative_when_inverted():
    assert mcc({"tp": 10, "fp": 0, "fn": 0, "tn": 90, "n": 100}) == pytest.approx(1.0)
    assert mcc({"tp": 0, "fp": 10, "fn": 90, "tn": 0, "n": 100}) < 0


# ---- the committed gold standard ----

def test_the_committed_gold_standard_parses_and_carries_the_two_sided_signal():
    """Load-bearing: this file IS the substrate. 1,075 genes, 68 conditionally essential."""
    records = load_labels()
    assert len(records) == 1075
    assert len(conditionally_essential_genes(records)) == 68


def test_the_supplement_flag_and_the_recomputed_subset_agree_on_the_real_data():
    """They do agree today (68 = 68). If a future edit breaks the parse, this catches it."""
    records = load_labels()
    assert sum(1 for r in records if r.conditionally_essential) == 68


def test_every_condition_has_both_essential_and_dispensable_genes():
    """A condition where every gene is dispensable would make its column uninformative."""
    records = load_labels()
    for c in CONDITIONS:
        calls = [r.experimental[c] for r in records]
        assert 0 < sum(calls) < len(calls), c


def test_the_papers_own_fba_columns_are_present_and_are_not_identical_to_the_labels():
    """The FBA columns are a reproduction GATE, so they must be loaded -- and they must differ from the
    experimental columns, or the gate would be vacuous."""
    records = load_labels()
    disagreements = sum(1 for r in records for c in CONDITIONS
                        if r.experimental[c] != r.paper_fba[c])
    assert disagreements > 0


# ---- null controls (the number is meaningless without them) ----

def test_constant_predictors_score_zero_exact_matches_but_a_high_per_cell_rate():
    """THE reason the switch metric needs a null. Predicting 'dispensable everywhere' gets 0 patterns
    right yet still agrees on ~56% of cells, because most conditionally-essential genes are essential in
    only 1-2 of the 4 media. Measured: 0.5588 -- against which the models' 0.57 is ~1 point of lift."""
    from dna_decode.fba.conditional_essentiality import constant_baselines

    nulls = constant_baselines(load_labels())
    assert nulls["always_dispensable"]["exact_set_match"] == 0
    assert nulls["always_essential"]["exact_set_match"] == 0
    assert nulls["always_dispensable"]["per_condition_agreement"] == pytest.approx(0.5588, abs=1e-3)
    assert nulls["always_essential"]["per_condition_agreement"] == pytest.approx(0.4412, abs=1e-3)


def test_the_two_nulls_are_complementary_per_cell():
    """Sanity pin on the metric itself: every cell is either essential or not, so the two constant
    predictors' per-cell rates must sum to 1."""
    from dna_decode.fba.conditional_essentiality import constant_baselines

    n = constant_baselines(load_labels())
    assert (n["always_essential"]["per_condition_agreement"]
            + n["always_dispensable"]["per_condition_agreement"]) == pytest.approx(1.0)


# ---- pattern shape: WHY the switch score is low ----

def test_pattern_distribution_flags_a_constant_predictor_as_constant():
    """The mechanism behind the low switch score. A model predicting the same call in every medium is
    'not switching at all', and that must be countable rather than inferred."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    r = _rec("b1", [True, False, False, False])
    got = pattern_distribution([r], {c: {"b1": False} for c in CONDITIONS})
    assert got["n_constant_pattern"] == 1
    assert got["constant_pattern_fraction"] == pytest.approx(1.0)
    assert got["patterns"] == {"....": 1}


def test_true_patterns_are_never_constant_by_definition():
    """A conditionally-essential gene has >=1 E and >=1 N, so its own pattern can never be constant.
    If this ever fires, the two-sided subset is being computed wrong."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    got = pattern_distribution(load_labels())
    assert got["n_constant_pattern"] == 0
    assert got["n_genes"] == 68


def test_the_papers_own_model_is_constant_on_91_percent_of_the_switching_genes():
    """The headline diagnostic, pinned: iJO1366 predicts '....' or 'EEEE' for 62 of the 68 genes whose
    essentiality actually depends on the medium."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    records = load_labels()
    got = pattern_distribution(records, {c: {r.gene_id: r.paper_fba[c] for r in records}
                                         for c in CONDITIONS})
    assert got["n_constant_pattern"] == 62
    assert got["constant_pattern_fraction"] == pytest.approx(0.9118, abs=1e-3)


def test_the_experimental_patterns_are_genuinely_diverse():
    """If the true patterns collapsed to one or two shapes the metric would be near-degenerate."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    assert pattern_distribution(load_labels())["n_distinct_patterns"] >= 10
