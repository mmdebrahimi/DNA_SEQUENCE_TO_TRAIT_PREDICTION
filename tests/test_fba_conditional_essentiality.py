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


# ---- continuous readout: is the signal absent, or is the cutoff discarding it? ----

def test_continuous_readout_recovers_a_perfect_ranking_the_binary_cutoff_would_miss():
    """The case the metric exists to detect: ratios that order the conditions correctly but never fall
    below the 1% cutoff. AUROC must see it; the deployed threshold must not."""
    from dna_decode.fba.conditional_essentiality import continuous_readout

    r = _rec("b1", [True, False, False, False])
    keys = sorted(CONDITIONS)
    ratios = {keys[0]: {"b1": 0.20}, keys[1]: {"b1": 0.90},
              keys[2]: {"b1": 0.95}, keys[3]: {"b1": 0.99}}
    got = continuous_readout([r], ratios)
    assert got["auroc"] == pytest.approx(1.0)          # ranking is perfect
    assert got["deployed_confusion"]["tp"] == 0        # ...and the deployed cutoff catches none of it


def test_continuous_readout_reports_auroc_half_when_the_ratio_is_flat():
    """A flat ratio carries no conditional information; AUROC must say so rather than flatter it."""
    from dna_decode.fba.conditional_essentiality import continuous_readout

    r = _rec("b1", [True, False, False, False])
    ratios = {c: {"b1": 0.5} for c in CONDITIONS}
    assert continuous_readout([r], ratios)["auroc"] == pytest.approx(0.5)


def test_the_oracle_threshold_is_labelled_as_an_upper_bound():
    """Load-bearing honesty rail, same shape as the Track B deltaG arm: the best threshold is fitted ON
    the evaluation set and must never be presented as a deployable number."""
    from dna_decode.fba.conditional_essentiality import continuous_readout

    r = _rec("b1", [True, False, False, False])
    keys = sorted(CONDITIONS)
    got = continuous_readout([r], {keys[0]: {"b1": 0.2}, keys[1]: {"b1": 0.9},
                                   keys[2]: {"b1": 0.9}, keys[3]: {"b1": 0.9}})
    assert "upper bound" in got["oracle_note"].lower()
    assert got["oracle_mcc"] >= got["deployed_mcc"]     # an oracle can never be worse


def test_continuous_readout_is_degenerate_safe():
    """All-essential or all-dispensable cells give no AUROC; report it, do not divide by zero."""
    from dna_decode.fba.conditional_essentiality import continuous_readout

    r = _rec("b1", [True, True, True, True])           # not two-sided -> no cells
    assert continuous_readout([r], {c: {"b1": 0.5} for c in CONDITIONS})["auroc"] is None


# ---- deployable threshold (the honest version of the oracle) ----

def test_deployable_threshold_splits_by_gene_not_by_cell():
    """Load-bearing against leakage: the four cells of one gene share its ratio profile, so splitting
    cells would put a gene's own data in both train and test."""
    from dna_decode.fba.conditional_essentiality import deployable_threshold

    keys = sorted(CONDITIONS)
    recs = [_rec(f"b{i}", [True, False, False, False]) for i in range(10)]
    ratios = {c: {f"b{i}": (0.2 if c == keys[0] else 0.9) for i in range(10)} for c in keys}
    got = deployable_threshold(recs, ratios, n_folds=5)
    assert got["n_cells"] == 40                     # 10 genes x 4 conditions
    assert got["n_folds"] == 5
    assert len(got["thresholds_per_fold"]) == 5


def test_deployable_threshold_recovers_a_clean_separable_signal_out_of_sample():
    """When the ratio separates the classes perfectly, a held-out fit must find it."""
    from dna_decode.fba.conditional_essentiality import deployable_threshold

    keys = sorted(CONDITIONS)
    recs = [_rec(f"b{i}", [True, False, False, False]) for i in range(10)]
    ratios = {c: {f"b{i}": (0.1 if c == keys[0] else 0.9) for i in range(10)} for c in keys}
    assert deployable_threshold(recs, ratios, n_folds=5)["held_out_mcc"] == pytest.approx(1.0)


def test_deployable_threshold_does_not_invent_signal_from_noise():
    """A ratio identical for essential and dispensable cells must NOT produce a positive held-out MCC."""
    from dna_decode.fba.conditional_essentiality import deployable_threshold

    recs = [_rec(f"b{i}", [True, False, False, False]) for i in range(10)]
    ratios = {c: {f"b{i}": 0.5 for i in range(10)} for c in CONDITIONS}
    assert deployable_threshold(recs, ratios, n_folds=5)["held_out_mcc"] <= 0.0


def test_deployable_threshold_is_degenerate_safe():
    from dna_decode.fba.conditional_essentiality import deployable_threshold

    assert deployable_threshold([], {}, n_folds=5)["held_out_mcc"] is None


# ---- generalisation to an arbitrary condition set (Fitness Browser, 2026-08-12) ----

def test_metrics_accept_a_condition_set_beyond_the_four_shipped_media():
    """The 4-media CONDITIONS is a default, not a hard-coding. The Fitness Browser path scores 25 carbon
    sources through the SAME functions."""
    from dna_decode.fba.conditional_essentiality import constant_baselines

    keys = ("carbA", "carbB", "carbC")
    r = GeneRecord("b1", "b1", {"carbA": True, "carbB": False, "carbC": False}, {}, True)
    pred = {c: {"b1": (c == "carbA")} for c in keys}
    assert switch_accuracy([r], pred, conditions=keys)["exact_set_match"] == 1
    assert constant_baselines([r], conditions=keys)["always_dispensable"]["exact_set_match"] == 0


def test_pattern_distribution_MUST_use_the_callers_conditions():
    """THE bug this pins. Defaulting to the 4 media against a 25-carbon-source prediction made every
    lookup miss, so every gene read as '....' and the run reported a FALSE '100% constant across 1 shape'
    -- which contradicted a positive exact-set count in the very same run. A contradiction between two
    metrics of one run is a bug signal, not a finding."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    keys = ("carbA", "carbB", "carbC")
    r = GeneRecord("b1", "b1", {"carbA": True, "carbB": False, "carbC": False}, {}, True)
    pred = {c: {"b1": (c == "carbA")} for c in keys}

    got = pattern_distribution([r], pred, conditions=keys)
    assert got["n_constant_pattern"] == 0            # it is NOT constant -- it matches exactly
    assert got["patterns"] == {"E..": 1}

    stale = pattern_distribution([r], pred)          # default 4-media keys -> every lookup misses
    assert stale["n_constant_pattern"] == 1          # the false "constant" the bug produced
    assert got["n_constant_pattern"] != stale["n_constant_pattern"]


def test_constant_detection_is_LENGTH_AGNOSTIC():
    """The SECOND hardcoded-4 bug in this one function. The literal ("....","EEEE") test meant a
    25-condition all-dispensable pattern matched NEITHER, so a 25-carbon-source run reported
    '0.0% constant' when the true figure was 184/217 = 84.8%. Generalising the condition KEYS was not
    enough -- the constant TEST was hardcoded too."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    keys = tuple(f"c{i}" for i in range(25))
    r = GeneRecord("b1", "b1", {k: (k == "c0") for k in keys}, {}, True)

    all_disp = {k: {"b1": False} for k in keys}          # 25 dots
    assert pattern_distribution([r], all_disp, conditions=keys)["n_constant_pattern"] == 1
    all_ess = {k: {"b1": True} for k in keys}            # 25 E's
    assert pattern_distribution([r], all_ess, conditions=keys)["n_constant_pattern"] == 1
    varying = {k: {"b1": (k == "c3")} for k in keys}
    assert pattern_distribution([r], varying, conditions=keys)["n_constant_pattern"] == 0


def test_the_four_media_constant_numbers_are_UNAFFECTED_by_the_fix():
    """Blast-radius check: at 4 conditions the old literal test and the new length-agnostic test agree,
    so the shipped 4-media result (94.0% constant) needs no revision."""
    from dna_decode.fba.conditional_essentiality import pattern_distribution

    records = load_labels()
    got = pattern_distribution(records, {c: {r.gene_id: r.paper_fba[c] for r in records}
                                         for c in CONDITIONS})
    assert got["n_constant_pattern"] == 62          # the published iJO1366 figure, unchanged
