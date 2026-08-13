"""Why the model does not switch — pure diagnosis logic (no solver, no feba.db)."""
from __future__ import annotations

from dna_decode.fba.conditional_essentiality import GeneRecord
from scripts.fba_constant_gene_diagnostic import classify_miss, diagnose, verdict_for

KEYS = ("a", "b", "c")


def _rec(gid, exp):
    return GeneRecord(gid, gid, dict(zip(KEYS, exp, strict=True)), {}, True)


def test_a_flat_ratio_is_a_model_problem_not_a_threshold_one():
    """Ratio 1.0 means the deletion changed NOTHING -- the model has an alternative route. No readout
    change can recover this."""
    assert classify_miss(1.0) == "flat"
    assert classify_miss(0.9999999) == "flat"


def test_a_materially_depressed_ratio_is_a_readout_problem():
    """The model SAW a large growth defect and the 1% binary cutoff threw it away."""
    assert classify_miss(0.5) == "material"
    assert classify_miss(0.05) == "near_threshold"      # nearly called it


def test_a_slight_defect_is_not_counted_as_readout_recoverable():
    """A 2% growth dip is not a missed essentiality call; counting it would inflate the recoverable
    fraction and overstate how much a graded metric buys."""
    assert classify_miss(0.98) == "slight"


def test_the_three_prediction_strata_partition_the_gene_set():
    subset = [_rec("gDisp", (True, False, False)),
              _rec("gEss", (True, False, False)),
              _rec("gCommit", (False, True, False))]
    calls = {"a": {"gDisp": False, "gEss": True, "gCommit": False},
             "b": {"gDisp": False, "gEss": True, "gCommit": True},
             "c": {"gDisp": False, "gEss": True, "gCommit": False}}
    ratios = {c: {"gDisp": 1.0, "gEss": 0.0, "gCommit": 1.0} for c in KEYS}
    d = diagnose(subset, calls, ratios, KEYS)
    assert d["n_predicted_all_dispensable"] == 1
    assert d["n_predicted_all_essential"] == 1
    assert d["n_committed"] == 1
    assert (d["n_predicted_all_dispensable"] + d["n_predicted_all_essential"]
            + d["n_committed"]) == len(subset)


def test_misses_are_counted_only_where_TRUTH_says_essential():
    """A dispensable-everywhere prediction is only WRONG in the conditions the gene is truly essential;
    counting its correct dispensable calls as misses would invent a deficit."""
    subset = [_rec("g1", (True, False, False))]           # truly essential in 'a' only
    calls = {c: {"g1": False} for c in KEYS}
    ratios = {"a": {"g1": 0.4}, "b": {"g1": 1.0}, "c": {"g1": 1.0}}
    d = diagnose(subset, calls, ratios, KEYS)
    assert d["missed_essential_cells"] == 1               # not 3
    assert d["missed_cells_by_cause"] == {"material": 1}


def test_a_gene_flat_in_every_true_essential_condition_is_flagged():
    subset = [_rec("gFlat", (True, True, False)), _rec("gMixed", (True, True, False))]
    calls = {c: {"gFlat": False, "gMixed": False} for c in KEYS}
    ratios = {"a": {"gFlat": 1.0, "gMixed": 1.0}, "b": {"gFlat": 1.0, "gMixed": 0.3},
              "c": {"gFlat": 1.0, "gMixed": 1.0}}
    d = diagnose(subset, calls, ratios, KEYS)
    assert d["n_genes_flat_in_EVERY_true_essential_condition"] == 1   # gFlat only
    assert d["example_flat_genes"] == ["gFlat"]


def test_overcalls_are_split_by_cause_on_the_conditions_truth_says_dispensable():
    subset = [_rec("gOver", (True, False, False))]
    calls = {c: {"gOver": True} for c in KEYS}            # predicted essential everywhere
    ratios = {"a": {"gOver": 0.0}, "b": {"gOver": None}, "c": {"gOver": 0.001}}
    d = diagnose(subset, calls, ratios, KEYS)
    assert d["overcalled_dispensable_cells_by_cause"] == {"infeasible": 1, "sub_threshold": 1}


def test_the_verdict_tracks_which_diagnosis_dominates():
    assert verdict_for({"readout_recoverable_fraction": 0.8}) == \
        "DEFICIT_IS_MOSTLY_A_READOUT_PROBLEM"
    assert verdict_for({"readout_recoverable_fraction": 0.3}) == \
        "DEFICIT_IS_MIXED_MODEL_AND_READOUT"
    assert verdict_for({"readout_recoverable_fraction": 0.05}) == \
        "DEFICIT_IS_A_MODEL_PROBLEM_NOT_A_READOUT_ONE"
    assert verdict_for({"readout_recoverable_fraction": None}) == "NO_MISSES_TO_EXPLAIN"
