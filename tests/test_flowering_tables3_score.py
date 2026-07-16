"""Tests for the flowering cell's Table S3 scorer.

Pins the honesty gates, not the headline number: the gates are the load-bearing part. Pure/offline tests
run everywhere; the real-data tests skip if the CC-BY Table S3 is absent.
"""
from __future__ import annotations

import pytest

from scripts.flowering_tables3_score import (
    ScoringError,
    TABLE_S3,
    _directional,
    confusion,
    load_table_s3,
    phenotype_attrition,
    predict,
    score,
)

pytestmark = pytest.mark.filterwarnings("ignore")


def _row(deleterious: str, ft16: str, group: str = "g1", allele: str = "a001") -> dict:
    return {"deleterious_allele": deleterious, "FT16_mean": ft16, "group": group, "allele_group": allele}


# ---- gate 2: the null baseline ---------------------------------------------------------------------------

def test_confusion_reports_the_best_constant_predictor_not_just_accuracy():
    # 8 late / 2 early; a constant 'late' predictor scores 0.8 -- any real call must beat THAT, not 0.5.
    pairs = [("late", "late")] * 8 + [("late", "early")] * 2
    m = confusion(pairs)
    assert m["accuracy"] == 0.8
    assert m["null_accuracy"] == 0.8 and m["null_call"] == "late"


def test_single_class_set_is_flagged_degenerate_and_never_scored():
    m = confusion([("late", "late"), ("early", "late")])
    assert m["degenerate"] is True
    assert m["specificity"] is None       # no early observed -> undefined, not silently 0.0


def test_verdict_is_fails_null_when_the_rule_does_not_beat_a_constant():
    # An anti-correlated cohort: functional-FRI accessions are early, LoF ones are late (the rule inverted).
    # 8 early / 2 late -> a constant 'early' scores 0.8 while the rule scores 0.0.
    rows = [_row("FALSE", "10")] * 8 + [_row("TRUE", "90")] * 2
    out = score(rows)
    assert out["pooled"]["accuracy"] == 0.0
    assert out["pooled"]["null_accuracy"] == 0.8
    assert out["verdict"] == "FAILS_NULL_BASELINE"


# ---- gate 5: non-random phenotype attrition -------------------------------------------------------------

def test_na_phenotype_rows_are_excluded_from_scoring_not_coerced():
    rows = [_row("TRUE", "10"), _row("FALSE", "90"), _row("FALSE", "NA"), _row("TRUE", "")]
    assert score(rows)["pooled"]["n"] == 2


def test_attrition_detects_genotype_skewed_dropout():
    # Every dropped row is functional-FRI while the cohort is half deleterious -> skewed.
    rows = [_row("TRUE", "10"), _row("TRUE", "12"), _row("FALSE", "NA"), _row("FALSE", "NA")]
    a = phenotype_attrition(rows)
    assert a["n_dropped_no_ft16"] == 2
    assert a["deleterious_rate_among_dropped"] == 0.0
    assert a["dropout_is_genotype_skewed"] is True


def test_attrition_reports_unskewed_dropout_as_such():
    rows = [_row("TRUE", "10"), _row("FALSE", "90"), _row("TRUE", "NA"), _row("FALSE", "NA")]
    assert phenotype_attrition(rows)["dropout_is_genotype_skewed"] is False


# ---- gate 1: the cell's habit vocabulary is mapped explicitly, never guessed ------------------------------

def test_predict_maps_cell_habits_to_scoring_labels():
    assert predict(_row("TRUE", "10")) == "early"      # FRI LoF   -> summer_annual_early
    assert predict(_row("FALSE", "90")) == "late"      # FRI intact -> winter_annual_late (FLC assumed strong)


def test_predict_raises_rather_than_bucket_an_unmappable_habit(monkeypatch):
    import scripts.flowering_tables3_score as mod

    class _Abstain:
        habit = "ABSTAIN"

    monkeypatch.setattr(mod, "call_flowering_habit", lambda *_: _Abstain())
    with pytest.raises(ScoringError, match="ABSTAIN"):
        predict(_row("TRUE", "10"))


# ---- gate 3: population structure ------------------------------------------------------------------------

def test_group_weighting_excludes_single_class_and_tiny_groups():
    rows = ([_row("TRUE", "10", "big")] * 10 + [_row("FALSE", "90", "big")] * 10
            + [_row("TRUE", "10", "oneclass")] * 12          # early only -> unscorable
            + [_row("TRUE", "10", "tiny")] * 2)              # n<10 -> unscorable
    out = score(rows)
    g = out["by_structure_group"]
    assert g["big"]["scorable"] is True
    assert g["oneclass"]["scorable"] is False and g["oneclass"]["degenerate"] is True
    assert g["tiny"]["scorable"] is False
    assert out["group_weighted"]["n_groups_scorable"] == 1


def test_group_weighting_gives_each_group_one_vote():
    """The whole point of gate 3: an over-sampled group must not carry the metric.

    'big' (n=100) the rule aces; 'small' (n=20) it gets exactly backwards. Pooled is carried by 'big' and
    looks strong; one-vote-per-group shows the rule is really 1-for-2 across ancestries.
    """
    rows = ([_row("TRUE", "10", "big")] * 50 + [_row("FALSE", "90", "big")] * 50
            + [_row("TRUE", "90", "small")] * 10 + [_row("FALSE", "10", "small")] * 10)
    out = score(rows)
    assert out["by_structure_group"]["big"]["accuracy"] == 1.0
    assert out["by_structure_group"]["small"]["accuracy"] == 0.0
    assert out["pooled"]["accuracy"] > 0.8                 # pooled: carried by the over-sampled group
    assert out["group_weighted"]["n_groups_scorable"] == 2
    assert out["group_weighted"]["mean_accuracy"] == 0.5    # one vote each: the honest figure
    assert out["group_weighted"]["mean_accuracy"] < out["pooled"]["accuracy"]


# ---- the directional diagnostic --------------------------------------------------------------------------

def test_directional_flags_false_positives_clustered_at_the_threshold():
    # FPs sitting right at the cut = a thresholding artifact, not a mechanism.
    rows = [_row("FALSE", "49")] * 5 + [_row("FALSE", "80")] * 5 + [_row("TRUE", "10")] * 5
    d = _directional(rows, 50.0)
    assert d["false_positive_spread"]["clustered_at_threshold"] is True


def test_directional_flags_false_positives_spread_far_below_threshold():
    rows = [_row("FALSE", "5")] * 5 + [_row("FALSE", "80")] * 5 + [_row("TRUE", "10")] * 5
    d = _directional(rows, 50.0)
    assert d["false_positive_spread"]["clustered_at_threshold"] is False


# ---- substrate integrity ---------------------------------------------------------------------------------

def test_load_refuses_a_table_that_disagrees_with_the_papers_own_counts(tmp_path):
    f = tmp_path / "bad.tsv"
    f.write_text("accession_id\tdeleterious_allele\tallele_group\tgroup\tFT16_mean\n1\tTRUE\ta001\tg\t50\n",
                 encoding="utf-8")
    with pytest.raises(ScoringError, match="stated counts"):
        load_table_s3(f)


def test_load_gives_an_actionable_error_when_the_browser_only_table_is_absent(tmp_path):
    with pytest.raises(ScoringError, match="browser-only"):
        load_table_s3(tmp_path / "missing.tsv")


# ---- real data (skips without the CC-BY table) -----------------------------------------------------------

@pytest.mark.skipif(not TABLE_S3.exists(), reason="Table S3 not present (browser-only fetch)")
def test_real_table_s3_matches_the_papers_stated_counts():
    rows = load_table_s3()
    assert len(rows) == 1017
    assert sum(1 for r in rows if r["deleterious_allele"] == "TRUE") == 245
    assert len({r["allele_group"] for r in rows}) == 103


@pytest.mark.skipif(not TABLE_S3.exists(), reason="Table S3 not present (browser-only fetch)")
def test_real_run_is_scored_and_the_asymmetry_holds():
    out = score(load_table_s3())
    assert out["verdict"] == "SCORED_BEATS_NULL"
    d = out["directional"]
    # The finding: FRI-LoF -> early is strong; FRI-functional -> late is weak (necessary, not sufficient).
    assert d["negative_direction"]["precision"] > 0.9
    assert d["positive_direction"]["precision"] < 0.75
    # ...and the FPs are NOT a thresholding artifact.
    assert d["false_positive_spread"]["clustered_at_threshold"] is False


@pytest.mark.skipif(not TABLE_S3.exists(), reason="Table S3 not present (browser-only fetch)")
def test_real_run_population_structure_shrinks_the_advantage():
    """The headline honesty claim: conditioning on ancestry collapses most of the pooled advantage."""
    out = score(load_table_s3())
    pooled_gain = out["pooled"]["accuracy"] - out["pooled"]["null_accuracy"]
    gw = out["group_weighted"]
    group_gain = gw["mean_accuracy"] - gw["mean_null_accuracy"]
    assert pooled_gain > 0.2          # pooled looks great...
    assert group_gain < 0.1           # ...within-group it is modest
    assert group_gain < pooled_gain / 2
