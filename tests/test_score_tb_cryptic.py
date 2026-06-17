"""Step 5 — v1b orchestrator: status gate, cluster_weighted_confusion reuse, M1, BLOCKED."""
from __future__ import annotations

from scripts import score_tb_cryptic as sc


def _cluster_two_lineages():
    # 4 isolates: lineage A (s1,s2 both R) + lineage B (s3,s4 both S)
    clusters = {"s1": 0, "s2": 0, "s3": 1, "s4": 1}
    labels = {"s1": "R", "s2": "R", "s3": "S", "s4": "S"}
    preds = {"s1": "R", "s2": "R", "s3": "S", "s4": "S"}
    return preds, labels, clusters


def test_full_cohort_earns_baseline_label_and_collapses():
    preds, labels, clusters = _cluster_two_lineages()
    out = sc.score_cohort(preds, labels, clusters, drug="rifampicin", cohort_complete=True)
    assert out["status"] == sc.BASELINE_LABEL
    # 4 isolates -> 2 lineage votes (1 R cluster, 1 S cluster)
    assert out["effective_lineage_n"] == {"R": 1, "S": 1}
    assert out["raw_to_lineage_shrinkage"] == {"R": [2, 1], "S": [2, 1]}
    assert out["lineage_collapsed"]["sens"] == 1.0 and out["lineage_collapsed"]["spec"] == 1.0
    assert "sens_wilson_ci" in out["lineage_collapsed"]


def test_partial_cohort_is_plumbing_not_baseline():
    preds, labels, clusters = _cluster_two_lineages()
    out = sc.score_cohort(preds, labels, clusters, drug="rifampicin", cohort_complete=False)
    assert out["status"] == sc.PLUMBING_LABEL
    assert "plumbing_note" in out


def test_all_singleton_clusters_is_blocked():
    # every isolate its own cluster -> no real collapse -> BLOCKED
    labels = {"s1": "R", "s2": "S"}
    preds = {"s1": "R", "s2": "S"}
    clusters = {"s1": 0, "s2": 1}
    out = sc.score_cohort(preds, labels, clusters, drug="rifampicin", cohort_complete=True)
    assert out["status"] == sc.BLOCKED_LABEL


def test_mixed_prediction_within_lineage_counted():
    # one lineage (cluster 0) with disagreeing member predictions
    clusters = {"s1": 0, "s2": 0, "s3": 1, "s4": 1}
    preds = {"s1": "R", "s2": "S", "s3": "S", "s4": "S"}
    assert sc.n_clusters_mixed_prediction(preds, clusters) == 1


def test_abstain_excluded_and_counted():
    clusters = {"s1": 0, "s2": 0, "s3": 1, "s4": 1}
    labels = {"s1": "R", "s2": "R", "s3": "S", "s4": "S"}
    preds = {"s1": "R", "s2": "R", "s3": "ABSTAIN", "s4": "S"}  # s3 uncallable
    out = sc.score_cohort(preds, labels, clusters, drug="isoniazid", cohort_complete=True)
    assert out["n_uncallable_abstain"] == 1
    # raw spec ignores the ABSTAIN: only s4 is a scored S -> tn=1
    assert out["raw"]["tn"] == 1


def test_raw_confusion_basic():
    preds = {"a": "R", "b": "S", "c": "R", "d": "S"}
    labels = {"a": "R", "b": "S", "c": "S", "d": "R"}
    raw = sc.raw_confusion(preds, labels)
    assert (raw["tp"], raw["fp"], raw["tn"], raw["fn"]) == (1, 1, 1, 1)
