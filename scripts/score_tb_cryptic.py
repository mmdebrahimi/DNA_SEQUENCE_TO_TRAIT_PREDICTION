"""Step 5 — v1b lineage-collapsed TB scoring orchestrator + cohort/callability gate (deliverable a).

Pure scoring (`score_cohort`) takes {preds, labels, clusters} and reuses the FROZEN
`clonality.cluster_weighted_confusion` (brainstorm C2a — NOT representative-dedup). It emits the
lineage-collapsed sens/spec + raw sens/spec + raw->lineage shrinkage + n_discordant +
n_clusters_mixed_prediction (M1) + Wilson CI + effective_lineage_n + n_uncallable, tagged
`WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

Status gate (C3 + honesty):
  - lineage assignment unavailable (no clusters / all UNASSIGNED-degenerate) ->
    `LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL` (never a raw-only headline).
  - scored set is a convenience subset (not the full prevalence-preserving per-drug cohort) ->
    `TB_SUBSET_PLUMBING` (metrics computed but NEVER the baseline label).
  - else -> `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

The actual cohort fetch (full per-drug masked + regeno VCFs, ~1.6 TB -> D:) is the CLI runtime path
(`main`); it is NOT exercised by unit tests, which drive the pure scoring with synthetic inputs.
"""
from __future__ import annotations

from dna_decode.eval.clonality import (
    cluster_members,
    cluster_weighted_confusion,
    effective_lineage_n,
    wilson_ci,
)

BASELINE_LABEL = "WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE"
PLUMBING_LABEL = "TB_SUBSET_PLUMBING"
BLOCKED_LABEL = "LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL"

_R, _S, _ABSTAIN = "R", "S", "ABSTAIN"


def n_clusters_mixed_prediction(preds: dict[str, str], clusters: dict[str, int]) -> int:
    """M1: same-cluster member predictions disagree (within-lineage determinant heterogeneity)."""
    n = 0
    for sids in cluster_members(clusters).values():
        calls = {preds.get(s, _ABSTAIN) for s in sids} & {_R, _S}
        if calls == {_R, _S}:
            n += 1
    return n


def raw_confusion(preds: dict[str, str], labels: dict[str, str]) -> dict:
    """Per-isolate confusion (ABSTAIN excluded) — the un-collapsed, clonality-inflated baseline."""
    tp = fp = tn = fn = 0
    for sid, lab in labels.items():
        p = preds.get(sid, _ABSTAIN)
        L = str(lab).upper()
        if p == _R and L == _R:
            tp += 1
        elif p == _R and L == _S:
            fp += 1
        elif p == _S and L == _S:
            tn += 1
        elif p == _S and L == _R:
            fn += 1
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def _lineage_available(clusters: dict[str, int]) -> bool:
    # degenerate iff every isolate is its own singleton (no real collapse happened)
    if not clusters:
        return False
    return len(set(clusters.values())) < len(clusters)


def score_cohort(
    preds: dict[str, str],
    labels: dict[str, str],
    clusters: dict[str, int],
    *,
    drug: str,
    cohort_complete: bool,
) -> dict:
    """Lineage-collapsed scoring + status gate. Returns a self-describing result dict."""
    n_abstain = sum(1 for s in labels if preds.get(s, _ABSTAIN) == _ABSTAIN)
    base = {
        "drug": drug,
        "rule_status": "deterministic",
        "n_isolates": len(labels),
        "n_uncallable_abstain": n_abstain,
    }

    if not _lineage_available(clusters):
        return {**base, "status": BLOCKED_LABEL,
                "reason": "no non-singleton lineage clusters — lineage assignment unavailable"}

    cw = cluster_weighted_confusion(preds, labels, clusters)
    raw = raw_confusion(preds, labels)
    eff_r = effective_lineage_n(clusters, labels, "R")
    eff_s = effective_lineage_n(clusters, labels, "S")
    raw_r = sum(1 for s, l in labels.items() if str(l).upper() == _R)
    raw_s = sum(1 for s, l in labels.items() if str(l).upper() == _S)

    result = {
        **base,
        "status": BASELINE_LABEL if cohort_complete else PLUMBING_LABEL,
        "honesty": ("In-distribution knowledge-baseline: the WHO catalogue was built partly from "
                    "CRyPTIC. NOT independent validation — see the separate post-2023 gold-set arm."),
        "lineage_collapsed": {
            "sens": cw["sens"], "spec": cw["spec"],
            "tp": cw["tp"], "fp": cw["fp"], "tn": cw["tn"], "fn": cw["fn"],
            "sens_wilson_ci": wilson_ci(cw["tp"], cw["tp"] + cw["fn"]),
            "spec_wilson_ci": wilson_ci(cw["tn"], cw["tn"] + cw["fp"]),
            "n_clusters_R": cw["n_clusters_R"], "n_clusters_S": cw["n_clusters_S"],
            "n_discordant": cw["n_discordant"], "n_cluster_abstain": cw["n_cluster_abstain"],
            "n_clusters_mixed_prediction": n_clusters_mixed_prediction(preds, clusters),
        },
        "raw": raw,
        "effective_lineage_n": {"R": eff_r, "S": eff_s},
        "raw_to_lineage_shrinkage": {"R": [raw_r, eff_r], "S": [raw_s, eff_s]},
    }
    if not cohort_complete:
        result["plumbing_note"] = ("convenience/partial cohort — metrics are NOT the baseline; "
                                   "fetch the full per-drug prevalence-preserving cohort to earn the label")
    return result
