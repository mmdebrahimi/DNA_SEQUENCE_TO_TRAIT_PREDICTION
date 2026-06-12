"""Clonality clustering + cluster-weighted metrics for the lineage-disclosure layer.

The provenance-disjoint report card reports RAW-ISOLATE sens/spec: every selected
isolate is one vote. But every SCORED R class is clonally dominated (<=17 effective
lineages at Mash 0.001; Klebsiella cipro R = 60 isolates but only ~3 lineages). The
raw number is therefore clonality-inflated — a single over-sampled clone can carry
the metric. This module computes the HONEST companion: lineage-effective N and
cluster-weighted sens/spec (one vote per lineage), each with a Wilson confidence
interval (the cluster-weighted N is tiny, so the CI is the point — never render a
weighted estimate without it).

Design decisions (from the plan's brainstorm C1-C3):
  - C2 GREEDY-REPRESENTATIVE dedup, NOT single-linkage. `phylogeny.cluster_by_ani`
    is union-find (single-linkage) and CHAINS: A~B, B~C with A!~C merges all three.
    Greedy-representative picks a representative, drops everything within `threshold`
    OF THE REPRESENTATIVE, repeats — membership is distance-to-representative only,
    so a chain cannot silently collapse distinct lineages.
  - C1 a cluster mixing R + S labels is DISCORDANT — it is NEVER majority-voted into
    one class (that would hide a real clone whose members disagree on phenotype).
    DISCORDANT clusters are excluded from sens/spec and surfaced as a count.
  - C3 cluster-weighted sens/spec MUST be rendered with a Wilson CI + effective-N.

Pure-logic (clustering math, aggregation, CI) is split from the Docker Mash call so
the math is unit-testable offline against synthetic distance matrices.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from dna_decode.eval.phylogeny import DistanceMatrix, compute_mash_distances


# 95% Wilson score interval.
_WILSON_Z = 1.959963984540054


def _normalize_label(label: object) -> str:
    """Coerce a label to canonical "R"/"S".

    Accepts the project's two conventions: ints (1=R, 0=S) and strings ("R"/"S",
    case-insensitive). Raises on anything else so a silently-mislabeled member
    cannot corrupt a cluster's class.
    """
    if isinstance(label, str):
        u = label.strip().upper()
        if u in ("R", "S"):
            return u
        raise ValueError(f"unrecognized string label {label!r} (expected 'R'/'S')")
    if label == 1:
        return "R"
    if label == 0:
        return "S"
    raise ValueError(f"unrecognized label {label!r} (expected 1/0 or 'R'/'S')")


def greedy_representative_clusters_from_matrix(
    distance_matrix: DistanceMatrix,
    threshold: float,
) -> dict[str, int]:
    """Greedy-representative clustering on a precomputed distance matrix.

    Algorithm (chaining-resistant):
      1. Sort strain ids (determinism — greedy picking is order-sensitive).
      2. Walk the sorted ids. The first unassigned id becomes a new cluster's
         REPRESENTATIVE (next cluster id).
      3. Assign every still-unassigned id whose distance TO THE REPRESENTATIVE is
         <= threshold to that cluster.
      4. Repeat until all ids are assigned.

    Returns strain_id -> cluster_id, ids 0..k-1 in representative-selection order.
    """
    sids = sorted(distance_matrix.strain_ids)
    if not sids:
        return {}

    idx = {sid: i for i, sid in enumerate(distance_matrix.strain_ids)}
    m = distance_matrix.matrix
    assigned: dict[str, int] = {}
    next_cluster = 0
    for rep in sids:
        if rep in assigned:
            continue
        cid = next_cluster
        next_cluster += 1
        ri = idx[rep]
        assigned[rep] = cid
        for other in sids:
            if other in assigned:
                continue
            # symmetric matrix; read [rep, other]
            if float(m[ri, idx[other]]) <= threshold:
                assigned[other] = cid
    return assigned


def greedy_representative_clusters(
    genomes: dict[str, Path],
    threshold: float,
    *,
    use_docker: bool = True,
) -> dict[str, int]:
    """Whole-cohort greedy-representative clustering from genome FASTAs.

    Runs ONE batched Mash sketch+dist over the whole cohort (via
    `phylogeny.compute_mash_distances`), then greedy-representative dedup at
    `threshold`. Whole-cohort (not pairwise-incremental) so the clustering sees
    every isolate at once.
    """
    dm = compute_mash_distances(genomes, use_docker=use_docker)
    return greedy_representative_clusters_from_matrix(dm, threshold)


def cluster_members(clusters: dict[str, int]) -> dict[int, list[str]]:
    """Invert strain->cluster into cluster_id -> sorted member strain_ids."""
    out: dict[int, list[str]] = {}
    for sid, cid in clusters.items():
        out.setdefault(cid, []).append(sid)
    for cid in out:
        out[cid].sort()
    return out


def cluster_class(member_labels) -> str:
    """Resolve a cluster's class from its members' true labels.

    Returns "R" (all members R), "S" (all members S), or "DISCORDANT" (members
    disagree). A discordant clone is NEVER majority-voted — it is surfaced as its
    own category (C1). Raises on an empty cluster.
    """
    norm = {_normalize_label(x) for x in member_labels}
    if not norm:
        raise ValueError("cannot classify an empty cluster")
    if norm == {"R"}:
        return "R"
    if norm == {"S"}:
        return "S"
    return "DISCORDANT"


def _cluster_prediction(member_preds) -> str:
    """Majority R/S prediction within a cluster; ABSTAIN if no R/S votes or a tie.

    A clone that collapses to one lineage-vote needs a single prediction. We take
    the majority of its members' R/S calls (ignoring ABSTAIN). An exact tie — or a
    cluster whose members all abstained — yields ABSTAIN (excluded from the
    confusion matrix) rather than an arbitrary call that could inflate the metric.
    """
    r = sum(1 for p in member_preds if str(p).upper() == "R")
    s = sum(1 for p in member_preds if str(p).upper() == "S")
    if r > s:
        return "R"
    if s > r:
        return "S"
    return "ABSTAIN"


def cluster_weighted_confusion(
    preds: dict[str, str],
    labels: dict[str, object],
    clusters: dict[str, int],
) -> dict:
    """Collapse each lineage to ONE vote, then compute confusion + sens/spec.

    Each cluster contributes a single (prediction, label) pair: the label is the
    cluster_class, the prediction is the within-cluster majority. DISCORDANT
    clusters are excluded from sens/spec and counted as `n_discordant`. Clusters
    that abstain (tie / all-abstain members) are excluded from the matrix and
    counted as `n_cluster_abstain`.

    Returns the same confusion shape as `independent_cohort_validate._conf`
    (tp/fp/tn/fn/sens/spec) plus lineage-level counts.
    """
    members = cluster_members(clusters)
    tp = fp = tn = fn = 0
    n_clusters_r = n_clusters_s = n_discordant = n_abstain = 0
    for cid, sids in members.items():
        member_labels = [labels[sid] for sid in sids]
        cls = cluster_class(member_labels)
        if cls == "DISCORDANT":
            n_discordant += 1
            continue
        if cls == "R":
            n_clusters_r += 1
        else:
            n_clusters_s += 1
        cpred = _cluster_prediction([preds.get(sid, "ABSTAIN") for sid in sids])
        if cpred == "ABSTAIN":
            n_abstain += 1
            continue
        if cpred == "R" and cls == "R":
            tp += 1
        elif cpred == "R" and cls == "S":
            fp += 1
        elif cpred == "S" and cls == "S":
            tn += 1
        elif cpred == "S" and cls == "R":
            fn += 1
    n = tp + fp + tn + fn
    return {
        "n_scored": n,
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "n_clusters_R": n_clusters_r,
        "n_clusters_S": n_clusters_s,
        "n_discordant": n_discordant,
        "n_cluster_abstain": n_abstain,
        "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "spec": round(tn / (tn + fp), 3) if (tn + fp) else None,
    }


def wilson_ci(k: int, n: int, z: float = _WILSON_Z) -> tuple[float, float]:
    """95% Wilson score interval for k successes in n trials.

    Returns (lo, hi), each clamped to [0, 1], rounded to 3 dp. n=0 -> (0.0, 1.0)
    (no information). The Wilson interval is well-behaved at the small n and
    extreme p (0 or 1) that the cluster-weighted metric produces, where the normal
    approximation would give a degenerate zero-width interval.
    """
    if n <= 0:
        return (0.0, 1.0)
    p = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2 * n)) / denom
    margin = (z / denom) * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))
    lo = max(0.0, center - margin)
    hi = min(1.0, center + margin)
    return (round(lo, 3), round(hi, 3))


def effective_lineage_n(
    clusters: dict[str, int],
    labels: dict[str, object],
    cls: str,
) -> int:
    """Count distinct same-label clusters of class `cls` ("R" or "S").

    This is the lineage-effective N for that class — the honest denominator behind
    the raw isolate count. DISCORDANT clusters count toward neither class.
    """
    want = _normalize_label(cls) if cls in ("R", "S", "r", "s", 0, 1) else cls
    members = cluster_members(clusters)
    return sum(
        1
        for sids in members.values()
        if cluster_class([labels[sid] for sid in sids]) == want
    )
