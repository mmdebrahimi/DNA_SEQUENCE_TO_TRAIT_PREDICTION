"""The open falsifier: does the router's 211/211 survive a LINEAGE-disjoint hold-out?

    uv run python scripts/identify_lineage_disjoint.py

WHY THIS EXISTS. The shipped leave-one-out is 211/211, but it holds out a GENOME, and 146 of 211
(69.2%) have a ~99.5%-ANI-or-closer twin still in the reference -- so for most folds a near-clone
stands in for the held-out genome and the fold is nearly free. The contract names this as the one open
falsifier. This runs it.

METHOD. Cluster the reference by Mash distance with the FROZEN chaining-resistant
`clonality.greedy_representative_clusters_from_matrix` (NOT single-linkage, which chains A~B~C into one
blob), then hold out WHOLE CLUSTERS: when scoring a genome, every member of its own cluster is removed
from the reference, not just itself. A genome whose entire lineage is held out must be identified from
a genuinely different lineage of the same organism, or abstain.

NO NEW CONTAINER WORK. Everything comes from the cached all-pairs matrix.

PRE-REGISTERED READING, fixed before the run so the result cannot be reinterpreted afterwards:
  * accuracy-on-called DROPS a lot and abstention RISES  -> the 211/211 was carried by near-twins, as
    suspected. That is the honest expectation and NOT a defect: abstaining when no same-lineage
    reference exists is the safe failure this router is built for.
  * accuracy holds                                       -> the router generalises across lineages of
    a supported organism, which is a STRONGER claim than anything currently in the contract.
  * genomes become unscorable because holding out their cluster empties their organism -> report them
    separately as UNSCORABLE, never as misses (the K. oxytoca n=1 precedent).
"""
from __future__ import annotations

import argparse
import collections
import json
from datetime import date
from pathlib import Path

import numpy as np

from dna_decode.eval.clonality import DistanceMatrix, greedy_representative_clusters_from_matrix
from dna_decode.identify.core import accession_from_reference_id, decide, parse_mash_dist
from dna_decode.identify.thresholds import frozen

REPO = Path(__file__).resolve().parent.parent
MATRIX_CACHE = Path("D:/dna_decode_cache/identify/allpairs.tsv")

#: Lineage thresholds to sweep. 0.005 ~ 99.5% ANI (the twin bar the clonality figure uses); the
#: coarser rungs merge progressively more of each organism into one lineage.
THRESHOLDS = (0.001, 0.005, 0.01, 0.02)


def load_matrix(text: str, labels: dict[str, str]) -> tuple[DistanceMatrix, dict[str, dict[str, float]]]:
    """Build a DistanceMatrix over labelled accessions + a nested dist lookup."""
    d: dict[str, dict[str, float]] = collections.defaultdict(dict)
    for h in parse_mash_dist(text):
        a = accession_from_reference_id(h.reference_id)
        b = accession_from_reference_id(h.query_id)
        if a in labels and b in labels:
            d[b][a] = h.distance
    ids = sorted(labels)
    m = np.zeros((len(ids), len(ids)), dtype=float)
    idx = {a: i for i, a in enumerate(ids)}
    for b, row in d.items():
        for a, v in row.items():
            m[idx[b], idx[a]] = v
    return DistanceMatrix(strain_ids=ids, matrix=m), d


def score_lineage_disjoint(dists, labels, clusters, thresholds) -> dict:
    """Hold out each genome's WHOLE cluster, then identify it."""
    per_org: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    confusions: collections.Counter = collections.Counter()
    unscorable: list[str] = []

    members = collections.defaultdict(set)
    for acc, cid in clusters.items():
        members[cid].add(acc)

    for acc, truth in labels.items():
        held = members[clusters[acc]]
        # does any genome of this organism survive outside the held-out cluster?
        survivors = [a for a, o in labels.items() if o == truth and a not in held]
        if not survivors:
            unscorable.append(acc)
            per_org[truth]["unscorable_organism_emptied"] += 1
            continue
        hits = sorted(
            (h for a, dv in [(a, dists[acc].get(a)) for a in labels if a not in held]
             if dv is not None
             for h in [type("H", (), {"reference_id": f"/x/{a}/genome.fna", "query_id": acc,
                                      "distance": dv, "p_value": 0.0, "shared_hashes": ""})()]),
            key=lambda h: h.distance)
        call = decide(hits, thresholds, lambda a: labels.get(a), lambda o: o)
        if call.abstained:
            per_org[truth]["abstain"] += 1
            per_org[truth][f"abstain_{call.reason.value}"] += 1
        elif call.organism == truth:
            per_org[truth]["hit"] += 1
        else:
            per_org[truth]["miss"] += 1
            confusions[f"{truth}->{call.organism}"] += 1

    hit = sum(c["hit"] for c in per_org.values())
    miss = sum(c["miss"] for c in per_org.values())
    ab = sum(c["abstain"] for c in per_org.values())
    n = hit + miss + ab
    return {
        "n_scored": n, "hit": hit, "miss": miss, "abstain": ab,
        "accuracy_on_called": round(hit / (hit + miss), 4) if (hit + miss) else None,
        "abstention_rate": round(ab / n, 4) if n else None,
        "per_organism": {o: dict(c) for o, c in sorted(per_org.items())},
        "confusions": dict(confusions.most_common()),
        "n_unscorable_organism_emptied": len(unscorable),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=None)
    a = ap.parse_args(argv)

    manifest = a.manifest or sorted((REPO / "wiki").glob("identify_reference_manifest_*.json"))[-1]
    labels = {acc: o for o, accs in
              json.loads(manifest.read_text(encoding="utf-8"))["accessions_by_organism"].items()
              for acc in accs}
    if not MATRIX_CACHE.exists():
        print(f"missing {MATRIX_CACHE}; run scripts/identify_threshold_analysis.py first")
        return 2

    dm, dists = load_matrix(MATRIX_CACHE.read_text(encoding="utf-8", errors="replace"), labels)
    th = frozen()
    out: dict = {"schema": "identify-lineage-disjoint-v1", "date": date.today().isoformat(),
                 "manifest": manifest.name,
                 "thresholds": {"max_distance": th.max_distance,
                                "ambiguity_margin": th.ambiguity_margin},
                 "method": "whole-CLUSTER hold-out using the frozen chaining-resistant "
                           "greedy-representative clustering; every member of a genome's own lineage "
                           "is removed from the reference before identifying it",
                 "baseline_genome_level_loo": {"n_scored": 211, "hit": 211,
                                               "accuracy_on_called": 1.0,
                                               "source": "wiki/identify_validation_2026-10-04.json"},
                 "sweep": {}}

    print(f"{'thr':>7} {'clusters':>9} {'largest':>8} {'n':>5} {'hit':>5} {'miss':>5} "
          f"{'abst':>5} {'acc_called':>11} {'abst_rate':>10}")
    for t in THRESHOLDS:
        clusters = greedy_representative_clusters_from_matrix(dm, t)
        sizes = collections.Counter(clusters.values())
        r = score_lineage_disjoint(dists, labels, clusters, th)
        r["n_clusters"] = len(sizes)
        r["largest_cluster_fraction"] = round(max(sizes.values()) / len(clusters), 4)
        out["sweep"][str(t)] = r
        print(f"{t:>7} {len(sizes):>9} {r['largest_cluster_fraction']:>8} {r['n_scored']:>5} "
              f"{r['hit']:>5} {r['miss']:>5} {r['abstain']:>5} "
              f"{str(r['accuracy_on_called']):>11} {str(r['abstention_rate']):>10}")

    dest = REPO / "wiki" / f"identify_lineage_disjoint_{date.today().isoformat()}.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nartifact: {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
