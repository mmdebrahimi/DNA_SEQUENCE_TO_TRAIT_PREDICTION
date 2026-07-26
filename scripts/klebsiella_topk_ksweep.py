"""Klebsiella cross-organism transfer — top-K accuracy (promiscuity metric) + k-mer-size sweep.

Extends the top-1 cross-organism number (0.453 clonality-corrected; wiki/klebsiella_crossorganism_result)
with the BIOLOGICALLY-CORRECT metric: depolymerases are promiscuous (one enzyme degrades several capsule
types; DpoTropiSearch reports top-10 hits), so a phage-therapy match wants a RANKED list of KL-types, not a
single guess. This computes clonality-corrected top-1/top-3/top-5 accuracy (the true KL-type is among the K
nearest neighbours' types) vs a top-K prior-frequency null, and sweeps k in {3,4,5,6} to find the best
depolymerase-domain resolution.

Data: DpoTropiSearch (Concha-Eloko et al., Nat Commun 2025; Zenodo 10.5281/zenodo.14065540). The data is
NON-COMMERCIAL-licensed + lives on D: (gitignored, NEVER redistributed) — this script is my own code doing
non-commercial research USE, and commits only the numeric result (no sequences).

Usage:  uv run python scripts/klebsiella_topk_ksweep.py --labels D:/dna_decode_cache/kleb/dpo_labels.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from dna_decode.phage.rbp_caller import kmer_similarity, protein_kmers


def load_stratified(labels_tsv: str, min_members: int, cap: int):
    rows = []
    with open(labels_tsv, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kl = (r.get("KL_type_LCA") or "").strip()
            seq = (r.get("domain_seq") or "").strip()
            if not seq or not kl or ";" in kl or "," in kl or kl.upper() in ("NAN", "NA", "NONE", ""):
                continue
            rows.append((kl, seq))
    by_kl = defaultdict(list)
    for kl, seq in rows:
        by_kl[kl].append(seq)
    sub = []
    for kl, seqs in by_kl.items():
        if len(seqs) >= min_members:
            sub.extend((kl, s) for s in seqs[:cap])
    return sub


def greedy_collapse(sub, k: int, thresh: float = 0.90):
    """Chaining-resistant greedy-representative clustering on domain k-mers (clonality correction)."""
    km = [protein_kmers(s, k) for _, s in sub]
    reps = []
    for i in range(len(sub)):
        if not any(kmer_similarity(km[i], km[r]) >= thresh for r in reps):
            reps.append(i)
    return [sub[i] for i in reps]


def topk_loo(reps, k: int, topk_vals=(1, 3, 5), min_similarity: float = 0.05):
    """For each rep, rank all others by domain k-mer similarity; the true KL-type is a HIT at level T if it
    appears among the KL-types of the T nearest CALLED neighbours."""
    labels = [kl for kl, _ in reps]
    km = [protein_kmers(s, k) for _, s in reps]
    n = len(reps)
    called = 0
    hits = {t: 0 for t in topk_vals}
    for i in range(n):
        sims = []
        for j in range(n):
            if j == i:
                continue
            s = kmer_similarity(km[i], km[j])
            if s >= min_similarity:
                sims.append((s, labels[j]))
        if not sims:
            continue
        called += 1
        sims.sort(key=lambda x: x[0], reverse=True)
        # top-T predictions = the KL-types of the T nearest neighbours, dedup preserving rank
        ranked_types = []
        for _, lab in sims:
            if lab not in ranked_types:
                ranked_types.append(lab)
        for t in topk_vals:
            if labels[i] in ranked_types[:t]:
                hits[t] += 1
    # top-K prior null: the T most-frequent KL-types
    prior = Counter(labels)
    most = [kl for kl, _ in prior.most_common()]
    null = {}
    for t in topk_vals:
        top_set = set(most[:t])
        null[t] = sum(1 for lab in labels if lab in top_set) / n
    return {"n_reps": n, "called": called,
            "topk_accuracy": {t: (hits[t] / called if called else None) for t in topk_vals},
            "topk_null": null,
            "topk_lift": {t: ((hits[t] / called) - null[t] if called else None) for t in topk_vals}}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default="D:/dna_decode_cache/kleb/dpo_labels.tsv")
    ap.add_argument("--min-members", type=int, default=5)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--date", default="2026-07-25")
    args = ap.parse_args()

    sub = load_stratified(args.labels, args.min_members, args.cap)
    result = {"cell": "klebsiella_depolymerase_topk_ksweep", "date": args.date,
              "source": "DpoTropiSearch (Concha-Eloko/Nat Commun 2025; Zenodo 10.5281/zenodo.14065540) — non-commercial-licensed data, D:-only, not redistributed",
              "subsample": f"KL-types >={args.min_members} members, cap {args.cap}/type",
              "n_subsample": len(sub), "by_k": {}}
    print(f"subsample n={len(sub)}", flush=True)
    for k in (3, 4, 5, 6):
        reps = greedy_collapse(sub, k=k, thresh=0.90)   # collapse at each k (k-specific k-mer space)
        r = topk_loo(reps, k=k)
        result["by_k"][str(k)] = r
        acc = r["topk_accuracy"]
        print(f"k={k}: reps={r['n_reps']} called={r['called']} | "
              f"top1={acc[1]:.3f} top3={acc[3]:.3f} top5={acc[5]:.3f} "
              f"(null1={r['topk_null'][1]:.3f})", flush=True)

    # headline = best-k top-K (clonality-corrected)
    best_k = max(result["by_k"], key=lambda kk: result["by_k"][kk]["topk_accuracy"][1])
    result["headline"] = {"best_k": int(best_k),
                          "top1": result["by_k"][best_k]["topk_accuracy"][1],
                          "top3": result["by_k"][best_k]["topk_accuracy"][3],
                          "top5": result["by_k"][best_k]["topk_accuracy"][5]}
    Path("wiki").mkdir(exist_ok=True)
    (Path("wiki") / f"klebsiella_topk_ksweep_{args.date}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("DONE ->", f"wiki/klebsiella_topk_ksweep_{args.date}.json", flush=True)


if __name__ == "__main__":
    main()
