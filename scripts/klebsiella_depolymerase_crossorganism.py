"""CROSS-ORGANISM test: does the RBP/depolymerase -> phenotype paradigm transfer to Klebsiella capsule?

The E. coli phage receptor cell decodes phage genome/RBP -> outer-membrane receptor. This asks the
breadth question: does the SAME deterministic sequence-homology-transfer idea predict a DIFFERENT host's
(Klebsiella pneumoniae) DIFFERENT phenotype (capsule KL-type) from the phage DEPOLYMERASE (the enzymatic
capsule-degrading tail-spike domain)? Depolymerases are MORE modular than tail fibers, so a homology caller
may transfer BETTER on capsule than the E. coli cross-lab RBP number (0.364).

Data: DpoTropiSearch (Concha-Eloko/Nat Comms 2025, "Unlocking data in Klebsiella lysogens...";
github.com/conchaeloko/DpoTropiSearch + Zenodo 10.5281/zenodo.14065540) — INDEPENDENT of LBNL/BASEL, a
DIFFERENT LAB + DIFFERENT ORGANISM. `KL_type_LCA` = the capsule type (prophage-host LCA-inferred);
`domain_seq` = the depolymerase enzymatic domain.

Method: exactly the phage RBP caller architecture (`protein_kmers` + Jaccard nearest-neighbour), retargeted:
domain k-mer nearest-neighbour transfer -> predicted KL-type, leave-one-out. Reported vs a prior-frequency
null. HONEST scope: labels are prophage-inferred (in-distribution, like the within-LBNL LOO), so this is the
paradigm-transfer analogue of that number, NOT an independent wet-lab number (the 63 exp_validated set is
the gold-standard follow-on).

Usage:  uv run python scripts/klebsiella_depolymerase_crossorganism.py \
          --labels D:/dna_decode_cache/kleb/dpo_labels.tsv --min-members 5 --cap 20
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from dna_decode.phage.rbp_caller import kmer_similarity, protein_kmers


def load(labels_tsv: str, min_members: int, cap: int, seq_col: str):
    rows = []
    with open(labels_tsv, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kl = (r.get("KL_type_LCA") or "").strip()
            seq = (r.get(seq_col) or "").strip()
            if not seq or not kl or ";" in kl or "," in kl or kl.upper() in ("NAN", "NA", "NONE", ""):
                continue
            rows.append((r.get("Protein_name") or "", kl, seq))
    by_kl = defaultdict(list)
    for pid, kl, seq in rows:
        by_kl[kl].append((pid, seq))
    # KL-types with enough members; cap per class (stratified subsample; DETERMINISTIC — first-cap)
    sub = []
    for kl, members in by_kl.items():
        if len(members) >= min_members:
            sub.extend((pid, kl, seq) for pid, seq in members[:cap])
    return sub


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default="D:/dna_decode_cache/kleb/dpo_labels.tsv")
    ap.add_argument("--seq-col", default="domain_seq", choices=["domain_seq", "seq"])
    ap.add_argument("--min-members", type=int, default=5)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--min-similarity", type=float, default=0.05)
    ap.add_argument("--date", default="2026-07-25")
    args = ap.parse_args()

    sub = load(args.labels, args.min_members, args.cap, args.seq_col)
    labels = [kl for _, kl, _ in sub]
    kmers = [protein_kmers(seq, k=args.k) for _, _, seq in sub]
    n = len(sub)
    n_kl = len(set(labels))
    prior = Counter(labels)
    # prior-frequency null: always guess the most common KL-type
    null_acc = prior.most_common(1)[0][1] / n

    correct = called = 0
    per_kl = defaultdict(lambda: [0, 0])
    for i in range(n):
        best_j, best_sim = -1, 0.0
        ki = kmers[i]
        for j in range(n):
            if j == i:
                continue
            s = kmer_similarity(ki, kmers[j])
            if s > best_sim:
                best_sim, best_j = s, j
        if best_j < 0 or best_sim < args.min_similarity:
            continue  # INDETERMINATE abstain
        called += 1
        pred = labels[best_j]
        ok = pred == labels[i]
        per_kl[labels[i]][1] += 1
        if ok:
            correct += 1
            per_kl[labels[i]][0] += 1

    acc = correct / called if called else None
    # well-powered KL-types (>=5 called)
    strong = {k: v for k, v in per_kl.items() if v[1] >= 5}
    out = {
        "cell": "klebsiella_depolymerase_KLtype_crossorganism", "date": args.date,
        "question": "does the RBP/depolymerase->phenotype paradigm transfer cross-organism to Klebsiella capsule KL-type?",
        "source": "DpoTropiSearch (Concha-Eloko/Nat Comms 2025; Zenodo 10.5281/zenodo.14065540) — different lab + different organism",
        "seq_col": args.seq_col, "k": args.k,
        "n_depolymerases": n, "n_KL_types": n_kl,
        "subsample": f"KL-types with >={args.min_members} members, capped {args.cap}/type (stratified)",
        "called": called, "correct": correct,
        "crossorganism_LOO_accuracy": acc,
        "prior_null_accuracy": null_acc,
        "lift_over_null": (acc - null_acc) if acc is not None else None,
        "per_KLtype_correct_called_wellpowered": {k: v for k, v in sorted(strong.items(), key=lambda x: -x[1][1])[:20]},
        "honest_scope": "labels are prophage-host-LCA-INFERRED (in-distribution, like the within-LBNL LOO), NOT "
                        "independent wet-lab. This is the paradigm-TRANSFER analogue. The 63 experimentally-"
                        "validated depolymerases are the independent gold-standard follow-on. domain_seq = the "
                        "modular enzymatic capsule-degrading domain (the cleaner sequence->function unit).",
    }
    Path("wiki").mkdir(exist_ok=True)
    (Path("wiki") / f"klebsiella_crossorganism_result_{args.date}.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_depolymerases", "n_KL_types", "called", "correct",
          "crossorganism_LOO_accuracy", "prior_null_accuracy", "lift_over_null")}, indent=2))


if __name__ == "__main__":
    main()
