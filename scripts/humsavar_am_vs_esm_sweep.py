#!/usr/bin/env python
"""Multi-protein AlphaMissense-vs-ESM head-to-head on humsavar (robustness across proteins).

Loops several proteins through the single-protein harness (humsavar_am_vs_esm) and reports the
per-protein AUROC for both predictors + the median across proteins. n=1 is anecdote; a median over
~7 proteins is the robust statement of "supervised beats zero-shot on clinical pathogenicity".

Usage: python scripts/humsavar_am_vs_esm_sweep.py --datasrc <dir> [--model ...] [--out <json>]
Expects in <dir>: humsavar.txt, <AC>.fasta per protein, mlh1_am.tsv (variant\\tam) + sweep_am.tsv
(uniprot\\tvariant\\tam) — as staged by the sweep fetch.
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from humsavar_am_vs_esm import auroc, esm_scores, load_humsavar, spearman  # noqa: E402

# (uniprot, fasta_basename, am_file_basename, am_file_has_uniprot_column)
PROTEINS = [
    ("P40692", "mlh1.fasta", "mlh1_am.tsv", False),
    ("P43246", "P43246.fasta", "sweep_am.tsv", True),
    ("P04070", "P04070.fasta", "sweep_am.tsv", True),
    ("P02730", "P02730.fasta", "sweep_am.tsv", True),
    ("P51608", "P51608.fasta", "sweep_am.tsv", True),
    ("Q9BXM7", "Q9BXM7.fasta", "sweep_am.tsv", True),
    ("Q99972", "Q99972.fasta", "sweep_am.tsv", True),
]


def load_am(path, uniprot, has_uni):
    am = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        try:
            if has_uni and len(p) >= 3 and p[0] == uniprot:
                am[p[1]] = float(p[2])
            elif not has_uni and len(p) >= 2:
                am[p[0]] = float(p[1])
        except ValueError:
            pass
    return am


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasrc", required=True)
    ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    d = args.datasrc
    hum = os.path.join(d, "humsavar.txt")

    rows = []
    for uni, fa, amf, hu in PROTEINS:
        seq = "".join(l.strip() for l in open(os.path.join(d, fa)) if not l.startswith(">"))
        lab = load_humsavar(hum, uni)
        am = load_am(os.path.join(d, amf), uni, hu)
        esm = esm_scores(seq, list(lab), args.model)
        sh = sorted(set(lab) & set(esm) & set(am))
        if len(sh) < 10:
            print(f"{uni}: too few shared ({len(sh)}) — skip")
            continue
        y = [lab[v] for v in sh]
        ae = auroc(y, [-esm[v] for v in sh])
        aa = auroc(y, [am[v] for v in sh])
        rows.append({"uniprot": uni, "n": len(sh), "n_path": sum(y), "n_benign": len(y) - sum(y),
                     "esm_auroc": round(ae, 3), "am_auroc": round(aa, 3), "gap": round(aa - ae, 3)})
        print(f"{uni}  n={len(sh):3d} (P={sum(y):3d} B={len(y)-sum(y):3d})  "
              f"ESM_AUROC={ae:.3f}  AM_AUROC={aa:.3f}  gap={aa-ae:+.3f}")

    esm_med = statistics.median([r["esm_auroc"] for r in rows])
    am_med = statistics.median([r["am_auroc"] for r in rows])
    am_wins = sum(1 for r in rows if r["gap"] > 0)
    print(f"\nMEDIAN over {len(rows)} proteins: ESM={esm_med:.3f}  AM={am_med:.3f}  "
          f"gap={am_med-esm_med:+.3f}  (AM wins {am_wins}/{len(rows)})")
    if args.out:
        json.dump({"model": args.model, "n_proteins": len(rows), "esm_median_auroc": esm_med,
                   "am_median_auroc": am_med, "median_gap": round(am_med - esm_med, 3),
                   "am_wins": am_wins, "per_protein": rows,
                   "task": "humsavar pathogenic-vs-benign", "date": "2026-07-08"},
                  open(args.out, "w"), indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
