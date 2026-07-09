#!/usr/bin/env python
"""Head-to-head: AlphaMissense vs ESM masked-marginals on humsavar pathogenic-vs-benign.

For one protein, join (a) humsavar P/B labels, (b) ESM2 masked-marginal scores, (c) AlphaMissense
am_pathogenicity, on the shared variant set, and report AUROC + Spearman for BOTH predictors. This is
the "does the supervised predictor beat the zero-shot LM?" comparison — the standard variant-effect
question — on a FREE independent label set.

Predictor orientation: pathogenic => LOWER ESM masked-marginal score (log P(mut) - log P(wt)) and
HIGHER AlphaMissense am_pathogenicity. AUROC is computed with pathogenicity oriented so higher = more
pathogenic for both.

Run:
  python scripts/humsavar_am_vs_esm.py --uniprot P40692 --fasta mlh1.fasta \
      --humsavar humsavar.txt --am-tsv mlh1_am.tsv --model facebook/esm2_t12_35M_UR50D
"""
from __future__ import annotations

import argparse
import re
import sys

AA = list("ACDEFGHIKLMNPQRSTVWY")
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


def auroc(labels, scores):
    """Rank-based AUROC (higher score => class 1); average ranks for ties."""
    pairs = sorted(zip(scores, labels))
    n = len(pairs)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and pairs[j][0] == pairs[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0  # average of ranks i+1..j
        for k in range(i, j):
            ranks[k] = avg
        i = j
    npos = sum(labels)
    nneg = len(labels) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    sum_ranks_pos = sum(r for r, (_, l) in zip(ranks, pairs) if l == 1)
    u = sum_ranks_pos - npos * (npos + 1) / 2.0
    return u / (npos * nneg)


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j < len(v) and v[order[j]] == v[order[i]]:
                j += 1
            avg = (i + j - 1) / 2.0
            for k in range(i, j):
                r[order[k]] = avg
            i = j
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def load_humsavar(path, uniprot):
    labels = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if " VAR_" not in line:
            continue
        p = line.split()
        if p[1] != uniprot:
            continue
        ft = [x for x in p if x.startswith("VAR_")]
        i = p.index(ft[0])
        m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", p[i + 1])
        cat = p[i + 2]
        if not m or cat not in ("LP/P", "LB/B"):
            continue
        wt, pos, mut = AA3.get(m.group(1)), int(m.group(2)), AA3.get(m.group(3))
        if wt and mut and wt != mut:
            labels[f"{wt}{pos}{mut}"] = 1 if cat == "LP/P" else 0
    return labels


def esm_scores(seq, variants, model_name, maxlen=1022):
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name).eval()
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}
    positions = sorted({int(re.match(r"[A-Z](\d+)[A-Z]", v).group(1)) for v in variants})
    positions = [p for p in positions if p <= maxlen]
    enc = tok(seq[:maxlen], return_tensors="pt")
    ids = enc["input_ids"]
    logp_at = {}
    with torch.no_grad():
        for pos in positions:
            masked = ids.clone()
            masked[0, pos] = mask_id  # +1 offset handled: token index = pos (CLS at 0)
            out = model(masked).logits[0, pos].log_softmax(-1)
            logp_at[pos] = {aa: out[i].item() for aa, i in aa_ids.items()}
    scores = {}
    for v in variants:
        m = re.match(r"([A-Z])(\d+)([A-Z])", v)
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        if pos in logp_at and seq[pos - 1] == wt:
            scores[v] = logp_at[pos][mut] - logp_at[pos][wt]
    return scores


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--humsavar", required=True)
    ap.add_argument("--am-tsv", required=True, help="filtered AlphaMissense rows: variant\\tam_path\\tam_class")
    ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--maxlen", type=int, default=1022)
    args = ap.parse_args(argv)

    seq = "".join(l.strip() for l in open(args.fasta) if not l.startswith(">"))
    labels = load_humsavar(args.humsavar, args.uniprot)
    am = {}
    for line in open(args.am_tsv):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            try:
                am[parts[0]] = float(parts[1])
            except ValueError:
                pass
    esm = esm_scores(seq, list(labels), args.model, args.maxlen)

    shared = sorted(set(labels) & set(esm) & set(am))
    if len(shared) < 10:
        print(f"too few shared variants ({len(shared)}) — humsavar={len(labels)} esm={len(esm)} am={len(am)}")
        return 1
    y = [labels[v] for v in shared]
    esm_path = [-esm[v] for v in shared]   # higher => more pathogenic
    am_path = [am[v] for v in shared]
    npos, nneg = sum(y), len(y) - sum(y)
    print(f"{args.uniprot}  shared variants={len(shared)} (P={npos} B={nneg})  model={args.model}")
    print(f"  ESM masked-marginal : AUROC={auroc(y, esm_path):.3f}  |rho|={abs(spearman(esm_path, y)):.3f}")
    print(f"  AlphaMissense       : AUROC={auroc(y, am_path):.3f}  |rho|={abs(spearman(am_path, y)):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
