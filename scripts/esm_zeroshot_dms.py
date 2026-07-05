#!/usr/bin/env python
"""esm_zeroshot_dms.py - J2: real ESM-2 protein model vs wet-lab DMS (SELF-CONTAINED).

Runs the REAL learned protein model (ESM-2, zero-shot masked-marginals) to score protein
variants, correlates its scores against wet-lab deep-mutational-scan (DMS) measurements
(ProteinGym), reports the median rank-correlation across assays + a shuffled control.

Target to beat: median |Spearman| ~= 0.48 (the published ESM2-650M ProteinGym number).
PASS bar (pre-registered): median |Spearman| >= 0.45 AND shuffled control ~0.

No AlphaMissense, no network at RUN time (weights download once on first load, then cache).
GPU strongly recommended for the 650M model. Weights are FREE (MIT license).

Inputs (a --data-dir holding the ProteinGym substitution benchmark):
  <data-dir>/pg_reference.csv                         (has DMS_id, DMS_filename, target_seq)
  <data-dir>/pg_dms/DMS_ProteinGym_substitutions/*.csv (each: mutant, DMS_score columns)

Run:
  python esm_zeroshot_dms.py --data-dir <DIR> --model facebook/esm2_t33_650M_UR50D --max-assays 40
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

AA = "ACDEFGHIKLMNPQRSTVWY"
PASS_BAR = 0.45
STRETCH = 0.48
SHUFFLE_MAX = 0.05


def spearman(x, y):
    """Spearman rho via rank + Pearson (no scipy)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 5:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


def parse_variant(m):
    """'M1A' -> ('M', 1, 'A')  (position is 1-indexed)."""
    return m[0], int(m[1:-1]), m[-1]


def load_reference(ref_csv):
    refs = {}
    with open(ref_csv, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            refs[row["DMS_id"]] = (row["DMS_filename"], row["target_seq"])
    return refs


def load_dms(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            m = row.get("mutant", "")
            if not m or ":" in m or ";" in m:   # single mutants only
                continue
            try:
                out[m] = float(row["DMS_score"])
            except (ValueError, KeyError, TypeError):
                continue
    return out


def score_assay(seq, dms, tok, model, device, mask_id, aa_ids, batch, maxlen, torch):
    """ESM-2 masked-marginals: score = logP(mut|context) - logP(wt|context), summed per variant."""
    if len(seq) > maxlen:
        return None, "too_long"
    enc = tok(seq, return_tensors="pt")
    ids = enc["input_ids"].to(device)          # [1, L+2]  (CLS at 0, residues 1..L, EOS at L+1)
    L = len(seq)

    variants, positions, mism = [], set(), 0
    for m, y in dms.items():
        try:
            wt, pos, mut = parse_variant(m)
        except (ValueError, IndexError):
            continue
        if pos < 1 or pos > L or wt not in AA or mut not in AA:
            continue
        if seq[pos - 1] != wt:                 # reference-offset mismatch -> skip honestly
            mism += 1
            continue
        variants.append((wt, pos, mut, y))
        positions.add(pos)
    if len(variants) < 20:
        return None, "too_few"

    positions = sorted(positions)
    logp_at = {}
    with torch.no_grad():
        for i in range(0, len(positions), batch):
            chunk = positions[i:i + batch]
            stack = ids.repeat(len(chunk), 1)          # [b, L+2]
            for r, p in enumerate(chunk):
                stack[r, p] = mask_id                  # token index p == 1-indexed residue p
            logits = model(stack).logits               # [b, L+2, V]
            for r, p in enumerate(chunk):
                logp_at[p] = torch.log_softmax(logits[r, p], dim=-1).float().cpu()

    xs, ys = [], []
    for wt, pos, mut, y in variants:
        lp = logp_at[pos]
        xs.append(float(lp[aa_ids[mut]] - lp[aa_ids[wt]]))
        ys.append(y)

    rho = spearman(xs, ys)
    rng = np.random.default_rng(0)
    yss = np.array(ys); rng.shuffle(yss)
    return {"n": len(xs), "rho": rho, "rho_shuf": spearman(xs, yss), "mism": mism}, None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True, help="dir with pg_reference.csv + pg_dms/DMS_ProteinGym_substitutions/")
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--max-assays", type=int, default=40)
    ap.add_argument("--batch", type=int, default=16, help="masked positions per forward pass")
    ap.add_argument("--maxlen", type=int, default=1022, help="ESM context cap (residues); longer proteins skipped")
    args = ap.parse_args(argv)

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"=== ESM-2 zero-shot vs wet-lab DMS ===")
    print(f"model={args.model}  device={device}  batch={args.batch}  maxlen={args.maxlen}")
    if device == "cpu":
        print("WARNING: no CUDA — 650M on CPU is slow. Use a smaller --model or a GPU box.")

    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model).to(device).eval()
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}

    refs = load_reference(os.path.join(args.data_dir, "pg_reference.csv"))
    dms_dir = os.path.join(args.data_dir, "pg_dms", "DMS_ProteinGym_substitutions")

    results, skipped_long, skipped_few = [], 0, 0
    for dms_id, (fn, seq) in refs.items():
        if len(results) >= args.max_assays:
            break
        path = os.path.join(dms_dir, fn)
        if not os.path.exists(path) or not seq:
            continue
        dms = load_dms(path)
        if not dms:
            continue
        r, why = score_assay(seq, dms, tok, model, device, mask_id, aa_ids, args.batch, args.maxlen, torch)
        if why == "too_long":
            skipped_long += 1; continue
        if why == "too_few":
            skipped_few += 1; continue
        if r:
            r["dms_id"] = dms_id
            results.append(r)
            print(f"  {dms_id[:44]:44s} n={r['n']:5d}  rho={r['rho']:+.3f}  shuffled={r['rho_shuf']:+.3f}")

    if not results:
        print("FAIL: no assays scored (check --data-dir paths).")
        return 1

    med = float(np.median([abs(r["rho"]) for r in results]))
    med_shuf = float(np.median([abs(r["rho_shuf"]) for r in results]))
    print(f"\nESM-2 vs DMS over {len(results)} assays: median |Spearman| = {med:.3f}  "
          f"(shuffled control = {med_shuf:.3f})")
    print(f"skipped: {skipped_long} too-long (>{args.maxlen} aa, ESM cap), {skipped_few} too-few-variants")

    ok = med >= PASS_BAR and med_shuf < SHUFFLE_MAX
    print()
    if ok:
        tag = "  (>= 0.48 stretch — matches the field!)" if med >= STRETCH else ""
        print(f"PASS: median |rho| {med:.3f} >= {PASS_BAR}{tag}. The REAL learned protein model captures")
        print("      the causal molecular signal -> J2 done; JEPA/CLIP protein direction is green-lit on our own model.")
    else:
        print(f"FAIL: median |rho| {med:.3f} < {PASS_BAR} or shuffled {med_shuf:.3f} >= {SHUFFLE_MAX}.")
        print("      (Long proteins are skipped, not windowed — a known follow-up that lifts the score.)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
