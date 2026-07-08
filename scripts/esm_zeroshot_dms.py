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


def window_for_position(pos, L, maxlen):
    """Long-protein windowing (Phase 2): for a protein longer than `maxlen`, return the maxlen-residue
    window CENTERED on 1-based residue `pos`, clamped to sequence ends. Returns
    (start0, end0, local_token_index): 0-based residue slice [start0:end0] + the 1-based token index of
    `pos` inside that window (CLS at token 0, first window residue at token 1). Pure position math — the
    bug-prone part, so it is unit-tested independently of any model."""
    if L <= maxlen:
        return 0, L, pos                      # whole-sequence path: token index == residue index
    half = maxlen // 2
    start0 = max(0, min(pos - 1 - half, L - maxlen))
    end0 = start0 + maxlen
    return start0, end0, pos - start0          # 1-based local token index of `pos`


def combine_variant_scores(per_model_scores):
    """Ensemble combine (Phase 2): given a list of {variant_key: score} dicts (one per model), z-score each
    model's vector over the variants ALL models share, then average. Returns {variant_key: mean_z}. Rank-
    based ensembling via z-scores is the classic ProteinGym top approach — free (just more inference). Pure
    + unit-tested; a constant-score model contributes 0 (its z-vector is all-zero)."""
    if not per_model_scores:
        return {}
    keys = set(per_model_scores[0])
    for d in per_model_scores[1:]:
        keys &= set(d)
    keys = sorted(keys)
    if not keys:
        return {}
    combined = {k: 0.0 for k in keys}
    n = len(per_model_scores)
    for d in per_model_scores:
        v = np.array([d[k] for k in keys], dtype=float)
        sd = v.std()
        z = (v - v.mean()) / sd if sd > 0 else np.zeros_like(v)
        for k, zz in zip(keys, z):
            combined[k] += float(zz) / n
    return combined


def self_test_report(refs, maxlen, exists_fn):
    """No-GPU pre-run guard: given the parsed reference {id:(dms_csv_path, target_seq)}, report how many
    assays have their per-assay CSV attached + how many exceed `maxlen` (recovered only with --long-mode
    window). Lets the user confirm the ProteinGym dataset is correctly attached BEFORE spending GPU minutes.
    Pure (exists_fn injected) → unit-tested."""
    total = len(refs)
    have_csv = long_n = scorable_skip = 0
    for _id, (path, seq) in refs.items():
        present = bool(exists_fn(path))
        L = len(seq)
        if present:
            have_csv += 1
            if L <= maxlen:
                scorable_skip += 1
        if L > maxlen:
            long_n += 1
    return {"total": total, "have_csv": have_csv, "long_gt_maxlen": long_n,
            "scorable_skip_mode": scorable_skip, "scorable_window_mode": have_csv}


def score_assay(seq, dms, tok, model, device, mask_id, aa_ids, batch, maxlen, torch,
                long_mode="skip", keep_scores=False):
    """ESM-2 masked-marginals: score = logP(mut|context) - logP(wt|context), summed per variant.

    long_mode: 'skip' (baseline — proteins > maxlen dropped, the published-0.48 behavior) or 'window'
    (Phase 2 — each mutated position scored in a maxlen window centered on it, recovering long assays).
    keep_scores: also return per-variant [key, x, y] rows (for the ensemble-merge path)."""
    L = len(seq)
    if L > maxlen and long_mode != "window":
        return None, "too_long"

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
        if L <= maxlen:
            enc = tok(seq, return_tensors="pt")
            ids = enc["input_ids"].to(device)          # [1, L+2]  (CLS at 0, residues 1..L, EOS at L+1)
            for i in range(0, len(positions), batch):
                chunk = positions[i:i + batch]
                stack = ids.repeat(len(chunk), 1)          # [b, L+2]
                for r, p in enumerate(chunk):
                    stack[r, p] = mask_id                  # token index p == 1-indexed residue p
                logits = model(stack).logits               # [b, L+2, V]
                for r, p in enumerate(chunk):
                    logp_at[p] = torch.log_softmax(logits[r, p], dim=-1).float().cpu()
        else:                                              # windowed long-protein path (Phase 2)
            for i in range(0, len(positions), batch):
                chunk = positions[i:i + batch]
                wins = [window_for_position(p, L, maxlen) for p in chunk]
                seqs = [seq[s:e] for (s, e, _) in wins]     # each exactly maxlen residues
                locs = [loc for (_, _, loc) in wins]
                enc = tok(seqs, return_tensors="pt", padding=True)
                stack = enc["input_ids"].to(device)
                for r, loc in enumerate(locs):
                    stack[r, loc] = mask_id
                logits = model(stack).logits
                for r, (p, loc) in enumerate(zip(chunk, locs)):
                    logp_at[p] = torch.log_softmax(logits[r, loc], dim=-1).float().cpu()

    xs, ys, var = [], [], []
    for wt, pos, mut, y in variants:
        lp = logp_at[pos]
        x = float(lp[aa_ids[mut]] - lp[aa_ids[wt]])
        xs.append(x); ys.append(y)
        if keep_scores:
            var.append([f"{wt}{pos}{mut}", x, y])

    rho = spearman(xs, ys)
    rng = np.random.default_rng(0)
    yss = np.array(ys); rng.shuffle(yss)
    out = {"n": len(xs), "rho": rho, "rho_shuf": spearman(xs, yss), "mism": mism}
    if keep_scores:
        out["var"] = var
    return out, None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True, help="dir with pg_reference.csv + pg_dms/DMS_ProteinGym_substitutions/")
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--max-assays", type=int, default=40)
    ap.add_argument("--batch", type=int, default=16, help="masked positions per forward pass")
    ap.add_argument("--maxlen", type=int, default=1022, help="ESM context cap (residues); longer proteins skipped")
    ap.add_argument("--dtype", choices=["float32", "float16"], default="float32",
                    help="float16 fits ESM2-3B on a free T4/P100 (Phase 2 — the certain lift over 650M)")
    ap.add_argument("--long-mode", choices=["skip", "window"], default="skip",
                    help="skip = baseline (drop >maxlen proteins); window = Phase 2 (centered-window scoring)")
    args = ap.parse_args(argv)

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"=== ESM-2 zero-shot vs wet-lab DMS ===")
    print(f"model={args.model}  device={device}  batch={args.batch}  maxlen={args.maxlen}  "
          f"dtype={args.dtype}  long_mode={args.long_mode}")
    if device == "cpu":
        print("WARNING: no CUDA — 650M on CPU is slow. Use a smaller --model or a GPU box.")

    tok = AutoTokenizer.from_pretrained(args.model)
    dt = torch.float16 if args.dtype == "float16" else torch.float32
    model = AutoModelForMaskedLM.from_pretrained(args.model, torch_dtype=dt).to(device).eval()
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
        r, why = score_assay(seq, dms, tok, model, device, mask_id, aa_ids, args.batch, args.maxlen, torch,
                             long_mode=args.long_mode)
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
