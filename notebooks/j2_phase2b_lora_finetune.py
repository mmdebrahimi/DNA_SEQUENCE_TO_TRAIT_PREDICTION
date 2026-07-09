#!/usr/bin/env python
"""J2 Phase 2b — LoRA fine-tune ESM-2 to BEAT its zero-shot number on HELD-OUT proteins (leakage-safe).

WHAT: the genuine-training step of the world-model bet. Phase 2 (fp16/3B/ensemble) beats 0.48 by scaling +
combining FROZEN open weights (inference only). Phase 2b asks the deeper question: does LoRA-adapting the
ESM-2 backbone on some proteins improve GENERALIZABLE variant-effect prediction on proteins it never saw?

HONEST BAR (cross-protein, leakage-controlled — the real world-model-improvement test, NOT within-assay CV):
  * Split PROTEINS into folds so NO protein's assays appear in both train and test (split_proteins, below).
  * LoRA-adapt ESM-2 on the TRAIN proteins' variants (pairwise ranking loss vs the wet-lab DMS score).
  * Evaluate masked-marginals on the HELD-OUT proteins with BOTH the base model AND the adapted model.
  * HEADLINE = finetuned_median - zeroshot_median on the SAME held-out fold. >0 (with a real margin) means
    fine-tuning taught the world model transferable biology; <=0 means it only memorized train proteins.
  Comparing to the global 0.48 would be dishonest here (the fold is a subset) — the base-vs-adapted delta on
  the identical held-out set is the fair comparator.

WHY within-assay CV is NOT the bar: training a head on 80% of an assay + predicting its own held-out 20%
mostly measures the readout regressor, not the backbone's knowledge. Cross-protein is the world-model claim.

SCAFFOLD SCOPE: the leakage-safe split + the pair sampler (the bug-prone pure logic) are unit-tested on CPU.
The LoRA training + eval loop is real but GPU-run (peft + torch lazy-imported) — upload this file alongside
j2_phase1_esm2_proteingym.py to a free Kaggle/Colab kernel. It reuses the Phase-1 drift-guarded scoring core
(spearman / score_assay) so eval is byte-identical to the zero-shot path.

Run (free T4/P100; needs the ProteinGym data attached — see notebooks/J2_PHASE1_RUNBOOK.md Step 0):
  pip install peft
  python j2_phase2b_lora_finetune.py --data-dir /kaggle/input/<slug> --nfolds 5 --fold 0 \
      --model facebook/esm2_t33_650M_UR50D --dtype float16 --rank 8 --epochs 2 --out j2b_fold0.json
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys

import numpy as np

# Reuse the Phase-1 drift-guarded pure core (upload both files together; same dir on Kaggle/Colab).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from j2_phase1_esm2_proteingym import (  # noqa: E402
    AA, PASS_BAR, STRETCH, _summary, load_dms, load_reference, parse_variant, score_assay, spearman,
)


# ---- leakage-safe protein split (pure, unit-tested — the load-bearing correctness piece) ----
def protein_key(dms_id, uniprot=""):
    """Group assays of the SAME protein so no protein leaks across the train/test boundary. Prefer the
    reference's UniProt_ID; fall back to the DMS_id prefix before the first '_' (the UniProt-ish token,
    e.g. 'A0A140D2T1_ZIKV_Sourisseau_2019' -> 'A0A140D2T1')."""
    return (uniprot or "").strip() or dms_id.split("_")[0]


def split_proteins(assay_keys, nfolds, seed):
    """assay_keys: {dms_id: protein_key}. Returns {dms_id: fold_index} with ALL assays of a protein in ONE
    fold (leakage-safe by construction). Deterministic given seed (seeded shuffle → round-robin)."""
    if nfolds < 2:
        raise ValueError(f"nfolds must be >= 2, got {nfolds}")
    proteins = sorted(set(assay_keys.values()))
    rng = random.Random(seed)
    rng.shuffle(proteins)
    fold_of = {p: i % nfolds for i, p in enumerate(proteins)}
    return {a: fold_of[k] for a, k in assay_keys.items()}


def sample_pairs(n, max_pairs, seed):
    """Deterministic index pairs (i, j), i<j, for the pairwise ranking loss over one assay's n variants.
    Returns up to max_pairs unique pairs (all pairs if the total <= max_pairs). Pure + unit-tested."""
    if n < 2:
        return []
    allp = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(allp) <= max_pairs:
        return allp
    rng = random.Random(seed)
    rng.shuffle(allp)
    return sorted(allp[:max_pairs])


def assay_keys_from_reference(refs, uniprot_by_id=None):
    """{dms_id: protein_key} from a load_reference() result. uniprot_by_id (optional) maps dms_id->UniProt."""
    uniprot_by_id = uniprot_by_id or {}
    return {dms_id: protein_key(dms_id, uniprot_by_id.get(dms_id, "")) for dms_id in refs}


def load_uniprot_map(data_dir):
    """{DMS_id: UniProt_ID} from the reference CSV, when the column exists (official ProteinGym has it)."""
    import csv
    for name in ("DMS_substitutions.csv", "pg_reference.csv"):
        p = os.path.join(data_dir, name)
        if os.path.exists(p):
            out = {}
            with open(p, encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    if row.get("DMS_id"):
                        out[row["DMS_id"]] = row.get("UniProt_ID", "")
            return out
    return {}


# ---- variant extraction (front-half of score_assay, reused for building training examples) ----
def assay_variants(seq, dms, maxlen):
    """[(wt, pos, mut, y)] for single mutants that match the reference sequence, within maxlen."""
    if len(seq) > maxlen:
        return []
    out = []
    for m, y in dms.items():
        try:
            wt, pos, mut = parse_variant(m)
        except (ValueError, IndexError):
            continue
        if pos < 1 or pos > len(seq) or wt not in AA or mut not in AA or seq[pos - 1] != wt:
            continue
        out.append((wt, pos, mut, y))
    return out


# ---- LoRA training (GPU-run; peft + torch lazy-imported) ----
def _assay_ranking_loss(model, tok, seq, variants, device, torch, mask_id, aa_ids, batch, max_pairs, seed):
    """Differentiable masked-marginal scores for one assay's variants + a pairwise margin ranking loss so
    predicted order matches the wet-lab DMS order. Reuses the Phase-1 scoring definition (logP(mut)-logP(wt))
    but keeps the graph for backprop."""
    positions = sorted({p for _w, p, _m, _y in variants})
    enc = tok(seq, return_tensors="pt")
    ids = enc["input_ids"].to(device)
    logp_at = {}
    for i in range(0, len(positions), batch):
        chunk = positions[i:i + batch]
        stack = ids.repeat(len(chunk), 1)
        for r, p in enumerate(chunk):
            stack[r, p] = mask_id
        logits = model(stack).logits
        for r, p in enumerate(chunk):
            logp_at[p] = torch.log_softmax(logits[r, p].float(), dim=-1)
    scores = torch.stack([logp_at[p][aa_ids[mut]] - logp_at[p][aa_ids[wt]]
                          for wt, p, mut, _y in variants])
    ys = [y for _w, _p, _m, y in variants]
    pairs = sample_pairs(len(variants), max_pairs, seed)
    if not pairs:
        return scores.sum() * 0.0
    losses = []
    for i, j in pairs:
        if ys[i] == ys[j]:
            continue
        sign = 1.0 if ys[i] > ys[j] else -1.0        # variant i should score higher iff its DMS is higher
        losses.append(torch.clamp(0.5 - sign * (scores[i] - scores[j]), min=0.0))   # margin=0.5
    if not losses:
        return scores.sum() * 0.0
    return torch.stack(losses).mean()


def train_lora(train_examples, base_model, rank, epochs, lr, device, torch, tok, mask_id, aa_ids,
               batch, max_pairs, seed):
    """LoRA-adapt ESM-2 on TRAIN proteins' variants. Returns the adapted model (adapters active).

    TRAINING DTYPE = float32 ALWAYS (not the eval --dtype): raw-fp16 LoRA training with AdamW NaNs/stalls
    — small (lr*grad) updates underflow in fp16, and there is no autocast/GradScaler here. ESM2-650M in
    fp32 fits a free T4 16 GB comfortably (~2.6 GB weights). For a 3B backbone, fp32 training is tight on a
    free T4 — prefer bf16 + autocast (a documented follow-up), NOT raw fp16. `tok` is passed explicitly
    (no module global)."""
    from peft import LoraConfig, get_peft_model
    from transformers import AutoModelForMaskedLM
    model = AutoModelForMaskedLM.from_pretrained(base_model, torch_dtype=torch.float32)
    cfg = LoraConfig(r=rank, lora_alpha=2 * rank, lora_dropout=0.05, bias="none",
                     target_modules=["query", "key", "value"])   # ESM-2 attention projections
    model = get_peft_model(model, cfg).to(device)
    model.train()
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=lr)
    for ep in range(epochs):
        random.Random(seed + ep).shuffle(train_examples)
        tot = 0.0
        for k, (seq, variants) in enumerate(train_examples):
            loss = _assay_ranking_loss(model, tok, seq, variants, device, torch, mask_id, aa_ids,
                                       batch, max_pairs, seed + k)
            opt.zero_grad(); loss.backward(); opt.step()
            tot += float(loss.detach())
        print(f"  epoch {ep + 1}/{epochs}: mean assay ranking loss = {tot / max(1, len(train_examples)):.4f}")
    model.eval()
    return model


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/kaggle/input/proteingym")
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--nfolds", type=int, default=5)
    ap.add_argument("--fold", type=int, default=0, help="held-out fold index (0..nfolds-1)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--rank", type=int, default=8, help="LoRA rank")
    ap.add_argument("--epochs", type=int, default=2)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--dtype", choices=["float32", "float16"], default="float16")
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--maxlen", type=int, default=1022)
    ap.add_argument("--max-pairs", type=int, default=2000, help="pairwise-ranking pairs sampled per assay")
    ap.add_argument("--max-train-assays", type=int, default=0, help="0 = all train-fold assays")
    ap.add_argument("--out", default="", help="write per-fold result JSON here")
    ap.add_argument("--dry-run", action="store_true",
                    help="no GPU: just report the leakage-safe train/test split for this fold, then exit")
    args = ap.parse_args(argv)

    refs = load_reference(args.data_dir)
    uni = load_uniprot_map(args.data_dir)
    folds = split_proteins(assay_keys_from_reference(refs, uni), args.nfolds, args.seed)
    test_ids = sorted([d for d, f in folds.items() if f == args.fold])
    train_ids = sorted([d for d, f in folds.items() if f != args.fold])
    n_test_prot = len({protein_key(d, uni.get(d, "")) for d in test_ids})
    n_train_prot = len({protein_key(d, uni.get(d, "")) for d in train_ids})
    leak = {protein_key(d, uni.get(d, "")) for d in test_ids} & {protein_key(d, uni.get(d, "")) for d in train_ids}
    print(f"fold {args.fold}/{args.nfolds}: train {len(train_ids)} assays / {n_train_prot} proteins | "
          f"test {len(test_ids)} assays / {n_test_prot} proteins | protein-leak={len(leak)} (must be 0)")
    assert not leak, "LEAKAGE: a protein appears in both train and test — split is broken"

    if args.dry_run:
        print("DRY-RUN: leakage-safe split verified; no GPU used.")
        if args.out:
            json.dump({"fold": args.fold, "nfolds": args.nfolds, "n_train_assays": len(train_ids),
                       "n_test_assays": len(test_ids), "n_train_proteins": n_train_prot,
                       "n_test_proteins": n_test_prot, "protein_leak": len(leak)},
                      open(args.out, "w", encoding="utf-8"), indent=2)
        return 0

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = torch.float16 if args.dtype == "float16" else torch.float32   # EVAL/inference dtype only
    if device == "cpu":
        print("WARNING: no CUDA — LoRA fine-tune on CPU is impractically slow. Use a free T4/P100 kernel.")
    tok = AutoTokenizer.from_pretrained(args.model)
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}
    dms_of = {}
    for d in refs:
        path, seq = refs[d]
        if os.path.exists(path) and seq:
            dms_of[d] = (path, seq)

    # build train examples (leakage-safe: train fold only)
    train_examples = []
    for d in train_ids:
        if d not in dms_of:
            continue
        path, seq = dms_of[d]
        variants = assay_variants(seq, load_dms(path), args.maxlen)
        if len(variants) >= 20:
            train_examples.append((seq, variants))
        if args.max_train_assays and len(train_examples) >= args.max_train_assays:
            break
    print(f"training on {len(train_examples)} train-fold assays (>=20 single mutants each)")

    # eval helper: median |Spearman| over the HELD-OUT test fold with a given model
    def eval_fold(model):
        rhos = []
        for d in test_ids:
            if d not in dms_of:
                continue
            path, seq = dms_of[d]
            r, why = score_assay(seq, load_dms(path), tok, model, device, mask_id, aa_ids,
                                 args.batch, args.maxlen, torch, long_mode="skip")
            if r:
                rhos.append(r["rho"])
        return _summary([{"rho": x, "rho_shuf": 0.0} for x in rhos])

    base = AutoModelForMaskedLM.from_pretrained(args.model, torch_dtype=dtype).to(device).eval()
    zs = eval_fold(base)
    print(f"ZERO-SHOT on held-out fold: median |Spearman| = {zs['median_abs_rho']:.3f} "
          f"({zs['n_assays']} assays)")
    del base

    adapted = train_lora(train_examples, args.model, args.rank, args.epochs, args.lr, device,
                         torch, tok, mask_id, aa_ids, args.batch, args.max_pairs, args.seed)
    ft = eval_fold(adapted)
    delta = ft["median_abs_rho"] - zs["median_abs_rho"]
    print(f"FINE-TUNED on held-out fold: median |Spearman| = {ft['median_abs_rho']:.3f}")
    print(f"\nDELTA (finetuned - zeroshot) on the SAME held-out proteins = {delta:+.3f}  "
          f"-> {'IMPROVED the world model' if delta > 0.005 else 'no cross-protein gain'}")

    if args.out:
        json.dump({"fold": args.fold, "nfolds": args.nfolds, "model": args.model, "rank": args.rank,
                   "epochs": args.epochs, "n_train_assays": len(train_examples),
                   "n_test_assays": zs["n_assays"], "zeroshot_median": zs["median_abs_rho"],
                   "finetuned_median": ft["median_abs_rho"], "delta": delta},
                  open(args.out, "w", encoding="utf-8"), indent=2)
        print(f"wrote {args.out}")
    return 0 if delta > 0.005 else 1


if __name__ == "__main__":
    raise SystemExit(main())
