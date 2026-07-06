#!/usr/bin/env python
"""esm_zeroshot_dms.py - J2: real ESM-2 protein model vs wet-lab DMS.

Runs the real learned protein model (ESM-2, zero-shot masked marginals) to
score protein variants, correlates those scores against wet-lab DMS
measurements (ProteinGym), and reports the median rank correlation across
assays plus a shuffled control.

Target to beat: median |Spearman| ~= 0.48 (published ESM2-650M ProteinGym
number). Pass bar: median |Spearman| >= 0.45 AND shuffled control ~0.

Inputs:
  <data-dir>/pg_reference.csv
  <data-dir>/pg_dms/DMS_ProteinGym_substitutions/*.csv
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
    """Spearman rho via rank + Pearson, without scipy."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 5:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else float("nan")


def parse_variant(mutant):
    """Parse M1A -> (M, 1, A). Position is 1-indexed."""
    return mutant[0], int(mutant[1:-1]), mutant[-1]


def load_reference(ref_csv):
    refs = {}
    with open(ref_csv, encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            refs[row["DMS_id"]] = (row["DMS_filename"], row["target_seq"])
    return refs


def load_dms(path):
    out = {}
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            mutant = row.get("mutant", "")
            if not mutant or ":" in mutant or ";" in mutant:
                continue
            try:
                out[mutant] = float(row["DMS_score"])
            except (ValueError, KeyError, TypeError):
                continue
    return out


def load_transformers_runner(model_name, auto_tokenizer, auto_model, device):
    """Load via transformers/Hugging Face."""
    tok = auto_tokenizer.from_pretrained(model_name)
    try:
        model = auto_model.from_pretrained(model_name)
        return {
            "kind": "transformers",
            "device": device,
            "model": model.to(device).eval(),
            "tokenizer": tok,
            "mask_id": tok.mask_token_id,
            "aa_ids": {aa: tok.convert_tokens_to_ids(aa) for aa in AA},
        }
    except OSError as exc:
        message = str(exc)
        if "403" not in message and "Forbidden" not in message:
            raise
        print("NOTE: safetensors fetch hit 403; retrying with pytorch_model.bin.")
        model = auto_model.from_pretrained(model_name, use_safetensors=False)
        return {
            "kind": "transformers",
            "device": device,
            "model": model.to(device).eval(),
            "tokenizer": tok,
            "mask_id": tok.mask_token_id,
            "aa_ids": {aa: tok.convert_tokens_to_ids(aa) for aa in AA},
        }


def load_fair_esm_runner(model_name, device):
    """Load via fair-esm checkpoints when Hugging Face payloads are blocked."""
    if "TORCH_HOME" not in os.environ and "HF_HOME" in os.environ:
        os.environ["TORCH_HOME"] = os.path.join(os.environ["HF_HOME"], "torch_home")

    try:
        import esm
    except ImportError as exc:
        raise ImportError(
            "fair-esm is required for the direct-checkpoint fallback. Install it with `pip install fair-esm`."
        ) from exc

    bare_model_name = model_name.split("/")[-1]
    loader = getattr(esm.pretrained, bare_model_name, None)
    if loader is None:
        raise ValueError(f"No fair-esm loader found for model {bare_model_name}")

    model, alphabet = loader()
    return {
        "kind": "fair-esm",
        "device": device,
        "model": model.to(device).eval(),
        "batch_converter": alphabet.get_batch_converter(),
        "mask_id": alphabet.mask_idx,
        "aa_ids": {aa: alphabet.get_idx(aa) for aa in AA},
    }


def load_runner(model_name, auto_tokenizer, auto_model, device):
    """Prefer transformers; fall back to fair-esm when weight downloads are blocked."""
    try:
        return load_transformers_runner(model_name, auto_tokenizer, auto_model, device)
    except OSError as exc:
        message = str(exc)
        if "403" not in message and "Forbidden" not in message:
            raise
        print("NOTE: Hugging Face weights remain blocked; switching to fair-esm checkpoint download.")
        return load_fair_esm_runner(model_name, device)


def score_assay(seq, dms, runner, batch, maxlen, torch):
    """ESM-2 masked-marginals score = logP(mut|context) - logP(wt|context)."""
    if len(seq) > maxlen:
        return None, "too_long"

    if runner["kind"] == "transformers":
        ids = runner["tokenizer"](seq, return_tensors="pt")["input_ids"].to(runner["device"])
    else:
        _, _, ids = runner["batch_converter"]([("seq", seq)])
        ids = ids.to(runner["device"])

    length = len(seq)
    variants = []
    positions = set()
    mismatches = 0

    for mutant, score in dms.items():
        try:
            wt, pos, mut = parse_variant(mutant)
        except (ValueError, IndexError):
            continue
        if pos < 1 or pos > length or wt not in AA or mut not in AA:
            continue
        if seq[pos - 1] != wt:
            mismatches += 1
            continue
        variants.append((wt, pos, mut, score))
        positions.add(pos)

    if len(variants) < 20:
        return None, "too_few"

    positions = sorted(positions)
    logp_at = {}
    with torch.no_grad():
        for start in range(0, len(positions), batch):
            chunk = positions[start:start + batch]
            stack = ids.repeat(len(chunk), 1)
            for row_index, pos in enumerate(chunk):
                stack[row_index, pos] = runner["mask_id"]
            if runner["kind"] == "transformers":
                logits = runner["model"](stack).logits
            else:
                logits = runner["model"](stack, repr_layers=[], return_contacts=False)["logits"]
            for row_index, pos in enumerate(chunk):
                logp_at[pos] = torch.log_softmax(logits[row_index, pos], dim=-1).float().cpu()

    xs, ys = [], []
    for wt, pos, mut, score in variants:
        lp = logp_at[pos]
        xs.append(float(lp[runner["aa_ids"][mut]] - lp[runner["aa_ids"][wt]]))
        ys.append(score)

    rho = spearman(xs, ys)
    rng = np.random.default_rng(0)
    shuffled = np.array(ys)
    rng.shuffle(shuffled)

    return {
        "n": len(xs),
        "rho": rho,
        "rho_shuf": spearman(xs, shuffled),
        "mism": mismatches,
    }, None


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    parser.add_argument("--max-assays", type=int, default=40)
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--maxlen", type=int, default=1022)
    args = parser.parse_args(argv)

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("=== ESM-2 zero-shot vs wet-lab DMS ===")
    print(f"model={args.model}  device={device}  batch={args.batch}  maxlen={args.maxlen}")
    if device == "cpu":
        print("WARNING: no CUDA - 650M on CPU is slow. Use a GPU box or smaller model.")

    runner = load_runner(args.model, AutoTokenizer, AutoModelForMaskedLM, device)

    refs = load_reference(os.path.join(args.data_dir, "pg_reference.csv"))
    dms_dir = os.path.join(args.data_dir, "pg_dms", "DMS_ProteinGym_substitutions")

    results = []
    skipped_long = 0
    skipped_few = 0
    for dms_id, (filename, seq) in refs.items():
        if len(results) >= args.max_assays:
            break
        path = os.path.join(dms_dir, filename)
        if not os.path.exists(path) or not seq:
            continue
        dms = load_dms(path)
        if not dms:
            continue
        result, why = score_assay(
            seq,
            dms,
            runner,
            args.batch,
            args.maxlen,
            torch,
        )
        if why == "too_long":
            skipped_long += 1
            continue
        if why == "too_few":
            skipped_few += 1
            continue
        if result:
            result["dms_id"] = dms_id
            results.append(result)
            print(
                f"  {dms_id[:44]:44s} n={result['n']:5d}  rho={result['rho']:+.3f}  "
                f"shuffled={result['rho_shuf']:+.3f}"
            )

    if not results:
        print("FAIL: no assays scored (check --data-dir paths).")
        return 1

    med = float(np.median([abs(r["rho"]) for r in results]))
    med_shuf = float(np.median([abs(r["rho_shuf"]) for r in results]))
    print(
        f"\nESM-2 vs DMS over {len(results)} assays: median |Spearman| = {med:.3f}  "
        f"(shuffled control = {med_shuf:.3f})"
    )
    print(f"skipped: {skipped_long} too-long (>{args.maxlen} aa), {skipped_few} too-few-variants")

    ok = med >= PASS_BAR and med_shuf < SHUFFLE_MAX
    print()
    if ok:
        tag = "  (>= 0.48 stretch - matches the field!)" if med >= STRETCH else ""
        print(f"PASS: median |rho| {med:.3f} >= {PASS_BAR}{tag}. Real learned protein model captures the signal.")
    else:
        print(f"FAIL: median |rho| {med:.3f} < {PASS_BAR} or shuffled {med_shuf:.3f} >= {SHUFFLE_MAX}.")
        print("      Long proteins are skipped, not windowed - known follow-up if the run underperforms.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
