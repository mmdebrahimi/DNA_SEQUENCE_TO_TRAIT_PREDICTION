#!/usr/bin/env python
"""J2 Phase 1 — ESM-2 zero-shot vs wet-lab DMS, SELF-CONTAINED for a FREE cloud GPU (Kaggle / Colab).

WHAT: run the REAL learned protein model (ESM2-650M, zero-shot masked-marginals) over the joinable
ProteinGym deep-mutational-scan assays, correlate its per-variant scores against the WET-LAB measured
effects (Spearman), report the median across assays + a shuffled negative control + leaderboard context.
This is the "get the real ~0.48 field number" step (J2 Phase 1). It is INFERENCE ONLY — no training,
free MIT weights, ~3-4 GB VRAM → fits a free T4/P100 (16 GB) with huge headroom.

WHY SELF-CONTAINED: this ONE file is meant to be uploaded to a fresh Kaggle/Colab kernel and run with
nothing else — no repo checkout, no PAT. The scoring core (spearman / parse_variant / load_dms /
score_assay masked-marginals) is byte-faithful to the canonical `scripts/esm_zeroshot_dms.py`; a drift
guard (`tests/test_j2_phase1_notebook.py`) pins the two together so this copy can't silently diverge.

DISTRIBUTE FOR SPEED: `--shard i/n` runs only assay-subset i of n (strided, disjoint, covers-all). Run
shard 0/2 on Kaggle + shard 1/2 on Colab AT THE SAME TIME, each writes a JSON, then merge:
    python j2_phase1_esm2_proteingym.py --data-dir DIR --shard 0/2 --out shard0.json     # Kaggle
    python j2_phase1_esm2_proteingym.py --data-dir DIR --shard 1/2 --out shard1.json     # Colab
    python j2_phase1_esm2_proteingym.py --merge shard0.json shard1.json                  # anywhere (no GPU)

DATA (`--data-dir` must hold the ProteinGym substitution benchmark):
    <data-dir>/DMS_substitutions.csv   (or pg_reference.csv)  -> columns DMS_id, DMS_filename, target_seq
    <data-dir>/DMS_ProteinGym_substitutions/*.csv            (or pg_dms/DMS_ProteinGym_substitutions/*.csv)
  Get it by EITHER (see notebooks/J2_PHASE1_RUNBOOK.md for the full recipe):
   (a) Kaggle "Add Data": upload your local D:/dna_decode_cache/proteingym folder as a private dataset; OR
   (b) --fetch : download the reference file from the official ProteinGym GitHub (the per-assay CSVs still
       need attaching/downloading — the reference alone is not enough to score).

Target: median |Spearman| ~= 0.48 (published ESM2-650M ProteinGym number). PASS (pre-registered) >= 0.45.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import urllib.request

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

AA = "ACDEFGHIKLMNPQRSTVWY"
PASS_BAR = 0.45
STRETCH = 0.48
SHUFFLE_MAX = 0.05
# Official ProteinGym reference file (small; the per-assay DMS CSVs are large + attached/downloaded separately).
_REF_URL = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv"


# ---- pure scoring core (byte-faithful to scripts/esm_zeroshot_dms.py; drift-guarded by the test) ----
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
    """'M1A' -> ('M', 1, 'A')  (position 1-indexed)."""
    return m[0], int(m[1:-1]), m[-1]


def load_dms(path):
    """{mutant: DMS_score} for SINGLE mutants only."""
    out = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            m = row.get("mutant", "")
            if not m or ":" in m or ";" in m:
                continue
            try:
                out[m] = float(row["DMS_score"])
            except (ValueError, KeyError, TypeError):
                continue
    return out


# ---- data location (tolerant of both the official + the local cached layout) ----
def _find(data_dir, *candidates):
    for c in candidates:
        p = os.path.join(data_dir, c)
        if os.path.exists(p):
            return p
    return None


def load_reference(data_dir):
    """{DMS_id: (dms_csv_path, target_seq)} from DMS_substitutions.csv OR pg_reference.csv."""
    ref = _find(data_dir, "DMS_substitutions.csv", "pg_reference.csv")
    if not ref:
        raise SystemExit(f"no reference file (DMS_substitutions.csv / pg_reference.csv) under {data_dir} "
                         f"— see notebooks/J2_PHASE1_RUNBOOK.md, or pass --fetch to grab the reference.")
    dms_root = None
    for c in ("DMS_ProteinGym_substitutions", os.path.join("pg_dms", "DMS_ProteinGym_substitutions"),
              "DMS_substitutions", "."):
        if os.path.isdir(os.path.join(data_dir, c)):
            dms_root = os.path.join(data_dir, c); break
    if dms_root is None:
        dms_root = data_dir
    out = {}
    with open(ref, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            fn, seq = row.get("DMS_filename", ""), row.get("target_seq", "")
            if fn and seq:
                out[row["DMS_id"]] = (os.path.join(dms_root, fn), seq)
    return out


def fetch_reference(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    dst = os.path.join(data_dir, "DMS_substitutions.csv")
    print(f"downloading reference -> {dst}")
    urllib.request.urlretrieve(_REF_URL, dst)
    print("done. NOTE: the per-assay DMS CSVs are NOT fetched by this (large) — attach them via a Kaggle "
          "dataset or the ProteinGym Resources page (see the runbook).")


# ---- sharding (pure, testable) ----
def shard_assays(dms_ids, shard, nshards):
    """Deterministic strided partition: shard i of n = sorted(ids)[i::n]. Disjoint + covers-all."""
    if not (0 <= shard < nshards):
        raise ValueError(f"shard {shard} out of range for nshards {nshards}")
    return sorted(dms_ids)[shard::nshards]


def _summary(results):
    if not results:
        return {"n_assays": 0, "median_abs_rho": float("nan"), "median_abs_shuf": float("nan")}
    return {
        "n_assays": len(results),
        "median_abs_rho": float(np.median([abs(r["rho"]) for r in results])),
        "median_abs_shuf": float(np.median([abs(r["rho_shuf"]) for r in results])),
    }


def merge(paths):
    """Combine per-shard JSON result files -> dedup by dms_id -> pooled summary."""
    by_id = {}
    for p in paths:
        obj = json.load(open(p, encoding="utf-8"))
        for r in obj.get("results", []):
            by_id[r["dms_id"]] = r          # last writer wins; shards are disjoint so no real collision
    results = list(by_id.values())
    s = _summary(results)
    s["ok"] = bool(s["n_assays"] and s["median_abs_rho"] >= PASS_BAR and s["median_abs_shuf"] < SHUFFLE_MAX)
    return {"merged_from": list(paths), "results": results, **s}


# ---- ESM-2 masked-marginals (copied faithfully from scripts/esm_zeroshot_dms.py::score_assay) ----
def score_assay(seq, dms, tok, model, device, mask_id, aa_ids, batch, maxlen, torch):
    """ESM-2 masked-marginals: score = logP(mut|ctx) - logP(wt|ctx) per variant."""
    if len(seq) > maxlen:
        return None, "too_long"
    enc = tok(seq, return_tensors="pt")
    ids = enc["input_ids"].to(device)
    L = len(seq)
    variants, positions, mism = [], set(), 0
    for m, y in dms.items():
        try:
            wt, pos, mut = parse_variant(m)
        except (ValueError, IndexError):
            continue
        if pos < 1 or pos > L or wt not in AA or mut not in AA:
            continue
        if seq[pos - 1] != wt:
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
            stack = ids.repeat(len(chunk), 1)
            for r, p in enumerate(chunk):
                stack[r, p] = mask_id
            logits = model(stack).logits
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/kaggle/input/proteingym",
                    help="dir with DMS_substitutions.csv (or pg_reference.csv) + the per-assay DMS CSVs")
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--max-assays", type=int, default=0, help="0 = all joinable assays in this shard")
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--maxlen", type=int, default=1022)
    ap.add_argument("--shard", default="0/1", help="i/n — score only assay-subset i of n (distribute across kernels)")
    ap.add_argument("--out", default="", help="write per-assay results JSON here (for --merge later)")
    ap.add_argument("--fetch", action="store_true", help="download the reference file, then continue")
    ap.add_argument("--merge", nargs="+", default=None, help="MERGE these shard JSONs (no GPU) and exit")
    args = ap.parse_args(argv)

    if args.merge:
        m = merge(args.merge)
        print(json.dumps({k: v for k, v in m.items() if k != "results"}, indent=2))
        tag = "  (>= 0.48 — matches the field!)" if m["median_abs_rho"] >= STRETCH else ""
        print(f"\nMERGED: {m['n_assays']} assays, median |Spearman| = {m['median_abs_rho']:.3f}{tag} "
              f"(shuffled {m['median_abs_shuf']:.3f}) -> {'PASS' if m['ok'] else 'below 0.45'}")
        return 0 if m["ok"] else 1

    if args.fetch:
        fetch_reference(args.data_dir)

    shard, nshards = (int(x) for x in args.shard.split("/"))
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"=== J2 Phase 1 — ESM-2 zero-shot vs wet-lab DMS ===")
    print(f"model={args.model} device={device} shard={shard}/{nshards} batch={args.batch} maxlen={args.maxlen}")
    if device == "cpu":
        print("WARNING: no CUDA — 650M on CPU is slow. Use a free T4/P100 kernel or a smaller --model.")

    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model).to(device).eval()
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}

    refs = load_reference(args.data_dir)
    ids = shard_assays(list(refs), shard, nshards)
    print(f"joinable assays total={len(refs)} | this shard={len(ids)}\n")

    results, skipped_long, skipped_few, missing = [], 0, 0, 0
    for dms_id in ids:
        if args.max_assays and len(results) >= args.max_assays:
            break
        path, seq = refs[dms_id]
        if not os.path.exists(path) or not seq:
            missing += 1
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

    s = _summary(results)
    print(f"\nESM-2 vs DMS over {s['n_assays']} assays (shard {shard}/{nshards}): "
          f"median |Spearman| = {s['median_abs_rho']:.3f} (shuffled {s['median_abs_shuf']:.3f})")
    print(f"skipped: {skipped_long} too-long, {skipped_few} too-few, {missing} missing-file")

    if args.out:
        json.dump({"shard": args.shard, "model": args.model, "results": results, **s},
                  open(args.out, "w", encoding="utf-8"), indent=2)
        print(f"wrote {args.out}  (merge shards with:  --merge shard0.json shard1.json ...)")

    ok = bool(s["n_assays"] and s["median_abs_rho"] >= PASS_BAR and s["median_abs_shuf"] < SHUFFLE_MAX)
    if nshards > 1:
        print("\n(this is ONE shard — run the other shard(s) + `--merge` for the full-cohort number)")
    elif ok:
        tag = "  (>= 0.48 stretch — matches the field!)" if s["median_abs_rho"] >= STRETCH else ""
        print(f"\nPASS: median |rho| {s['median_abs_rho']:.3f} >= {PASS_BAR}{tag} — the REAL learned protein "
              f"model captures the causal molecular signal (J2 Phase 1 done).")
    else:
        print(f"\nbelow bar: median |rho| {s['median_abs_rho']:.3f} < {PASS_BAR} (or shuffled too high).")
    return 0 if (nshards > 1 or ok) else 1


if __name__ == "__main__":
    sys.exit(main())
