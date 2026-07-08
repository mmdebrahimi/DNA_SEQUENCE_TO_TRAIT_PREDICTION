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


def ensemble_merge(paths):
    """Ensemble across MODELS (Phase 2): each input JSON is a full run of a DIFFERENT model written with
    --keep-scores (per-assay `var` = [key, x, y] rows). For each assay present in ALL models, z-score-average
    the per-variant scores (combine_variant_scores) then recompute rho vs the wet-lab y. Needs no GPU."""
    objs = [json.load(open(p, encoding="utf-8")) for p in paths]
    per_model = [{r["dms_id"]: r for r in o.get("results", []) if "var" in r} for o in objs]
    shared = set(per_model[0])
    for d in per_model[1:]:
        shared &= set(d)
    results = []
    for dms_id in sorted(shared):
        rows_by_model = [{k: x for k, x, _y in d[dms_id]["var"]} for d in per_model]
        y_by_key = {k: y for k, _x, y in per_model[0][dms_id]["var"]}
        combined = combine_variant_scores(rows_by_model)
        if len(combined) < 20:
            continue
        keys = sorted(combined)
        xs = [combined[k] for k in keys]
        ys = [y_by_key[k] for k in keys]
        rng = np.random.default_rng(0)
        yss = np.array(ys); rng.shuffle(yss)
        results.append({"dms_id": dms_id, "n": len(xs), "rho": spearman(xs, ys),
                        "rho_shuf": spearman(xs, yss)})
    s = _summary(results)
    s["ok"] = bool(s["n_assays"] and s["median_abs_rho"] >= PASS_BAR and s["median_abs_shuf"] < SHUFFLE_MAX)
    return {"ensemble_of": list(paths), "results": results, **s}


# ---- J2 Phase-2 shared pure helpers (drift-guarded against scripts/esm_zeroshot_dms.py) ----
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


# ---- ESM-2 masked-marginals (copied faithfully from scripts/esm_zeroshot_dms.py::score_assay) ----
def score_assay(seq, dms, tok, model, device, mask_id, aa_ids, batch, maxlen, torch,
                long_mode="skip", keep_scores=False):
    """ESM-2 masked-marginals: score = logP(mut|ctx) - logP(wt|ctx) per variant.

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
        if L <= maxlen:
            enc = tok(seq, return_tensors="pt")
            ids = enc["input_ids"].to(device)
            for i in range(0, len(positions), batch):
                chunk = positions[i:i + batch]
                stack = ids.repeat(len(chunk), 1)
                for r, p in enumerate(chunk):
                    stack[r, p] = mask_id                       # token index p == 1-based residue p
                logits = model(stack).logits
                for r, p in enumerate(chunk):
                    logp_at[p] = torch.log_softmax(logits[r, p], dim=-1).float().cpu()
        else:                                                  # windowed long-protein path (Phase 2)
            for i in range(0, len(positions), batch):
                chunk = positions[i:i + batch]
                wins = [window_for_position(p, L, maxlen) for p in chunk]
                seqs = [seq[s:e] for (s, e, _) in wins]         # each exactly maxlen residues
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
    ap.add_argument("--dtype", choices=["float32", "float16"], default="float32",
                    help="float16 fits ESM2-3B on a free T4/P100 (Phase 2 — the certain lift over 650M)")
    ap.add_argument("--long-mode", choices=["skip", "window"], default="skip",
                    help="skip = baseline (drop >maxlen proteins, the published-0.48 behavior); "
                         "window = Phase 2 (score long proteins in a centered maxlen window)")
    ap.add_argument("--keep-scores", action="store_true",
                    help="write per-variant [key,x,y] rows into --out (needed for --ensemble-merge)")
    ap.add_argument("--ensemble-merge", nargs="+", default=None,
                    help="ENSEMBLE across model runs (each written with --keep-scores); no GPU; then exit")
    ap.add_argument("--self-test", action="store_true",
                    help="no-GPU pre-run check: is the ProteinGym data attached + how many assays scorable?")
    args = ap.parse_args(argv)

    if args.self_test:
        refs = load_reference(args.data_dir)
        rep = self_test_report(refs, args.maxlen, os.path.exists)
        print(json.dumps(rep, indent=2))
        print(f"\nSELF-TEST: {rep['have_csv']}/{rep['total']} assays have their per-assay CSV attached; "
              f"{rep['long_gt_maxlen']} exceed maxlen={args.maxlen} (recovered only with --long-mode window). "
              + ("READY to run on a GPU kernel." if rep["have_csv"]
                 else "NO PER-ASSAY DATA — attach the ProteinGym CSVs (see the runbook Step 0)."))
        return 0 if rep["have_csv"] else 1

    if args.merge:
        m = merge(args.merge)
        print(json.dumps({k: v for k, v in m.items() if k != "results"}, indent=2))
        tag = "  (>= 0.48 — matches the field!)" if m["median_abs_rho"] >= STRETCH else ""
        print(f"\nMERGED: {m['n_assays']} assays, median |Spearman| = {m['median_abs_rho']:.3f}{tag} "
              f"(shuffled {m['median_abs_shuf']:.3f}) -> {'PASS' if m['ok'] else 'below 0.45'}")
        return 0 if m["ok"] else 1

    if args.ensemble_merge:
        m = ensemble_merge(args.ensemble_merge)
        print(json.dumps({k: v for k, v in m.items() if k != "results"}, indent=2))
        tag = "  (>= 0.48 — matches/beats the field!)" if m["median_abs_rho"] >= STRETCH else ""
        print(f"\nENSEMBLE: {m['n_assays']} assays, median |Spearman| = {m['median_abs_rho']:.3f}{tag} "
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
    dt = torch.float16 if args.dtype == "float16" else torch.float32
    model = AutoModelForMaskedLM.from_pretrained(args.model, torch_dtype=dt).to(device).eval()
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
        r, why = score_assay(seq, dms, tok, model, device, mask_id, aa_ids, args.batch, args.maxlen, torch,
                             long_mode=args.long_mode, keep_scores=args.keep_scores)
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
