"""Windowed-NT arm of the EP-4 representation bake-off (ExPEC vs EPEC).

Tests the foundation-model premise against the classical bar set by the
full-spectrum k-mer diagnostic (2026-06-02):
    top-10k counts        = 0.514 (chance; under-powered by core-genome bias)
    full-spectrum counts  = 0.604
    full-spectrum presence/absence = 0.729  <- BEST CLASSICAL = the bar to beat
    -> per CLAUDE.md ">=3pp over classical", NT must reach >= ~0.76 LOSO to justify FM.

Design (per user decision 2026-06-02): whole-genome WINDOWED NT (no Bakta), tile
each genome into max_context windows, embed each window, then pool across windows
TWO ways:
  - mean-pool : soft-average; the analogue of the count representation that
                under-performed (0.604) -> PREDICT ~0.60 (dilution).
  - max-pool  : presence-aligned (fires on the most-activating window) -> the
                analogue of presence/absence (0.729) -> PREDICT higher.
If max >> mean, that corroborates "signal is localized gene presence" and tells
us mean-pooling (the cheap default) was the wrong aggregation.

Per-window embeddings are CACHED to C: (data/nt_windows/, gitignored) keyed by
accession + window params, so pooling variants and restarts never re-pay GPU.
The 24 genomes are the SAME cached ENA assemblies as the k-mer arm (load_strains).

CAVEAT: genomes are contig-concatenated with 100-N separators; windows spanning a
contig boundary contain a short N-run (negligible fraction). Windows that are
>50% N are skipped. N=24, study-confounded (Salipante vs Hazen) -> a win is
necessary-not-sufficient; lineage/ST-aware splits still needed to claim biology.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score

from scripts.pathotype_kmer_bakeoff import load_strains
from dna_decode.models.foundation import model_factory
from dna_decode.models.classifiers import predict_proba, train_xgboost_classifier

try:
    import torch as _torch  # for periodic empty_cache() to fight 4 GB fragmentation
except ImportError:
    _torch = None

# Attention memory is O(batch * heads * tokens^2). NT v2 tokenizes ~6 bp/token, BUT
# N-runs at contig boundaries break 6-mer alignment so some windows tokenize toward
# the model max (2048 tokens) despite the 6144 bp size — making per-batch memory
# spiky. On a 4 GB GTX 860M, CUDA context + fp32 weights already eat ~1.4 GB, leaving
# ~2.5 GB; batch 8 OOMs on a long-window batch. Batch 2 survives even a worst-case
# 2048-token window (~0.7 GB attention+hidden) and batching gives little speedup on
# this GPU anyway. empty_cache() between batches fights fragmentation.
WINDOW = 6144           # bp per window (well under 12288 max_context)
STRIDE = 6144           # non-overlapping tiles
# Batch 2 is the only TDR-safe config: the GTX 860M is ALSO the display GPU, so any
# CUDA kernel running >2 s triggers Windows Timeout Detection & Recovery (driver reset).
# Batch 8 forward kernels exceed 2 s; batch 2 @ ~1024 tokens stays under it (proven for
# 500+ windows). ~5-8 hr for 24 genomes; restartable via per-genome .npy cache.
GPU_BATCH = 2           # TDR-safe on the display GPU (kernel < 2 s)
EMPTY_CACHE_EVERY = 40  # forward passes between torch.cuda.empty_cache() calls
MAX_N_FRAC = 0.01       # skip windows with >1% N: removes contig-boundary windows whose N-runs
                        # break 6-mer tokenization into long, memory-spiky token sequences.
                        # Survivors are clean intra-contig ACGT -> uniform token length -> batchable.
NCACHE = Path("C:/Users/Farshad/PythonProjects/dna_decode/data/nt_windows")
BEST_CLASSICAL = 0.7291666666666665  # full-spectrum presence/absence (bar to beat)
FM_BAR = 0.76                        # best-classical + ~3pp (CLAUDE.md FM premise)


def tile(seq: str, window: int, stride: int) -> list[str]:
    out = []
    for start in range(0, max(1, len(seq) - window + 1), stride):
        w = seq[start:start + window]
        if w.upper().count("N") <= MAX_N_FRAC * len(w):   # keep clean intra-contig windows only
            out.append(w)
    if not out:  # very short / N-heavy genome fallback: take the whole thing
        out = [seq[:window]]
    return out


def genome_windows_embedding(model, acc: str, seq: str) -> np.ndarray:
    """(n_windows, dim) per-window NT embeddings, cached to C: as .npy."""
    NCACHE.mkdir(parents=True, exist_ok=True)
    cache = NCACHE / f"{acc}_w{WINDOW}_s{STRIDE}.npy"
    if cache.exists():
        arr = np.load(cache)
        if arr.ndim == 2 and arr.shape[0] > 0:
            return arr
    windows = tile(seq, WINDOW, STRIDE)
    rows = []
    for b, i in enumerate(range(0, len(windows), GPU_BATCH)):
        rows.append(model.embed_batch(windows[i:i + GPU_BATCH]))
        if _torch is not None and (b + 1) % EMPTY_CACHE_EVERY == 0:
            _torch.cuda.empty_cache()
        if (i // GPU_BATCH) % 25 == 0:
            print(f"    {acc}: {min(i + GPU_BATCH, len(windows))}/{len(windows)} windows", flush=True)
    arr = np.vstack(rows).astype(np.float32)
    np.save(cache, arr)
    return arr


def loso(X: np.ndarray, y: np.ndarray, drug: str) -> tuple[float, np.ndarray]:
    n = len(y)
    proba = np.zeros(n, dtype=float)
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        clf = train_xgboost_classifier(X[tr], y[tr], drug_name=drug, calibrate=False)
        proba[i] = float(predict_proba(clf, X[i:i + 1])[0])
    auroc = float(roc_auc_score(y, proba)) if len(set(y.tolist())) == 2 else float("nan")
    return auroc, proba


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--limit", type=int, default=0, help="smoke: embed only first N genomes")
    args = ap.parse_args(argv)

    seqs_by_strain, labels_by_strain, ids = load_strains()
    if args.limit:
        ids = ids[:args.limit]
    npos = sum(labels_by_strain[s] for s in ids); nneg = len(ids) - npos
    print(f"[nt] {len(ids)} strains ({nneg} ExPEC / {npos} EPEC); window={WINDOW} stride={STRIDE}")

    print("[nt] loading NT v2 100M...", flush=True)
    model = model_factory("nucleotide_transformer",
                          config_path=str(REPO / "config" / "datasources.yaml"),
                          device=args.device)

    Xmean, Xmax, y = [], [], []
    for sid in ids:
        print(f"[nt] embed {sid} ...", flush=True)
        win = genome_windows_embedding(model, sid.split("|")[0], seqs_by_strain[sid])
        Xmean.append(win.mean(axis=0)); Xmax.append(win.max(axis=0))
        y.append(int(labels_by_strain[sid]))
        print(f"[nt]   {sid}: {win.shape[0]} windows -> dim {win.shape[1]}", flush=True)
    Xmean = np.vstack(Xmean).astype(np.float32)
    Xmax = np.vstack(Xmax).astype(np.float32)
    y = np.array(y, dtype=int)

    if args.limit:
        print(f"[nt] SMOKE ok: embedded {len(ids)} genomes, Xmean={Xmean.shape}. Skipping LOSO.")
        return 0

    results = {}
    for name, X in (("mean_pool", Xmean), ("max_pool", Xmax)):
        auroc, _ = loso(X, y, drug=f"expec_vs_epec_nt_{name}")
        delta_cls = auroc - BEST_CLASSICAL
        print(f"[nt] {name}: AUROC={auroc:.4f} | vs best-classical(0.729)={delta_cls:+.4f} "
              f"| clears FM bar(0.76)={auroc >= FM_BAR}", flush=True)
        results[name] = {"auroc": auroc, "delta_vs_best_classical": delta_cls,
                         "clears_fm_bar": bool(auroc >= FM_BAR)}

    best = max(r["auroc"] for r in results.values())
    if best >= FM_BAR:
        verdict = "FM_JUSTIFIED: NT clears best-classical+3pp -> foundation-model premise holds on this contrast (subject to confound resolution)"
    elif best >= BEST_CLASSICAL:
        verdict = "FM_MARGINAL: NT matches/slightly beats classical but not by 3pp -> not worth the FM complexity yet"
    else:
        verdict = "FM_FAILS: NT below best classical (0.729) -> pooled NT dilutes the localized signal (tet precedent confirmed); pivot to gene-targeted detection"
    mean_lt_max = results["max_pool"]["auroc"] > results["mean_pool"]["auroc"]
    print(f"[nt] VERDICT: {verdict}")
    print(f"[nt] dilution prediction (max>mean): {mean_lt_max}")

    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "n": len(ids),
           "n_expec": nneg, "n_epec": npos, "model": "nucleotide_transformer_v2_100m",
           "representation": "whole_genome_windowed", "window_bp": WINDOW, "stride_bp": STRIDE,
           "cv": "loso", "best_classical_bar": BEST_CLASSICAL, "fm_bar_3pp": FM_BAR,
           "arms": results, "best_auroc": best, "verdict": verdict,
           "dilution_prediction_max_gt_mean": bool(mean_lt_max),
           "caveat": "N=24, study-confounded; contig-concat 100-N separators; non-overlapping tiles; >50%-N windows skipped",
           "ids": ids}
    out = REPO / "research_outputs/pathotype_nt_windowed_arm_2026-06-02.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[nt] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
