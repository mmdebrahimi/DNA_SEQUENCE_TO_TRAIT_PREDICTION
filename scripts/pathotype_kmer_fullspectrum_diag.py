"""Full-spectrum k-mer diagnostic for EP-4 (ExPEC vs EPEC).

WHY: the canonical bake-off (scripts/pathotype_kmer_bakeoff.py) uses top_n=10_000
MOST-FREQUENT 8-mers and landed AUROC=0.514 (chance). top-N-by-frequency selects
the conserved CORE genome — exactly the features that can't separate two E. coli
pathotypes. The discriminating signal (virulence loci) lives in the rare ACCESSORY
genome, which frequency-ranking discards. This script removes the top-N cap to
disambiguate:

  full-spectrum >= ~0.7  -> signal IS in the accessory genome (the 10k baseline was
                           under-powered by construction) -> pooled-NT arm warranted,
                           held to an ABSOLUTE >=0.65 bar (not +3pp over chance).
  full-spectrum ~ chance -> no genome-wide k-mer signal -> pooled representations
                           (k-mer OR NT mean-pool, tet precedent) likely futile ->
                           pivot to gene-targeted virulence detection.

Two arms: full-spectrum COUNTS and full-spectrum PRESENCE/ABSENCE. Same cached
per-genome counts, same seeded XGBoost, same LOSO folds as the canonical run.
"""
from __future__ import annotations
import json, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score

from scripts.pathotype_kmer_bakeoff import (
    load_strains,
    precompute_counts,
    run_kmer_loso_cached,
)

K = 8
FULL_TOP_N = 4 ** K  # 65536 -> captures ALL observed 8-mers (no top-N cap)
BASELINE_TOP10K_AUROC = 0.5138888888888888  # canonical bake-off result, for delta


def _score(yt: np.ndarray, ys: np.ndarray) -> tuple[float, int, bool]:
    auroc = float(roc_auc_score(yt, ys)) if len(set(yt.tolist())) == 2 else float("nan")
    uniq = sorted(set(round(float(s), 4) for s in ys))
    return auroc, len(uniq), len(uniq) <= 2


def main() -> int:
    seqs_by_strain, labels_by_strain, ids = load_strains()
    npos = sum(labels_by_strain.values()); nneg = len(ids) - npos
    print(f"[diag] {len(ids)} usable ({nneg} ExPEC / {npos} EPEC)")
    if npos < 3 or nneg < 3:
        print("[diag] ABORT: too few per class"); return 1

    print(f"[diag] precomputing per-genome k-mer counts (k={K}, once)...", flush=True)
    cache = precompute_counts(seqs_by_strain, k=K)
    distinct = sorted({km for d in cache.values() for km in d})
    print(f"[diag] full 8-mer space observed: {len(distinct)} / {FULL_TOP_N}")

    arms = {}
    for name, binary in (("fullspectrum_counts", False), ("fullspectrum_presence", True)):
        print(f"[diag] arm={name} (binary={binary}) running LOSO...", flush=True)
        yt, ys = run_kmer_loso_cached(seqs_by_strain, labels_by_strain, ids,
                                      drug=f"expec_vs_epec_{name}", k=K,
                                      top_n=FULL_TOP_N, cache=cache, binary=binary)
        auroc, ndistinct, degenerate = _score(yt, ys)
        delta = auroc - BASELINE_TOP10K_AUROC
        print(f"[diag] arm={name}: AUROC={auroc:.4f} (delta vs top10k={delta:+.4f}) "
              f"distinct={ndistinct} degenerate={degenerate}")
        arms[name] = {"auroc": auroc, "distinct_scores": ndistinct,
                      "degenerate": degenerate, "delta_vs_top10k": delta}

    best = max(a["auroc"] for a in arms.values() if not np.isnan(a["auroc"]))
    if best >= 0.70:
        verdict = "ACCESSORY_SIGNAL: full-spectrum recovers signal -> top10k was under-powered -> run NT arm (absolute >=0.65 bar)"
    elif best >= 0.60:
        verdict = "WEAK_SIGNAL: some accessory signal -> NT arm plausible but marginal; weigh against tet precedent"
    else:
        verdict = "NO_POOLED_SIGNAL: full-spectrum still ~chance -> pooled representations likely futile -> pivot to gene-targeted virulence detection"
    print(f"[diag] VERDICT: {verdict}")

    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "n": len(ids),
           "n_expec": nneg, "n_epec": npos, "k": K, "cv": "loso",
           "distinct_8mers_observed": len(distinct), "full_kmer_space": FULL_TOP_N,
           "baseline_top10k_auroc": BASELINE_TOP10K_AUROC,
           "arms": arms, "best_auroc": best, "verdict": verdict,
           "confound": "study-confounded (ExPEC=Salipante, EPEC=Hazen); 0.514 top10k already ruled out trivial batch separability",
           "ids": ids}
    out = REPO / "research_outputs/pathotype_kmer_fullspectrum_diag_2026-06-02.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[diag] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
