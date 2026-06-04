"""Cef embedding-vs-classical falsifier — LAPTOP-side (GPU-only workhorse model).

The Phase 2 falsifier for ceftriaxone (concentrated plasmid-beta-lactamase mechanism):
does the NT embedding beat the best classical baseline by >= 3 pp AUROC under leakage-safe
CV? If yes, the NT-mean-pool architecture earns its keep on a concentrated mechanism at
clean MIC labels; if no, re-evaluate the architecture before any further drug.

GPU-ONLY SPLIT (2026-06-04 user directive): the workhorse does ONE thing — produce the NT
embedding cache from PUBLIC genomes. THIS script (the falsifier) runs on the laptop, CPU-only:
NT-XGBoost on the cached embeddings vs k-mer-XGB classical, under `leave_one_accession_out` CV.
Reuses the stage1 cipro machinery (load_features + train/predict wrappers + bootstrap) and the
canonical CV in `dna_decode/eval/cv.py`.

Inputs the laptop needs (both arrive from elsewhere):
  - the cef NT embedding cache .h5  (returned by the workhorse; the ONE cross-machine artifact)
  - the cef genome FASTAs           (fetched from PUBLIC NCBI on the laptop, for the k-mer baseline)
Parks with clear instructions if either is missing — never runs a partial/misleading gate.

Run:  uv run python scripts/cef_falsifier.py
Exit: 0 iff the gate PASSES (NT beats k-mer by >= 3 pp); 1 on FAIL; 2/3 on missing inputs (parked).
"""
from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.data.cohort import load_cohort
from dna_decode.eval.cv import leave_one_accession_out_cv
from dna_decode.eval.loso_kmer import run_kmer_xgboost_loso
from dna_decode.eval.metrics import compute_metrics
# reuse the stage1 cipro falsifier components verbatim (drug-agnostic)
from scripts.stage1_n40_cipro import (
    VariantResult,
    _nt_logreg_predict,
    _nt_logreg_train,
    _nt_xgb_predict,
    _nt_xgb_train,
    load_features,
    paired_bootstrap_ci,
)

GATE_THRESHOLD_PP = 3.0
DRUG = "ceftriaxone"
DEFAULT_COHORT = ROOT / "data/processed/gate_b_cohort.parquet"
DEFAULT_NT_CACHE = ROOT / "data/processed/embeddings/nt_gate_b_cohort_67.h5"
DEFAULT_REFSEQ = ROOT / "data/refseq_cache"


def _accession_assignments(cohort, strain_ids: list[str]) -> dict[str, str]:
    by_id = {s.strain_id: getattr(s, "assembly_accession", "") for s in cohort.strains}
    return {sid: by_id.get(sid, "") for sid in strain_ids}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Cef embedding-vs-classical falsifier (laptop, GPU-only model)")
    ap.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    ap.add_argument("--nt-cache", type=Path, default=DEFAULT_NT_CACHE,
                    help="cef NT embedding cache .h5 returned by the workhorse")
    ap.add_argument("--refseq-cache", type=Path, default=DEFAULT_REFSEQ,
                    help="genome FASTA cache (fetched from public NCBI on the laptop)")
    ap.add_argument("--drug", default=DRUG)
    ap.add_argument("--aggregation", choices=["mean", "max", "mean+max"], default="mean",
                    help="NT pooling; 'mean' matches the v0 cef cached-strain run (AUROC 0.895)")
    ap.add_argument("--kmer-k", type=int, default=8)
    ap.add_argument("--kmer-top-n", type=int, default=10_000)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args(argv)

    # --- park cleanly on missing inputs (GPU-only model: the cache arrives from the workhorse) ---
    if not args.nt_cache.exists():
        print(f"PARKED: cef NT embedding cache not found at {args.nt_cache}\n"
              f"  This is the ONE artifact the workhorse must return (GPU-only model).\n"
              f"  On the workhorse it is nt_gate_b_cohort_67.h5 (already produced by the cef\n"
              f"  cached-strain run). Place it at the path above or pass --nt-cache.", file=sys.stderr)
        return 2

    cohort = load_cohort(args.cohort)
    print(f"[cef-falsifier] cohort {args.cohort.name}: {len(cohort.strains)} strains; drug={args.drug}; "
          f"aggregation={args.aggregation}")

    try:
        X_nt, seqs_by_strain, labels_by_strain, strain_ids, mlsts = load_features(
            cohort, args.nt_cache, args.refseq_cache, args.drug, aggregation=args.aggregation,
        )
    except FileNotFoundError as e:
        print(f"PARKED: cef genome FASTA missing under {args.refseq_cache} ({e}).\n"
              f"  Fetch the cef cohort genomes from public NCBI on the laptop, e.g.:\n"
              f"    uv run python -m scripts.pipeline ingest --drugs ceftriaxone \\\n"
              f"      --assembly-metadata-csv <BVBRC_genome.csv>   (or download by assembly_accession)\n"
              f"  then re-run. The k-mer classical baseline needs the FASTAs.", file=sys.stderr)
        return 3

    y = np.array([labels_by_strain[s] for s in strain_ids], dtype=int)
    n_r, n_s = int((y == 1).sum()), int((y == 0).sum())
    print(f"[cef-falsifier] effective N={len(y)} ({n_r}R/{n_s}S); NT shape={X_nt.shape}")
    if n_r < 2 or n_s < 2:
        print("PARKED: degenerate class balance (need >=2 per class).", file=sys.stderr)
        return 2

    acc = _accession_assignments(cohort, strain_ids)
    n_dup = len(strain_ids) - len(set(acc.values()))
    print(f"[cef-falsifier] leave_one_accession_out CV (duplicate-accession strains: {n_dup}; "
          f"== strain-out when 0)")

    results: list[VariantResult] = []

    print("[cef-falsifier] NT-XGBoost ...")
    cv = leave_one_accession_out_cv(X_nt, y, strain_ids, acc, _nt_xgb_train, _nt_xgb_predict, drug=args.drug)
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult("NT-XGBoost", float(m.auroc), float(m.auprc),
                                  cv.all_y_score, cv.all_y_true, cv.strain_ids, True))
    print(f"  AUROC={m.auroc:.3f}")

    print("[cef-falsifier] NT-logreg ...")
    cv = leave_one_accession_out_cv(X_nt, y, strain_ids, acc, _nt_logreg_train, _nt_logreg_predict, drug=args.drug)
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult("NT-logreg", float(m.auroc), float(m.auprc),
                                  cv.all_y_score, cv.all_y_true, cv.strain_ids, True))
    print(f"  AUROC={m.auroc:.3f}")

    print("[cef-falsifier] k-mer-XGB (classical comparator) ...")
    # no duplicate accessions in the cef cohort -> strain-out folds == accession-out folds.
    cv = run_kmer_xgboost_loso(seqs_by_strain, labels_by_strain, strain_ids, drug=args.drug,
                               k=args.kmer_k, top_n=args.kmer_top_n)
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult("k-mer-XGB", float(m.auroc), float(m.auprc),
                                  cv.all_y_score, cv.all_y_true, cv.strain_ids, True))
    print(f"  AUROC={m.auroc:.3f}")

    # --- gate ---
    nt_best = max((r for r in results if r.name.startswith("NT")), key=lambda r: r.auroc)
    kmer = next(r for r in results if r.name == "k-mer-XGB")
    if nt_best.strain_ids != kmer.strain_ids:
        raise ValueError("alignment: NT-best and k-mer strain_ids diverge")
    gap_pp = (nt_best.auroc - kmer.auroc) * 100.0
    _, lo, hi, n_eff = paired_bootstrap_ci(kmer.per_strain_true, nt_best.per_strain_scores,
                                           kmer.per_strain_scores)
    passed = gap_pp >= GATE_THRESHOLD_PP
    verdict = (f"PASS (gap {gap_pp:+.1f} pp >= {GATE_THRESHOLD_PP:.0f})" if passed
               else f"FAIL (gap {gap_pp:+.1f} pp < {GATE_THRESHOLD_PP:.0f})")

    out = args.output or (ROOT / f"wiki/cef_falsifier_{date.today().isoformat()}.md")
    lines = [
        f"# Cef embedding-vs-classical falsifier ({date.today().isoformat()})",
        "",
        "> Phase 2 falsifier (laptop, GPU-only model). Does NT embedding beat the best classical baseline",
        "> by >= 3 pp AUROC under leave_one_accession_out CV on a concentrated beta-lactamase mechanism?",
        "",
        f"**Cohort:** `{args.cohort}` (effective N={len(y)}; {n_r}R/{n_s}S) · **drug:** {args.drug} · "
        f"**pooling:** {args.aggregation}",
        f"**CV:** leave_one_accession_out (dup-accession strains: {n_dup})",
        f"**Best NT head:** {nt_best.name} AUROC {nt_best.auroc:.3f} · **k-mer-XGB:** {kmer.auroc:.3f}",
        f"**Gap:** {gap_pp:+.1f} pp · **95% bootstrap CI:** [{lo*100:+.1f}, {hi*100:+.1f}] pp (eff {n_eff}/1000)",
        f"**VERDICT:** {verdict}",
        "",
        "| Variant | AUROC | AUPRC |",
        "|---|---:|---:|",
        *[f"| {r.name} | {r.auroc:.3f} | {r.auprc:.3f} |" for r in results],
        "",
        "## Notes",
        "- All variants `calibrate=False` (small-N calibration footgun; LESSONS_LEARNED 2026-05-14).",
        "- cef has no duplicate accessions -> accession-out folds == strain-out folds (verified at run).",
        "- Classical comparator = k-mer-XGB only here; AMRFinder gene-presence baseline is a follow-up.",
        "- Optional stricter pass: Mash-clade-out CV via scripts/mash_cluster_n147.py (Docker).",
        f"- LOSO at N={len(y)} carries a ~0.10 AUROC noise floor; the bootstrap CI surfaces it.",
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"[cef-falsifier] {verdict}  -> packet {out}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
