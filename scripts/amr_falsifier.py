"""Drug-agnostic AMR embedding-vs-classical falsifier (laptop, GPU-only model).

Generalizes the cef-specific runner. For ANY drug + cohort + NT embedding cache, asks: does the
frozen-NT mean-pool embedding beat the best classical baseline by >= 3 pp AUROC under leakage-safe
CV — on a cohort that is actually DE-CONFOUNDED?

Two preconditions, both hard:
  1. COHORT DE-CONFOUND GATE (dna_decode/eval/cohort_deconfound): refuse a promotable verdict on a
     CONFOUNDED cohort (R and S near-separable by lineage/geography → a verdict measures batch, not
     biology). This is the precondition the cef run lacked — pathotype + cef both died on it.
  2. INPUTS PRESENT: NT cache (one .h5, from the workhorse) + genome FASTAs (public NCBI, local).

CI-AWARE VERDICT (fixes the cef brainstorm finding): a point gap >= 3 pp is necessary but NOT
sufficient at small N — PASS requires the paired-bootstrap CI lower bound > 0 too; else NOISY.

Run (cipro N=147 clean substrate):
  uv run python scripts/amr_falsifier.py --drug ciprofloxacin \
    --cohort data/processed/stage2_n150_cipro_cohort.parquet \
    --nt-cache data/processed/embeddings/nt_n147_cipro.h5
Exit: 0 PASS(promotable) · 1 FAIL · 2/3 missing inputs (parked) · 4 CONFOUNDED substrate (blocked) ·
5 NOISY · 6 NON-PROMOTABLE (WARN screen or --allow-confounded diagnostic — never a clean PASS).
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
from dna_decode.eval.cohort_deconfound import (
    CONFOUNDED, DE_CONFOUNDED, confound_report_for_cohort, render_report,
)
from dna_decode.eval.cv import leave_one_accession_out_cv
from dna_decode.eval.loso_kmer import run_kmer_xgboost_loso
from dna_decode.eval.metrics import compute_metrics
from scripts.stage1_n40_cipro import (
    VariantResult, _nt_logreg_predict, _nt_logreg_train, _nt_xgb_predict, _nt_xgb_train,
    load_features, paired_bootstrap_ci,
)

GATE_THRESHOLD_PP = 3.0


def _accession_assignments(cohort, strain_ids):
    by_id = {s.strain_id: getattr(s, "assembly_accession", "") for s in cohort.strains}
    return {sid: by_id.get(sid, "") for sid in strain_ids}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Drug-agnostic AMR embedding-vs-classical falsifier")
    ap.add_argument("--drug", required=True)
    ap.add_argument("--cohort", type=Path, required=True)
    ap.add_argument("--nt-cache", type=Path, required=True)
    ap.add_argument("--refseq-cache", type=Path, default=ROOT / "data/refseq_cache")
    ap.add_argument("--aggregation", choices=["mean", "max", "mean+max"], default="mean")
    ap.add_argument("--kmer-k", type=int, default=8)
    ap.add_argument("--kmer-top-n", type=int, default=10_000)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--allow-confounded", action="store_true",
                    help="override the de-confound gate (emits a non-promotable diagnostic only)")
    args = ap.parse_args(argv)

    cohort = load_cohort(args.cohort)

    # --- PRECONDITION 1: cohort de-confound gate (BEFORE touching the GPU cache or genomes) ---
    rep = confound_report_for_cohort(cohort, args.drug)
    print(f"[falsifier] de-confound gate: {render_report(rep)}")
    if rep["verdict"] == CONFOUNDED and not args.allow_confounded:
        out = args.output or (ROOT / f"wiki/{args.drug}_falsifier_{date.today().isoformat()}.md")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            f"# {args.drug} falsifier — BLOCKED: CONFOUNDED_SUBSTRATE ({date.today().isoformat()})\n\n"
            f"The de-confound gate refused this cohort: **{rep['reason']}**\n\n"
            f"`{render_report(rep)}`\n\n"
            "An embedding-vs-classical verdict here would measure lineage/geography/batch, not biology "
            "(the pathotype/cef study==class trap). Rebuild a cohort where R and S co-occur within "
            "multiple lineages AND geographies, or pass --allow-confounded for a NON-PROMOTABLE "
            "diagnostic only. No verdict computed.\n", encoding="utf-8")
        print(f"[falsifier] BLOCKED (CONFOUNDED_SUBSTRATE) -> {out}; no verdict computed.")
        return 4

    # promotability contract: ONLY a DE_CONFOUNDED screen (and no --allow-confounded override) may
    # back a promotable PASS. WARN or --allow-confounded → diagnostic run, NON-PROMOTABLE (exit 6).
    promotable = bool(rep.get("promotable")) and not args.allow_confounded
    if not promotable:
        print(f"[falsifier] NON-PROMOTABLE run (cohort screen={rep['verdict']}"
              f"{'; --allow-confounded override' if args.allow_confounded else ''}) — diagnostic only.")

    # --- PRECONDITION 2: inputs present ---
    if not args.nt_cache.exists():
        print(f"PARKED: NT embedding cache not found at {args.nt_cache} (GPU-only: due from the workhorse).",
              file=sys.stderr)
        return 2
    try:
        X_nt, seqs, labels_by, strain_ids, mlsts = load_features(
            cohort, args.nt_cache, args.refseq_cache, args.drug, aggregation=args.aggregation)
    except FileNotFoundError as e:
        print(f"PARKED: genome FASTA missing under {args.refseq_cache} ({e}); fetch public NCBI genomes first.",
              file=sys.stderr)
        return 3

    y = np.array([labels_by[s] for s in strain_ids], dtype=int)
    n_r, n_s = int((y == 1).sum()), int((y == 0).sum())
    print(f"[falsifier] effective N={len(y)} ({n_r}R/{n_s}S); NT shape={X_nt.shape}")
    if n_r < 2 or n_s < 2:
        print("PARKED: degenerate class balance.", file=sys.stderr)
        return 2

    acc = _accession_assignments(cohort, strain_ids)
    results = []
    for name, tr, pr, gb in [("NT-XGBoost", _nt_xgb_train, _nt_xgb_predict, True),
                             ("NT-logreg", _nt_logreg_train, _nt_logreg_predict, True)]:
        cv = leave_one_accession_out_cv(X_nt, y, strain_ids, acc, tr, pr, drug=args.drug)
        m = compute_metrics(cv.all_y_true, cv.all_y_score)
        results.append(VariantResult(name, float(m.auroc), float(m.auprc),
                                     cv.all_y_score, cv.all_y_true, cv.strain_ids, gb))
        print(f"[falsifier] {name} AUROC={m.auroc:.3f}")
    cvk = run_kmer_xgboost_loso(seqs, labels_by, strain_ids, drug=args.drug, k=args.kmer_k, top_n=args.kmer_top_n)
    mk = compute_metrics(cvk.all_y_true, cvk.all_y_score)
    results.append(VariantResult("k-mer-XGB", float(mk.auroc), float(mk.auprc),
                                 cvk.all_y_score, cvk.all_y_true, cvk.strain_ids, True))
    print(f"[falsifier] k-mer-XGB AUROC={mk.auroc:.3f}")

    nt_best = max((r for r in results if r.name.startswith("NT")), key=lambda r: r.auroc)
    kmer = next(r for r in results if r.name == "k-mer-XGB")
    gap = (nt_best.auroc - kmer.auroc) * 100.0
    # ALIGN per-strain scores by strain_id before the paired bootstrap. NT uses
    # leave_one_accession_out_cv (fold order = accession grouping) while k-mer uses plain LOSO
    # (fold order = strain_ids input) — same strains, different EMISSION ORDER. AUROC is
    # order-independent (the 0.914/0.824 stand); only the paired CI needs a common order.
    nt_by = {s: (sc, t) for s, sc, t in zip(nt_best.strain_ids, nt_best.per_strain_scores, nt_best.per_strain_true)}
    km_by = {s: (sc, t) for s, sc, t in zip(kmer.strain_ids, kmer.per_strain_scores, kmer.per_strain_true)}
    common = [s for s in nt_best.strain_ids if s in km_by]
    if set(nt_by) != set(km_by):
        only = (set(nt_by) ^ set(km_by))
        print(f"[falsifier] WARN: NT/k-mer strain sets differ by {len(only)} strain(s); "
              f"pairing on the {len(common)}-strain intersection.")
    y_al = np.array([nt_by[s][1] for s in common], dtype=int)
    nt_al = np.array([nt_by[s][0] for s in common], dtype=np.float32)
    km_al = np.array([km_by[s][0] for s in common], dtype=np.float32)
    _, lo, hi, n_eff = paired_bootstrap_ci(y_al, nt_al, km_al)

    # CI-aware verdict (brainstorm fix): point >= 3pp is necessary; ci_lo > 0 makes it a PASS.
    # A non-promotable run NEVER returns a clean PASS, regardless of the gap (exit 6).
    if not promotable:
        verdict, rc = (f"NON-PROMOTABLE diagnostic (cohort screen {rep['verdict']}): "
                       f"gap {gap:+.1f}pp, CI [{lo*100:+.1f},{hi*100:+.1f}]pp — not a promotable verdict"), 6
    elif gap >= GATE_THRESHOLD_PP and lo * 100 > 0:
        verdict, rc = f"PASS (gap {gap:+.1f}pp, CI lo {lo*100:+.1f}>0)", 0
    elif gap >= GATE_THRESHOLD_PP:
        verdict, rc = f"NOISY (gap {gap:+.1f}pp but CI lo {lo*100:+.1f}<=0 — not separable from noise)", 5
    else:
        verdict, rc = f"FAIL (gap {gap:+.1f}pp < {GATE_THRESHOLD_PP:.0f})", 1

    out = args.output or (ROOT / f"wiki/{args.drug}_falsifier_{date.today().isoformat()}.md")
    lines = [
        f"# {args.drug} embedding-vs-classical falsifier ({date.today().isoformat()})", "",
        f"**De-confound gate:** `{render_report(rep)}`",
        f"**Cohort:** `{args.cohort}` (N={len(y)}; {n_r}R/{n_s}S) · pooling {args.aggregation} · CV leave_one_accession_out",
        f"**Best NT:** {nt_best.name} {nt_best.auroc:.3f} · **k-mer-XGB:** {kmer.auroc:.3f} · **gap {gap:+.1f}pp**",
        f"**95% bootstrap CI on gap:** [{lo*100:+.1f}, {hi*100:+.1f}]pp (eff {n_eff}/1000)",
        f"**VERDICT:** {verdict}", "",
        "| Variant | AUROC | AUPRC |", "|---|---:|---:|",
        *[f"| {r.name} | {r.auroc:.3f} | {r.auprc:.3f} |" for r in results], "",
        "## Notes",
        "- De-confound gate is a PRECONDITION (CONFOUNDED cohort → blocked, no verdict).",
        "- CI-aware verdict: point gap >= 3pp AND bootstrap CI lower bound > 0 for a PASS (else NOISY).",
        "- Classical comparator = k-mer-XGB; AMRFinder gene-presence baseline is a follow-up (stronger for concentrated mechanisms).",
        "- calibrate=False; verify_complete cache integrity = follow-up.",
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"[falsifier] {verdict} -> {out}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
