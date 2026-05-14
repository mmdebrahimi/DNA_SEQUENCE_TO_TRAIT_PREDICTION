"""Stage 1 Phase 2 engineering screen — N=40 (effective 38) cipro cohort.

Locked criterion (2026-05-14 user lock):
  PASS iff max(NT-XGBoost, NT-logreg) AUROC >= k-mer-XGB AUROC + 3 pp under LOSO.

4-experiment matrix (per /brainstorm 2026-05-14):
  - NT-XGBoost (primary gate-bearing head; matches 12-strain smoke)
  - NT-logreg  (sanity-check baseline; H13 plumbing check)
  - k-mer-XGB  (classical comparator)
  - NT+k-mer-fusion-logreg (DIAGNOSTIC ONLY -- NOT gate-bearing)

All four use `calibrate=False` so AUROC measures representation quality, not
calibration-wrapper behavior (small-N calibration footgun; see LESSONS_LEARNED
2026-05-14).

Per `plans/Stage1_Refactor_And_Test_Hardening_Plan.md`:
  - LOSO over NT variants reuses `leave_one_strain_out_cv` from `dna_decode/eval/cv.py`
  - k-mer + fusion LOSO factored into `dna_decode/eval/loso_kmer.py`
  - Logreg via `_train_baseline_logreg(..., calibrate=False)` (explicit calibration-skip
    branch, NOT the default CalibratedClassifierCV wrapping)
  - Verdict semantics frozen as a pure function of point gap; `stage2_action` is the
    decision-layer (BURST_STAGE_2 / HOLD_STAGE_2_CI_DEGENERATE / ALTERNATIVE_POOLING_RERUN /
    PIVOT_TO_BAKTA) computed from (verdict, ci_lo, fusion behavior).
  - `compute_gate_outcome` validates `strain_ids` alignment: raises on NT-vs-k-mer
    mismatch; suppresses fusion note (without raising) on fusion mismatch.
  - `ClassifierTrainingError` from k-mer/fusion paths re-raises (no silent mean-fallback).

Diagnostic appendix in the result packet:
  - Unique MLST count + per-MLST R/S split
  - Per-strain LOSO predictions table with MLST labels
  - Paired bootstrap CI (B=1000, with n_effective surfacing) on the gap NT-best - k-mer-XGB
  - 3-bucket verdict (pure function of point gap) + stage2_action (decision-layer output)

NOT a powered statistical comparison -- Stage 1 is an engineering screen for
go/no-go on Stage 2 Databricks burst budget. Stage 2 (N=150, Mash-clade-out CV,
Option-C threshold >=5 pp + biology check) is the real ship gate.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Literal

import numpy as np

from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path
from dna_decode.eval.cv import CVResult, leave_one_strain_out_cv
from dna_decode.eval.loso_kmer import run_fusion_loso, run_kmer_xgboost_loso
from dna_decode.eval.metrics import compute_metrics
from dna_decode.models.cache import EmbeddingCache
from dna_decode.models.classical_baselines import CONTIG_SEPARATOR, _train_baseline_logreg
from dna_decode.models.classifiers import (
    aggregate_strain_features,
    predict_proba,
    train_xgboost_classifier,
)


GATE_THRESHOLD_PP = 3.0
BOOTSTRAP_ITERATIONS = 1000
BOOTSTRAP_SEED = 42


Stage2Action = Literal[
    "BURST_STAGE_2",
    "HOLD_STAGE_2_CI_DEGENERATE",
    "ALTERNATIVE_POOLING_RERUN",
    "PIVOT_TO_BAKTA",
]


@dataclass
class VariantResult:
    name: str
    auroc: float
    auprc: float
    per_strain_scores: np.ndarray
    per_strain_true: np.ndarray
    strain_ids: list[str]
    is_gate_bearing: bool
    extra: dict = field(default_factory=dict)


def _load_fasta_contigs(fasta_p: Path) -> list[str]:
    from Bio import SeqIO

    return [str(record.seq).upper() for record in SeqIO.parse(str(fasta_p), "fasta")]


# --- training / prediction wrappers --------------------------------------------------

def _nt_xgb_train(X: np.ndarray, y: np.ndarray) -> object:
    return train_xgboost_classifier(X, y, drug_name="cipro_stage1_nt_xgb", calibrate=False)


def _nt_xgb_predict(clf, X: np.ndarray) -> np.ndarray:
    return predict_proba(clf, X)


def _nt_logreg_train(X: np.ndarray, y: np.ndarray):
    return _train_baseline_logreg(X, y, drug_name="cipro_stage1_nt_logreg", calibrate=False)


def _nt_logreg_predict(clf, X: np.ndarray) -> np.ndarray:
    return clf.model.predict_proba(X)[:, 1].astype(np.float32)


# --- feature loading -----------------------------------------------------------------

def load_features(
    cohort,
    nt_cache: Path,
    refseq_root: Path,
    drug: str,
) -> tuple[np.ndarray, dict[str, str], dict[str, int], list[str], list[str]]:
    """Load per-strain NT embeddings + N-concatenated genome strings.

    Returns:
        (X_nt, seqs_by_strain, labels_by_strain, strain_ids, mlsts)
        where X_nt[i] corresponds to strain_ids[i] (alignment contract).
        Stage 1's cohort is the NT-cache-present effective subset.

    Raises:
        ValueError: if any cohort strain is missing its MLST attribute. Stage 1's
        lineage diagnostic depends on MLST; collapsing missing values to "unknown"
        would inflate the largest_mlst_group calculation.
    """
    cache = EmbeddingCache(
        nt_cache,
        model_name="nucleotide_transformer",
        model_version="InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
        embedding_dim=512,
    )
    drug_lower = drug.lower()
    nt_rows: list[np.ndarray] = []
    seqs_by_strain: dict[str, str] = {}
    labels_by_strain: dict[str, int] = {}
    strain_ids: list[str] = []
    mlsts: list[str] = []
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        gene_ids = cache.list_genes(s.strain_id)
        if not gene_ids:
            continue
        gene_matrix = cache.bulk_get([(s.strain_id, g) for g in gene_ids])
        nt_rows.append(aggregate_strain_features(gene_matrix, "mean"))
        fp = fasta_path(s.assembly_accession, refseq_root)
        seqs_by_strain[s.strain_id] = CONTIG_SEPARATOR.join(_load_fasta_contigs(fp))
        labels_by_strain[s.strain_id] = int(s.ast_labels[drug_lower])
        strain_ids.append(s.strain_id)
        mlst = getattr(s, "mlst", None)
        if mlst is None:
            raise ValueError(
                f"strain {s.strain_id} missing MLST -- cohort metadata incomplete; "
                f"fix before Stage 1 (lineage diagnostic depends on MLST)"
            )
        mlsts.append(str(mlst))

    X_nt = np.stack(nt_rows)
    return X_nt, seqs_by_strain, labels_by_strain, strain_ids, mlsts


# --- bootstrap CI --------------------------------------------------------------------

def paired_bootstrap_ci(
    y_true: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    n_iterations: int = BOOTSTRAP_ITERATIONS,
    ci: float = 0.95,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[float, float, float, int]:
    """Paired bootstrap of (AUROC_a - AUROC_b). Returns (mean_gap, lo, hi, n_effective).

    n_effective is the count of non-degenerate resamples (those that contained both
    classes, so AUROC is defined). If n_effective << n_iterations, CI honesty is
    degraded -- caller should surface this.
    """
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    n = len(y_true)
    gaps: list[float] = []
    for _ in range(n_iterations):
        idx = rng.integers(0, n, size=n)
        ys = y_true[idx]
        if len(set(ys.tolist())) < 2:
            continue
        try:
            a = roc_auc_score(ys, scores_a[idx])
            b = roc_auc_score(ys, scores_b[idx])
            gaps.append(float(a - b))
        except ValueError:
            continue
    if not gaps:
        return float("nan"), float("nan"), float("nan"), 0
    arr = np.array(gaps)
    alpha = (1 - ci) / 2
    lo = float(np.quantile(arr, alpha))
    hi = float(np.quantile(arr, 1 - alpha))
    return float(arr.mean()), lo, hi, len(gaps)


# --- verdict + decision layer --------------------------------------------------------

def verdict_label(gap_pp: float, threshold_pp: float = GATE_THRESHOLD_PP) -> str:
    """Pure function of point gap. Frozen 3-bucket semantics per /brainstorm Round 1.

    The decision-layer (`stage2_action`) is computed separately in `decide_stage2_action`.
    """
    if gap_pp >= 5.0:
        return f"CLEAN PASS (gap {gap_pp:+.1f} pp >= 5 pp)"
    if gap_pp >= threshold_pp:
        return f"NOISY PASS (gap {gap_pp:+.1f} pp in [{threshold_pp:.0f}, 5) pp -- diagnostics flagged)"
    return f"FAIL (gap {gap_pp:+.1f} pp < {threshold_pp:.0f} pp threshold)"


def _verdict_bucket(gap_pp: float, threshold_pp: float = GATE_THRESHOLD_PP) -> str:
    """Internal: short bucket label for stage2_action lookup."""
    if gap_pp >= 5.0:
        return "CLEAN"
    if gap_pp >= threshold_pp:
        return "NOISY"
    return "FAIL"


def decide_stage2_action(
    verdict_bucket: str,
    ci_lo_pp: float,
    fusion_outperforms_primary: bool = False,
) -> Stage2Action:
    """Decision-layer: maps (verdict_bucket, ci_lo_pp) to one of 4 actions.

    Per `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md` Verdict-Time
    Pre-Commitments table (locked 2026-05-14). `fusion_outperforms_primary` is
    informational only -- never alters the action.
    """
    if verdict_bucket == "CLEAN":
        return "BURST_STAGE_2"  # ci_lo > 0 or <=0 -- both -> BURST (annotation differs in packet)
    if verdict_bucket == "NOISY":
        if ci_lo_pp > 0:
            return "BURST_STAGE_2"
        return "HOLD_STAGE_2_CI_DEGENERATE"
    return "ALTERNATIVE_POOLING_RERUN"  # FAIL


# --- gate outcome --------------------------------------------------------------------

def compute_gate_outcome(results: list[VariantResult]) -> dict:
    """Pure function: variant results -> gate outcome dict.

    Computes paired bootstrap CI, verdict bucket (frozen point-gap function),
    stage2_action (decision-layer output), and fusion alignment status.

    Raises:
        ValueError: if NT-best and k-mer-XGB variants have mismatched strain_ids.
            Gate-bearing; must fail loudly so misaligned comparison can't be hidden.
        ValueError: if any VariantResult has internal length inconsistency
            (len(strain_ids) != len(per_strain_scores) != len(per_strain_true)).
            Defensive guard per /brainstorm 2026-05-14 follow-up; can't trigger
            from current code paths but pins the contract before any downstream
            AUROC / bootstrap computation trusts the arrays.

    Fusion alignment mismatch: NOT a raise -- the fusion variant is diagnostic only,
    so a mismatch suppresses the fusion-outperforms note but lets the gate proceed.
    """
    # Per-variant length-consistency invariant. Cheap; pins the VariantResult contract.
    for r in results:
        if not (len(r.strain_ids) == len(r.per_strain_scores) == len(r.per_strain_true)):
            raise ValueError(
                f"VariantResult {r.name!r} internal length mismatch: "
                f"strain_ids={len(r.strain_ids)}, scores={len(r.per_strain_scores)}, "
                f"true={len(r.per_strain_true)}"
            )

    def find(name: str) -> VariantResult | None:
        return next((r for r in results if r.name == name), None)

    nt_xgb = find("NT-XGBoost")
    nt_lr = find("NT-logreg")
    kmer = find("k-mer-XGB")
    fusion = find("NT+k-mer-fusion-logreg")

    if nt_xgb is None or nt_lr is None or kmer is None:
        raise ValueError(
            "compute_gate_outcome requires NT-XGBoost, NT-logreg, k-mer-XGB variants"
        )

    # Best NT-only head
    if nt_xgb.auroc >= nt_lr.auroc:
        nt_best_name = "NT-XGBoost"
        nt_best_auroc = nt_xgb.auroc
        nt_best_scores = nt_xgb.per_strain_scores
        nt_best_strain_ids = nt_xgb.strain_ids
    else:
        nt_best_name = "NT-logreg"
        nt_best_auroc = nt_lr.auroc
        nt_best_scores = nt_lr.per_strain_scores
        nt_best_strain_ids = nt_lr.strain_ids

    # Gate-bearing alignment: NT-best vs k-mer-XGB. MUST match.
    if nt_best_strain_ids != kmer.strain_ids:
        raise ValueError(
            "Stage 1 alignment: NT-best and k-mer-XGB strain_ids diverge "
            f"(NT-best n={len(nt_best_strain_ids)}, k-mer n={len(kmer.strain_ids)}); "
            f"refactor regression -- the order-explicit shared API contract was violated."
        )

    gap_pp = (nt_best_auroc - kmer.auroc) * 100
    bucket = _verdict_bucket(gap_pp)
    verdict_bucket_label = verdict_label(gap_pp)

    # Paired bootstrap (uses k-mer's per_strain_true; identical to nt_best's by alignment check above)
    mean_gap, lo, hi, n_eff = paired_bootstrap_ci(
        y_true=kmer.per_strain_true,
        scores_a=nt_best_scores,
        scores_b=kmer.per_strain_scores,
    )

    # Fusion: diagnostic-only alignment check (permissive)
    fusion_alignment_valid = True
    fusion_outperforms_primary = False
    fusion_note = ""
    if fusion is not None:
        if fusion.strain_ids != nt_best_strain_ids:
            fusion_alignment_valid = False
            fusion_note = "fusion alignment mismatch -- diagnostic suppressed"
        else:
            primary_best = max(nt_xgb.auroc, nt_lr.auroc)
            fusion_diff_pp = (fusion.auroc - primary_best) * 100
            if fusion_diff_pp >= 3.0:
                fusion_outperforms_primary = True
                fusion_note = (
                    f"Stage 2 architecture note: fusion outperformed both NT-only heads "
                    f"by {fusion_diff_pp:+.1f} pp -- revisit at Stage 2"
                )

    stage2_action = decide_stage2_action(bucket, lo * 100, fusion_outperforms_primary)

    return {
        "nt_best_name": nt_best_name,
        "nt_best_auroc": nt_best_auroc,
        "nt_best_scores": nt_best_scores,
        "kmer_auroc": kmer.auroc,
        "gap_pp": gap_pp,
        "ci_mean_pp": mean_gap * 100,
        "ci_lo_pp": lo * 100,
        "ci_hi_pp": hi * 100,
        "ci_n_effective": n_eff,
        "verdict_bucket": verdict_bucket_label,  # full label string from verdict_label
        "verdict_bucket_short": bucket,           # "CLEAN" / "NOISY" / "FAIL"
        "stage2_action": stage2_action,
        "fusion_note": fusion_note,
        "fusion_outperforms_primary": fusion_outperforms_primary,
        "fusion_alignment_valid": fusion_alignment_valid,
    }


# --- per-MLST breakdown --------------------------------------------------------------

def per_mlst_breakdown(strain_ids: list[str], mlsts: list[str], y: np.ndarray) -> list[dict]:
    from collections import defaultdict

    by_mlst: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for sid, ml, lbl in zip(strain_ids, mlsts, y):
        by_mlst[ml].append((sid, int(lbl)))
    rows = []
    for ml, strains in sorted(by_mlst.items()):
        rows.append({
            "mlst": ml,
            "n": len(strains),
            "r": sum(1 for _, lbl in strains if lbl == 1),
            "s": sum(1 for _, lbl in strains if lbl == 0),
        })
    return rows


# --- result packet -------------------------------------------------------------------

def write_packet(
    results: list[VariantResult],
    strain_ids: list[str],
    mlsts: list[str],
    y: np.ndarray,
    output_path: Path,
    cohort_path: Path,
    drug: str,
) -> dict:
    outcome = compute_gate_outcome(results)

    today = date.today().isoformat()
    lines = [
        f"# Stage 1 -- N={len(y)} cipro engineering screen ({today})",
        "",
        "> **Engineering screen, NOT a powered statistical comparison.** Stage 1's role is go/no-go for spending Stage 2 N=150 Databricks burst budget. Stage 2 is the real ship gate (>=5 pp AUROC + biology check on gyrA/parC/parE attribution).",
        "",
        f"**Cohort:** `{cohort_path}` (effective N={len(y)}; balance {int((y==1).sum())}R/{int((y==0).sum())}S)",
        f"**Drug:** {drug}",
        f"**Gate threshold:** >={GATE_THRESHOLD_PP:.0f} pp (max(NT-XGBoost, NT-logreg) AUROC - k-mer-XGB AUROC)",
        f"**Best NT-only head:** {outcome['nt_best_name']} (AUROC {outcome['nt_best_auroc']:.3f})",
        f"**Gap vs k-mer-XGB:** {outcome['gap_pp']:+.1f} pp",
        f"**Paired bootstrap 95% CI on gap:** [{outcome['ci_lo_pp']:+.1f}, {outcome['ci_hi_pp']:+.1f}] pp "
        f"(B={BOOTSTRAP_ITERATIONS} effective {outcome['ci_n_effective']}, mean {outcome['ci_mean_pp']:+.1f} pp)",
        f"**Verdict:** {outcome['verdict_bucket']}",
        f"**Stage 2 action:** `{outcome['stage2_action']}`",
        "",
        "## Per-variant LOSO results",
        "",
        "| Variant | AUROC | AUPRC | Gate-bearing? |",
        "|---|---:|---:|:---:|",
    ]
    for r in results:
        gb = "yes" if r.is_gate_bearing else "diagnostic"
        lines.append(f"| {r.name} | {r.auroc:.3f} | {r.auprc:.3f} | {gb} |")

    lines.extend(["", "## Gate analysis", ""])
    lines.append(f"- best NT-only: **{outcome['nt_best_name']}** at **{outcome['nt_best_auroc']:.3f}**")
    lines.append(f"- k-mer-XGB AUROC: **{outcome['kmer_auroc']:.3f}**")
    lines.append(f"- Gap: **{outcome['gap_pp']:+.1f} pp** (best NT-only - k-mer-XGB)")
    lines.append(
        f"- Paired bootstrap 95% CI on gap: **[{outcome['ci_lo_pp']:+.1f}, {outcome['ci_hi_pp']:+.1f}] pp** "
        f"(mean {outcome['ci_mean_pp']:+.1f} pp, B={BOOTSTRAP_ITERATIONS} effective {outcome['ci_n_effective']})"
    )
    if outcome['ci_n_effective'] < 800:
        lines.append(f"- CI honesty degraded ({outcome['ci_n_effective']}/{BOOTSTRAP_ITERATIONS} effective resamples); investigate cohort imbalance.")
    lines.append(f"- Verdict (point-gap function): **{outcome['verdict_bucket']}**")
    lines.append(f"- Stage 2 action (decision-layer): **`{outcome['stage2_action']}`**")
    if outcome['fusion_note']:
        lines.append(f"- {outcome['fusion_note']}")

    # Lineage diagnostic
    mlst_rows = per_mlst_breakdown(strain_ids, mlsts, y)
    unique_mlst = len(mlst_rows)
    largest = max(mlst_rows, key=lambda r: r["n"])
    lines.extend([
        "",
        "## Lineage diagnostic",
        "",
        f"- Unique MLSTs: **{unique_mlst}** of {len(strain_ids)} strains (uniqueness fraction {unique_mlst/len(strain_ids):.2%})",
        f"- Largest MLST group: `{largest['mlst']}` with N={largest['n']} ({largest['r']}R/{largest['s']}S)",
        f"- LOMO note: at N={len(strain_ids)} with this MLST cardinality, most LOMO folds are size-1 -> degenerate. Reporting LOSO only; per-strain table below substitutes for LOMO diagnostics.",
        "",
        "### Per-strain LOSO predictions",
        "",
        "| Strain | MLST | True | NT-best score | k-mer score | NT-best correct? | k-mer correct? |",
        "|---|---|:---:|---:|---:|:---:|:---:|",
    ])
    kmer = next(r for r in results if r.name == "k-mer-XGB")
    for sid, ml, t, ntb, km in zip(
        strain_ids, mlsts, y, outcome["nt_best_scores"], kmer.per_strain_scores
    ):
        nt_pred = 1 if ntb >= 0.5 else 0
        km_pred = 1 if km >= 0.5 else 0
        lines.append(
            f"| {sid} | {ml} | {int(t)} | {ntb:.3f} | {km:.3f} | "
            f"{'OK' if nt_pred == t else 'X'} | "
            f"{'OK' if km_pred == t else 'X'} |"
        )

    lines.extend([
        "",
        "## Notes",
        "",
        "- All variants ran with `calibrate=False` (uniform calibration discipline; calibration is small-N footgun per LESSONS_LEARNED 2026-05-14).",
        "- k-mer + fusion use within-fold vocab rebuild from training-set sequences only (no held-out leakage).",
        "- Gene-presence variant NOT included -- RefSeq GFF3 carries `gene=` for ~11% of CDSs -> INDETERMINATE_IDENTIFIER_OOV on this annotation source. See `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md`.",
        "- AMRFinderPlus POINT* SNP-table baseline NOT included -- deferred to Stage 2 per Phase2_Decision_Gate D6. 'Best classical' here is bounded by what was run; gyrA/parC/parE point-mutation features are NOT part of the comparator.",
        "- LOSO at N=38 has +-0.10 noise floor on AUROC; >=3 pp is INSIDE the noise. The bootstrap CI surfaces this honestly.",
        "- Verdict semantics are frozen as a pure function of point gap. The `stage2_action` field is the decision-layer; see `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md` Verdict-Time Pre-Commitments.",
        "",
        "## Next action by stage2_action",
        "",
        "- **`BURST_STAGE_2`** -> proceed to Stage 2 Databricks burst with N=150 cohort build.",
        "- **`HOLD_STAGE_2_CI_DEGENERATE`** -> do NOT spend Stage 2 burst budget; next escalation is `ALTERNATIVE_POOLING_RERUN`.",
        "- **`ALTERNATIVE_POOLING_RERUN`** -> Stage 1b with `mean+max` aggregation; if still <3 pp, escalate to `PIVOT_TO_BAKTA`.",
        "- **`PIVOT_TO_BAKTA`** -> Bakta re-annotation + gene-presence comparator pathway per `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md` follow-up.",
    ])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return {
        "verdict": outcome["verdict_bucket"],
        "stage2_action": outcome["stage2_action"],
        "gap_pp": outcome["gap_pp"],
        "nt_best_name": outcome["nt_best_name"],
        "nt_best_auroc": outcome["nt_best_auroc"],
        "kmer_auroc": outcome["kmer_auroc"],
        "ci_lo_pp": outcome["ci_lo_pp"],
        "ci_hi_pp": outcome["ci_hi_pp"],
    }


# --- main ----------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 1 N=40 cipro engineering screen")
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_n40_cipro_cohort.parquet"))
    parser.add_argument("--nt-cache", type=Path, default=Path("D:/dna_decode_cache/embeddings/nt_n40_cipro.h5"))
    parser.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--kmer-k", type=int, default=8)
    parser.add_argument("--kmer-top-n", type=int, default=10_000)
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/stage1_n40_cipro_{date.today().isoformat()}.md")

    cohort = load_cohort(args.cohort)
    print(f"[stage1] cohort: {len(cohort.strains)} strains; drug={args.drug}")

    X_nt, seqs_by_strain, labels_by_strain, strain_ids, mlsts = load_features(
        cohort, args.nt_cache, args.refseq_cache, args.drug
    )
    y = np.array([labels_by_strain[s] for s in strain_ids], dtype=int)
    print(f"[stage1] effective N={len(y)}; balance {int((y==1).sum())}R/{int((y==0).sum())}S; NT shape={X_nt.shape}")

    results: list[VariantResult] = []

    print("[stage1] Running NT-XGBoost LOSO...")
    cv = leave_one_strain_out_cv(X_nt, y, strain_ids, _nt_xgb_train, _nt_xgb_predict, drug="cipro_stage1")
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult(
        name="NT-XGBoost",
        auroc=float(m.auroc),
        auprc=float(m.auprc),
        per_strain_scores=cv.all_y_score,
        per_strain_true=cv.all_y_true,
        strain_ids=cv.strain_ids,
        is_gate_bearing=True,
    ))
    print(f"  AUROC={m.auroc:.3f}")

    print("[stage1] Running NT-logreg LOSO...")
    cv = leave_one_strain_out_cv(X_nt, y, strain_ids, _nt_logreg_train, _nt_logreg_predict, drug="cipro_stage1")
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult(
        name="NT-logreg",
        auroc=float(m.auroc),
        auprc=float(m.auprc),
        per_strain_scores=cv.all_y_score,
        per_strain_true=cv.all_y_true,
        strain_ids=cv.strain_ids,
        is_gate_bearing=True,
    ))
    print(f"  AUROC={m.auroc:.3f}")

    print("[stage1] Running k-mer-XGB LOSO...")
    cv = run_kmer_xgboost_loso(seqs_by_strain, labels_by_strain, strain_ids, drug="cipro_stage1", k=args.kmer_k, top_n=args.kmer_top_n)
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult(
        name="k-mer-XGB",
        auroc=float(m.auroc),
        auprc=float(m.auprc),
        per_strain_scores=cv.all_y_score,
        per_strain_true=cv.all_y_true,
        strain_ids=cv.strain_ids,
        is_gate_bearing=True,
    ))
    print(f"  AUROC={m.auroc:.3f}")

    print("[stage1] Running NT+k-mer-fusion-logreg LOSO (diagnostic only)...")
    cv = run_fusion_loso(X_nt, seqs_by_strain, labels_by_strain, strain_ids, drug="cipro_stage1", k=args.kmer_k, top_n=args.kmer_top_n)
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    results.append(VariantResult(
        name="NT+k-mer-fusion-logreg",
        auroc=float(m.auroc),
        auprc=float(m.auprc),
        per_strain_scores=cv.all_y_score,
        per_strain_true=cv.all_y_true,
        strain_ids=cv.strain_ids,
        is_gate_bearing=False,
    ))
    print(f"  AUROC={m.auroc:.3f} (diagnostic; does NOT count for gate)")

    print(f"[stage1] Writing result packet to {args.output}")
    summary = write_packet(results, strain_ids, mlsts, y, args.output, args.cohort, args.drug)
    print(f"[stage1] verdict: {summary['verdict']}")
    print(f"[stage1] stage2_action: {summary['stage2_action']}")
    print(
        f"[stage1] best NT-only: {summary['nt_best_name']} @ {summary['nt_best_auroc']:.3f}; "
        f"k-mer-XGB @ {summary['kmer_auroc']:.3f}; gap {summary['gap_pp']:+.1f} pp; "
        f"CI [{summary['ci_lo_pp']:+.1f}, {summary['ci_hi_pp']:+.1f}] pp"
    )
    return 0 if summary["stage2_action"] == "BURST_STAGE_2" else 1


if __name__ == "__main__":
    sys.exit(main())
