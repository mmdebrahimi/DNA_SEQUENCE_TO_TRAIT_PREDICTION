"""Step 11.5 — Quantization-fidelity check (selective Phase 1 addition).

Phase 1 plan assumed 4-bit Evo quantization via bitsandbytes on RTX 4090.
Actual hardware is GTX 860M (CC=5.0) → bitsandbytes structurally unavailable
(requires CC ≥ 7.0), so this fidelity check is unreachable on the project's
real machine. Kept as scaffolding for future compute-upgrade scenarios. If
quantization were available + distorted attribution maps materially, every
Phase 1 attribution-precision number would be quantization-conditional and
the headline results would become unreliable.

This script runs gene-level ISM on a small subset (5-10 strains) at BOTH
4-bit and full-precision; computes concordance via top-K=20 set
intersection + Spearman rank correlation on prediction-delta values.

One-time cost (~1 rented A100 hour). Output:
`data/processed/quantize_fidelity_report.md` with GO/NO-GO recommendation.

This script is a SCAFFOLD — it sets up the comparison workflow. Real
invocation requires:
  1. Two trained classifiers (4-bit + full-precision) sharing the cohort.
  2. Cached embeddings at both precisions (separate HDF5 files).
  3. ISM run via existing dna_decode.interp.mutagenesis pipeline.

The threshold for "fidelity GO" is ≥0.7 Spearman concordance on top-K=20
attribution-delta rank. Below that → Phase 1 results are flagged
quantization-conditional in documentation; rerun at full precision for
the headline numbers.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


SPEARMAN_GO_THRESHOLD = 0.7
TOP_K_INTERSECTION_GO_THRESHOLD = 0.6  # ≥12 of top-20 must overlap


@dataclass
class FidelityResult:
    """Per-strain (or aggregated) fidelity comparison."""

    strain_id: str
    drug: str
    top_k: int
    intersection_size: int
    spearman: float | None
    fidelity_score: float  # 0..1; combines intersection-fraction + Spearman


@dataclass
class FidelityReport:
    """Aggregate report across strains."""

    per_strain: list[FidelityResult]
    overall_intersection_mean: float
    overall_spearman_mean: float
    go_no_go: str
    failure_reasons: list[str]


def compare_top_k_attribution(
    full_precision_table: pd.DataFrame,
    quantized_table: pd.DataFrame,
    top_k: int = 20,
) -> tuple[int, float | None]:
    """Compare two GeneEffectTables — top-K intersection + Spearman.

    Args:
        full_precision_table: GeneEffectTable from full-precision run.
        quantized_table: GeneEffectTable from 4-bit-quantized run.
        top_k: number of top genes to compare.

    Returns:
        (intersection_size, spearman_correlation_or_None)
    """
    if len(full_precision_table) == 0 or len(quantized_table) == 0:
        return 0, None

    full_top = set(full_precision_table.head(top_k)["gene_id"].tolist())
    quant_top = set(quantized_table.head(top_k)["gene_id"].tolist())
    intersection = full_top & quant_top

    # Spearman rank correlation across the union of top-K genes
    union_genes = sorted(full_top | quant_top)
    if len(union_genes) < 2:
        return len(intersection), None

    # Map each gene to its rank (1-based; rank for absent gene = top_k + 1)
    def _rank_lookup(table: pd.DataFrame) -> dict[str, int]:
        out: dict[str, int] = {}
        for i, gene in enumerate(table["gene_id"].head(top_k).tolist()):
            out[gene] = i + 1
        return out

    full_ranks = _rank_lookup(full_precision_table)
    quant_ranks = _rank_lookup(quantized_table)

    full_vector = [full_ranks.get(g, top_k + 1) for g in union_genes]
    quant_vector = [quant_ranks.get(g, top_k + 1) for g in union_genes]

    spearman: float | None
    try:
        from scipy.stats import spearmanr
        rho, _ = spearmanr(full_vector, quant_vector)
        spearman = float(rho) if not np.isnan(rho) else None
    except ImportError:
        # Pure-numpy fallback: Spearman is Pearson on ranks. Ranks computed via
        # argsort + reverse-permutation; ties get distinct integer ranks (good
        # enough for our discrete top-K-rank inputs).
        full_arr = np.asarray(full_vector, dtype=float)
        quant_arr = np.asarray(quant_vector, dtype=float)
        full_ranks = full_arr.argsort().argsort().astype(float)
        quant_ranks = quant_arr.argsort().argsort().astype(float)
        if full_ranks.std() == 0 or quant_ranks.std() == 0:
            spearman = None
        else:
            rho = float(np.corrcoef(full_ranks, quant_ranks)[0, 1])
            spearman = rho if not np.isnan(rho) else None

    return len(intersection), spearman


def build_fidelity_report(
    per_strain_results: list[FidelityResult],
    top_k: int,
) -> FidelityReport:
    """Aggregate per-strain results + decide GO/NO-GO."""
    if not per_strain_results:
        return FidelityReport(
            per_strain=[],
            overall_intersection_mean=0.0,
            overall_spearman_mean=0.0,
            go_no_go="NO-GO",
            failure_reasons=["no per-strain results to aggregate"],
        )

    intersection_fractions = [r.intersection_size / top_k for r in per_strain_results]
    spearmans = [r.spearman for r in per_strain_results if r.spearman is not None]

    mean_intersection = float(np.mean(intersection_fractions))
    mean_spearman = float(np.mean(spearmans)) if spearmans else 0.0

    failure_reasons: list[str] = []
    if mean_intersection < TOP_K_INTERSECTION_GO_THRESHOLD:
        failure_reasons.append(
            f"mean top-{top_k} intersection-fraction {mean_intersection:.2f} "
            f"< threshold {TOP_K_INTERSECTION_GO_THRESHOLD}"
        )
    if mean_spearman < SPEARMAN_GO_THRESHOLD:
        failure_reasons.append(
            f"mean Spearman {mean_spearman:.2f} < threshold {SPEARMAN_GO_THRESHOLD}"
        )

    return FidelityReport(
        per_strain=per_strain_results,
        overall_intersection_mean=mean_intersection,
        overall_spearman_mean=mean_spearman,
        go_no_go="GO" if not failure_reasons else "NO-GO",
        failure_reasons=failure_reasons,
    )


def write_report(report: FidelityReport, output_path: Path, drug: str) -> Path:
    """Write the fidelity report as markdown."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "# Quantization-Fidelity Check",
        "",
        f"**Drug:** {drug}",
        f"**Verdict:** `{report.go_no_go}`",
        f"**Mean top-K intersection:** {report.overall_intersection_mean:.2f}",
        f"**Mean Spearman:** {report.overall_spearman_mean:.2f}",
        "",
        "## Thresholds",
        f"- top-K intersection fraction ≥ {TOP_K_INTERSECTION_GO_THRESHOLD}",
        f"- Spearman rank correlation ≥ {SPEARMAN_GO_THRESHOLD}",
        "",
        "## Per-strain results",
        "",
        f"| Strain | top-K intersect | Spearman | Score |",
        "|---|---|---|---|",
    ]
    for r in report.per_strain:
        spearman_str = f"{r.spearman:.2f}" if r.spearman is not None else "—"
        lines.append(
            f"| {r.strain_id} | {r.intersection_size}/{r.top_k} | {spearman_str} | "
            f"{r.fidelity_score:.2f} |"
        )

    if report.failure_reasons:
        lines += ["", "## Failure reasons", ""]
        for reason in report.failure_reasons:
            lines.append(f"- {reason}")
        lines += [
            "",
            "**Remediation:** Re-run Phase 1 headline numbers at full precision on a "
            "rented A100. Document attribution numbers as quantization-conditional in "
            "the Phase 1 results section.",
        ]

    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase 1 quantization-fidelity check (Step 11.5)"
    )
    parser.add_argument(
        "--full-precision-attributions",
        required=True,
        help="JSON file: list of {strain_id, table_csv_path} for full-precision runs.",
    )
    parser.add_argument(
        "--quantized-attributions",
        required=True,
        help="JSON file: list of {strain_id, table_csv_path} for 4-bit runs.",
    )
    parser.add_argument("--drug", required=True)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument(
        "--output",
        default="data/processed/quantize_fidelity_report.md",
        help="Output markdown path.",
    )
    args = parser.parse_args(argv)

    full_manifest = json.loads(Path(args.full_precision_attributions).read_text())
    quant_manifest = json.loads(Path(args.quantized_attributions).read_text())

    # Index by strain_id
    full_by_strain = {entry["strain_id"]: entry for entry in full_manifest}
    quant_by_strain = {entry["strain_id"]: entry for entry in quant_manifest}
    common_strains = sorted(full_by_strain.keys() & quant_by_strain.keys())
    if not common_strains:
        print("[quant-fidelity] no overlapping strains", file=sys.stderr)
        return 2

    per_strain: list[FidelityResult] = []
    for sid in common_strains:
        full_table = pd.read_csv(full_by_strain[sid]["table_csv_path"], sep="\t")
        quant_table = pd.read_csv(quant_by_strain[sid]["table_csv_path"], sep="\t")
        intersection_size, spearman = compare_top_k_attribution(
            full_table, quant_table, top_k=args.top_k
        )
        intersection_frac = intersection_size / args.top_k
        spearman_score = spearman if spearman is not None else 0.0
        # Composite score: average of intersection-fraction and Spearman
        fidelity_score = (intersection_frac + max(spearman_score, 0.0)) / 2.0
        per_strain.append(
            FidelityResult(
                strain_id=sid,
                drug=args.drug,
                top_k=args.top_k,
                intersection_size=intersection_size,
                spearman=spearman,
                fidelity_score=fidelity_score,
            )
        )

    report = build_fidelity_report(per_strain, top_k=args.top_k)
    write_report(report, Path(args.output), drug=args.drug)
    print(f"[quant-fidelity] verdict: {report.go_no_go}")
    print(f"[quant-fidelity] mean intersection: {report.overall_intersection_mean:.2f}")
    print(f"[quant-fidelity] mean Spearman: {report.overall_spearman_mean:.2f}")
    if report.failure_reasons:
        for reason in report.failure_reasons:
            print(f"[quant-fidelity]   - {reason}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
