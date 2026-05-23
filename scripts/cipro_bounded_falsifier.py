"""Cipro bounded-falsifier runner — Claude DRAFT for Codex review.

Authoring discipline (per `wiki/cipro_bounded_falsifier_coordination_plan_2026-05-22.md`
Section 1 role-split):

- Codex (Precision 7780) has authority on RUNNER MECHANICS (model invocation,
  ISM call ordering, CV-fold reuse, exception handling).
- Claude (this file, GTX 860M laptop) has authority on the DIAGNOSTIC EXPORTS
  SCHEMA + the per-bucket verdict logic + the subset-JSON contract.

This draft exists so Codex can diff against a working skeleton rather than
implement from scratch. Either adopt with edits, or replace; do NOT preserve
this verbatim if your runtime / model-loading idioms differ.

Inputs:
  --cohort       parquet cohort path (default: stage2_n150_cipro_cohort.parquet)
  --model        trained NT-XGBoost pickle
  --cache        HDF5 NT embedding cache
  --refseq-cache RefSeq local cache (GFF3 + FASTA per accession)
  --subset       wiki/cipro_bounded_falsifier_subset_2026-05-22.json
  --leakage-check-json  optional gating; if loso_leakage_present=True, abort

Outputs:
  wiki/cipro_bounded_falsifier_results_<DATE>.json   (machine-readable)
  wiki/cipro_bounded_falsifier_results_<DATE>.md     (narrative)

Verdict matrix (per coordination plan Section 2 + subset JSON `verdict_matrix`):
  Bucket A ERS pass  + Bucket B ELX pass  + Bucket C handled  -> PASS
  Bucket A pass      + Bucket B fail      + any              -> FAIL
  Bucket A fail                                              -> RUNNER_REGRESSION
  Method change breaks Bucket A's positive deltas            -> REVERT
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date as _date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from dna_decode.data.annotations import extract_cds_sequences, parse_gff3
from dna_decode.data.cohort import load_cohort
from dna_decode.data.mic_tiers import (
    classify_gene_symbol,
    loci_by_mechanism_for,
    primary_mechanisms_for,
)
from dna_decode.data.refseq import fasta_path, gff_path
from dna_decode.interp.mutagenesis import gene_level_mutagenesis
from dna_decode.models.cache import EmbeddingCache
from dna_decode.models.classifiers import (
    aggregate_strain_features,
    predict_proba,
)


# ---------------------------------------------------------------------------
# Diagnostic-exports schema (Claude-authored; coordination plan Section 4).
# ---------------------------------------------------------------------------

SATURATION_PROBA_THRESHOLD = 0.95
SATURATION_MAX_ABS_DELTA_THRESHOLD = 0.01


@dataclass
class StrainResult:
    """Per-strain diagnostic record. JSON schema is stable across Codex runs."""

    strain_id: str
    accession: str
    label: str
    bucket: str  # "A_ERS" / "B_ELX" / "C_NEGATIVE"
    n_cached_genes: int
    baseline_proba_R: float
    baseline_logit: float
    max_abs_delta_all_genes: float
    saturation_flag: bool
    n_known_loci_hits: int
    best_known_locus_rank_pos_delta: int | None  # NEW positive-only ranking
    best_known_locus_rank_abs_delta: int | None  # baseline (existing audit logic)
    per_known_locus: list[dict[str, Any]] = field(default_factory=list)
    top_10_positive_delta: list[dict[str, Any]] = field(default_factory=list)
    indeterminate_reason: str | None = None  # e.g. "ALL_NEGATIVE_DELTA"


# ---------------------------------------------------------------------------
# Per-bucket pass criteria (coordination plan Section 2 + subset JSON).
# ---------------------------------------------------------------------------


def bucket_A_pass(results: list[StrainResult]) -> bool:
    """ERS control: >=3 of 4 still rank a QRDR locus in top-10 positive-only delta."""
    qrdr = {"gyrA", "gyrB", "parC", "parE"}
    n_pass = 0
    for r in results:
        if any(
            h.get("alias") in qrdr and h.get("rank_pos_delta") is not None and h["rank_pos_delta"] <= 10
            for h in r.per_known_locus
        ):
            n_pass += 1
    return n_pass >= 3


def bucket_B_pass(results: list[StrainResult]) -> bool:
    """ELX failure: >=2 of 4 move into top-10 AND median rank shift >=100x."""
    n_in_top10 = 0
    rank_ratios: list[float] = []
    for r in results:
        if r.best_known_locus_rank_pos_delta is not None and r.best_known_locus_rank_pos_delta <= 10:
            n_in_top10 += 1
        if (
            r.best_known_locus_rank_abs_delta is not None
            and r.best_known_locus_rank_pos_delta is not None
            and r.best_known_locus_rank_pos_delta > 0
        ):
            rank_ratios.append(r.best_known_locus_rank_abs_delta / r.best_known_locus_rank_pos_delta)
    if n_in_top10 < 2:
        return False
    median_ratio = float(np.median(rank_ratios)) if rank_ratios else 0.0
    return median_ratio >= 100.0


def bucket_C_handled(results: list[StrainResult]) -> tuple[bool, str]:
    """Bucket C: either recover a hit into top-50 OR explicitly INDETERMINATE w/ saturation diagnostic.

    Returns (handled, descriptor).
    """
    recovered = sum(
        1
        for r in results
        if r.best_known_locus_rank_pos_delta is not None and r.best_known_locus_rank_pos_delta <= 50
    )
    indeterminate_saturated = sum(
        1 for r in results if r.indeterminate_reason == "ALL_NEGATIVE_DELTA" and r.saturation_flag
    )
    if recovered >= 1:
        return True, f"{recovered}/{len(results)} recovered into top-50"
    if indeterminate_saturated == len(results):
        return True, f"all {len(results)} flagged INDETERMINATE_ALL_NEGATIVE_DELTA with saturation_flag=True"
    return False, f"{recovered} recovered; {indeterminate_saturated} indeterminate-saturated of {len(results)}"


def compute_verdict(
    bucket_A: list[StrainResult],
    bucket_B: list[StrainResult],
    bucket_C: list[StrainResult],
) -> tuple[str, dict[str, Any]]:
    """Apply the verdict matrix from coordination plan Section 2."""
    A = bucket_A_pass(bucket_A)
    B = bucket_B_pass(bucket_B)
    C_ok, C_desc = bucket_C_handled(bucket_C)
    if not A:
        verdict = "RUNNER_REGRESSION"
        rationale = "Bucket A (ERS control) failed; positive-only Δ ranking should not break working strains."
    elif not B:
        verdict = "FAIL"
        rationale = "Bucket A passed but Bucket B (ELX failure cases) did not improve; ranking is not the bottleneck."
    elif not C_ok:
        verdict = "REVERT"
        rationale = f"Bucket A + B passed but Bucket C unresolved: {C_desc}; method change ambiguous on negative-delta strains."
    else:
        verdict = "PASS"
        rationale = "All 3 buckets pass; continue method refinement on full N=67 + Mash-cluster."
    details = {
        "bucket_A_pass": A,
        "bucket_B_pass": B,
        "bucket_C_handled": C_ok,
        "bucket_C_descriptor": C_desc,
        "rationale": rationale,
    }
    return verdict, details


# ---------------------------------------------------------------------------
# Per-strain falsifier execution.
# ---------------------------------------------------------------------------


def _logit(p: float) -> float:
    eps = 1e-9
    p_clamped = min(max(p, eps), 1.0 - eps)
    return math.log(p_clamped / (1.0 - p_clamped))


def _ranked_by(deltas: list[float], descending: bool = True) -> list[int]:
    """Return rank assignments (1-indexed) for the given deltas. Larger = better when descending."""
    order = np.argsort(-np.array(deltas) if descending else np.array(deltas))
    ranks = np.empty_like(order)
    for rank, idx in enumerate(order, start=1):
        ranks[idx] = rank
    return ranks.tolist()


def run_strain(
    strain_id: str,
    accession: str,
    label: str,
    bucket: str,
    classifier,
    cache: EmbeddingCache,
    refseq_cache: Path,
    drug: str,
) -> StrainResult:
    """Run gene-level ISM + diagnostic exports for one strain."""
    gene_ids = cache.list_genes(strain_id)
    if not gene_ids:
        return StrainResult(
            strain_id=strain_id,
            accession=accession,
            label=label,
            bucket=bucket,
            n_cached_genes=0,
            baseline_proba_R=float("nan"),
            baseline_logit=float("nan"),
            max_abs_delta_all_genes=0.0,
            saturation_flag=False,
            n_known_loci_hits=0,
            best_known_locus_rank_pos_delta=None,
            best_known_locus_rank_abs_delta=None,
            indeterminate_reason="MISSING_FROM_CACHE",
        )

    gene_matrix = cache.bulk_get([(strain_id, g) for g in gene_ids])
    gene_emb = {g: gene_matrix[i] for i, g in enumerate(gene_ids)}

    # baseline diagnostic exports
    baseline_feat = aggregate_strain_features(gene_matrix, "mean").reshape(1, -1)
    baseline_p = float(predict_proba(classifier, baseline_feat)[0])
    baseline_l = _logit(baseline_p)

    # symbol map for known-locus matching
    gff = gff_path(accession, refseq_cache)
    annotations = parse_gff3(gff) if gff.exists() else None
    symbol_map: dict[str, str] = {}
    if annotations is not None and "gene_symbol" in annotations.columns:
        for _, row in annotations.iterrows():
            gid = str(row.get("gene_id", "") or "")
            sym = str(row.get("gene_symbol", "") or "")
            if gid and sym:
                symbol_map[gid] = sym

    # gene-level ISM (this is the per-strain ~95 s call on a 5000-gene strain)
    effects = gene_level_mutagenesis(classifier, gene_emb, annotations=annotations)
    if effects.empty:
        return StrainResult(
            strain_id=strain_id,
            accession=accession,
            label=label,
            bucket=bucket,
            n_cached_genes=len(gene_ids),
            baseline_proba_R=baseline_p,
            baseline_logit=baseline_l,
            max_abs_delta_all_genes=0.0,
            saturation_flag=baseline_p >= SATURATION_PROBA_THRESHOLD,
            n_known_loci_hits=0,
            best_known_locus_rank_pos_delta=None,
            best_known_locus_rank_abs_delta=None,
            indeterminate_reason="ISM_RETURNED_EMPTY",
        )

    deltas = effects["prediction_delta"].astype(float).tolist()
    max_abs = float(max(abs(d) for d in deltas))
    saturated = (baseline_p >= SATURATION_PROBA_THRESHOLD) and (max_abs < SATURATION_MAX_ABS_DELTA_THRESHOLD)

    # ranking #1: abs-Δ descending (baseline / existing audit logic)
    ranks_abs = _ranked_by([abs(d) for d in deltas], descending=True)
    # ranking #2: positive-only Δ descending — negative deltas pushed to bottom
    pos_for_ranking = [d if d > 0 else float("-inf") for d in deltas]
    # Replace -inf with a single rank-tie value AFTER the positives; use 0.0 fallback for ranking only.
    pos_replaced = [d if d > 0 else 0.0 for d in pos_for_ranking]
    ranks_pos = _ranked_by(pos_replaced, descending=True)
    # If a strain has ZERO positive deltas, all positive ranks are tied — flag.
    any_pos = any(d > 0 for d in deltas)

    # match known loci
    known_locus_set = set().union(*loci_by_mechanism_for(drug).values())
    primary_set = primary_mechanisms_for(drug)
    per_known: list[dict[str, Any]] = []
    for i, row in effects.reset_index(drop=True).iterrows():
        gid = str(row.get("gene_id", ""))
        sym = symbol_map.get(gid, "")
        if not sym:
            continue
        mech = classify_gene_symbol(drug, sym)
        if mech in {None, ""}:
            continue
        per_known.append(
            {
                "alias": sym,
                "gene_id": gid,
                "mechanism": mech,
                "is_primary_mechanism": mech in primary_set,
                "prediction_delta": float(row["prediction_delta"]),
                "baseline_probability": float(row.get("baseline_probability", baseline_p)),
                "knockout_probability": float(row.get("knockout_probability", float("nan"))),
                "rank_abs_delta": int(ranks_abs[i]),
                "rank_pos_delta": int(ranks_pos[i]) if row["prediction_delta"] > 0 else None,
            }
        )

    best_abs = (
        min((h["rank_abs_delta"] for h in per_known), default=None) if per_known else None
    )
    best_pos = (
        min((h["rank_pos_delta"] for h in per_known if h["rank_pos_delta"] is not None), default=None)
        if per_known
        else None
    )

    # top-10 positive-Δ aliases (regardless of known-locus status)
    top10_pos: list[dict[str, Any]] = []
    pos_indices = [i for i, d in enumerate(deltas) if d > 0]
    pos_indices.sort(key=lambda i: -deltas[i])
    for idx in pos_indices[:10]:
        gid = str(effects.iloc[idx].get("gene_id", ""))
        top10_pos.append(
            {
                "rank": ranks_pos[idx],
                "gene_id": gid,
                "gene_symbol": symbol_map.get(gid, ""),
                "prediction_delta": deltas[idx],
            }
        )

    indeterminate = None
    if not any_pos and per_known:
        indeterminate = "ALL_NEGATIVE_DELTA"

    return StrainResult(
        strain_id=strain_id,
        accession=accession,
        label=label,
        bucket=bucket,
        n_cached_genes=len(gene_ids),
        baseline_proba_R=baseline_p,
        baseline_logit=baseline_l,
        max_abs_delta_all_genes=max_abs,
        saturation_flag=saturated,
        n_known_loci_hits=len(per_known),
        best_known_locus_rank_pos_delta=best_pos,
        best_known_locus_rank_abs_delta=best_abs,
        per_known_locus=per_known,
        top_10_positive_delta=top10_pos,
        indeterminate_reason=indeterminate,
    )


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--refseq-cache", type=Path, required=True)
    parser.add_argument("--subset", type=Path,
                        default=Path("wiki/cipro_bounded_falsifier_subset_2026-05-22.json"))
    parser.add_argument("--leakage-check-json", type=Path, default=None,
                        help="Optional pre-gate: if loso_leakage_present=True, abort.")
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output-prefix", type=Path,
                        default=Path(f"wiki/cipro_bounded_falsifier_results_{_date.today().isoformat()}"))
    args = parser.parse_args(argv)

    # ---- gating: leakage check ----
    if args.leakage_check_json is not None and args.leakage_check_json.exists():
        leak = json.loads(args.leakage_check_json.read_text())
        if leak.get("loso_leakage_present"):
            sys.stderr.write(
                f"[falsifier] ABORT: leakage check reports same-genome LOSO leakage: {leak.get('strain_ids_sharing_accession')}\n"
                f"[falsifier] DEDUP cohort + retrain before interpreting falsifier; see {args.leakage_check_json}\n"
            )
            return 3

    # ---- load subset ----
    subset = json.loads(args.subset.read_text())
    bucket_specs: dict[str, list[dict[str, str]]] = {
        "A_ERS": [{"strain_id": s["strain_id"], "accession": s["accession"]} for s in subset["bucket_A_control_ERS"]],
        "B_ELX": [{"strain_id": s["strain_id"], "accession": s["accession"]} for s in subset["bucket_B_ELX_failure"]],
        "C_NEGATIVE": [{"strain_id": s["strain_id"], "accession": s["accession"]} for s in subset["bucket_C_all_negative_delta"]],
    }

    # ---- load cohort + classifier + cache ----
    cohort = load_cohort(args.cohort)
    import pickle
    with open(args.model, "rb") as f:
        model_obj = pickle.load(f)
    classifier = model_obj["classifier"] if isinstance(model_obj, dict) and "classifier" in model_obj else model_obj
    cache = EmbeddingCache(
        args.cache,
        model_name="nucleotide_transformer",
        model_version="InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
        embedding_dim=512,
    )

    # ---- run buckets ----
    label_lookup = {s.strain_id: s.ast_labels.get(args.drug.lower()) for s in cohort.strains}
    all_results: list[StrainResult] = []
    per_bucket: dict[str, list[StrainResult]] = {"A_ERS": [], "B_ELX": [], "C_NEGATIVE": []}
    for bucket_name, items in bucket_specs.items():
        for item in items:
            sid = item["strain_id"]
            acc = item["accession"]
            label = "R" if label_lookup.get(sid) == 1 else ("S" if label_lookup.get(sid) == 0 else "?")
            print(f"[falsifier] {bucket_name} {sid} ({acc}) ...")
            r = run_strain(
                sid, acc, label, bucket_name,
                classifier=classifier,
                cache=cache,
                refseq_cache=args.refseq_cache,
                drug=args.drug,
            )
            all_results.append(r)
            per_bucket[bucket_name].append(r)

    # ---- verdict ----
    verdict, details = compute_verdict(
        per_bucket["A_ERS"], per_bucket["B_ELX"], per_bucket["C_NEGATIVE"]
    )

    # ---- emit JSON ----
    out_json = args.output_prefix.with_suffix(".json")
    out_md = args.output_prefix.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run_date": _date.today().isoformat(),
        "drug": args.drug,
        "subset_path": str(args.subset),
        "model_path": str(args.model),
        "cache_path": str(args.cache),
        "cohort_path": str(args.cohort),
        "leakage_check_json": str(args.leakage_check_json) if args.leakage_check_json else None,
        "verdict": verdict,
        "verdict_details": details,
        "saturation_threshold_proba": SATURATION_PROBA_THRESHOLD,
        "saturation_threshold_max_abs_delta": SATURATION_MAX_ABS_DELTA_THRESHOLD,
        "results": [r.__dict__ for r in all_results],
    }
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[falsifier] wrote {out_json}")

    # ---- emit MD ----
    lines = [
        f"# Cipro bounded-falsifier results ({_date.today().isoformat()})",
        "",
        f"**Verdict:** {verdict}",
        f"**Rationale:** {details['rationale']}",
        "",
        f"- Bucket A (ERS control) pass: {details['bucket_A_pass']}",
        f"- Bucket B (ELX failure)  pass: {details['bucket_B_pass']}",
        f"- Bucket C (all-negative) handled: {details['bucket_C_handled']} ({details['bucket_C_descriptor']})",
        "",
        "## Per-strain diagnostic exports",
        "",
        "| bucket | strain_id | accession | label | n_genes | baseline_P(R) | baseline_logit | max_abs_Δ_all | sat_flag | hits | best_rank_abs | best_rank_pos | indeterminate |",
        "|---|---|---|---|---:|---:|---:|---:|:---:|---:|---:|---:|---|",
    ]
    for r in all_results:
        lines.append(
            f"| {r.bucket} | {r.strain_id} | {r.accession} | {r.label} | {r.n_cached_genes} | "
            f"{r.baseline_proba_R:.4f} | {r.baseline_logit:.3f} | {r.max_abs_delta_all_genes:.4f} | "
            f"{'YES' if r.saturation_flag else 'no'} | {r.n_known_loci_hits} | "
            f"{r.best_known_locus_rank_abs_delta} | {r.best_known_locus_rank_pos_delta} | "
            f"{r.indeterminate_reason or ''} |"
        )
    lines += [
        "",
        "## Per-strain known-locus hits (positive-only Δ ranking)",
        "",
    ]
    for r in all_results:
        lines.append(f"### {r.strain_id} ({r.accession}, bucket={r.bucket}, label={r.label})")
        if not r.per_known_locus:
            lines.append(f"  (no known-locus hits; indeterminate={r.indeterminate_reason or 'none'})")
            continue
        lines.append("")
        lines.append("| alias | mechanism | primary? | Δ | rank_abs | rank_pos |")
        lines.append("|---|---|:---:|---:|---:|---:|")
        for h in r.per_known_locus:
            lines.append(
                f"| {h['alias']} | {h['mechanism']} | {'Y' if h['is_primary_mechanism'] else 'n'} | "
                f"{h['prediction_delta']:.5f} | {h['rank_abs_delta']} | {h['rank_pos_delta']} |"
            )
        lines.append("")
    lines += [
        "",
        "## Decision branch",
        "",
        "Per `wiki/cipro_bounded_falsifier_coordination_plan_2026-05-22.md` Section 6:",
        "",
        f"- **{verdict}** → " + {
            "PASS": "diagnostic exports + Mash-cluster N=147 + ship v0 standard.",
            "FAIL": "ship v0 with scope-limit doc (Section 7 template); ranking is not the bottleneck.",
            "RUNNER_REGRESSION": "halt; debug runner against the 2026-05-21 audit baseline before interpreting.",
            "REVERT": "revert method change; ship v0 with current attribution + scope-limit doc.",
        }[verdict],
    ]
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"[falsifier] wrote {out_md}")

    # exit codes mirror verdict
    return {"PASS": 0, "FAIL": 1, "REVERT": 1, "RUNNER_REGRESSION": 2}[verdict]


if __name__ == "__main__":
    sys.exit(main())
