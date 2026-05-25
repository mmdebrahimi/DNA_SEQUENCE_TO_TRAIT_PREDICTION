"""BV-BRC strict-MIC 4-drug feasibility census.

Phase 2 entry per `project_state/dna-decode-2026-05-11.md` Candidate-next-actions
row 1. Counts per-drug feasibility of building N=150 strict-MIC cohorts from
BV-BRC AST data, gating each candidate strain through:

  Stage 1: E. coli broth-microdilution AST rows for the drug
  Stage 2: distinct genome IDs
  Stage 3: with parseable MIC value
  Stage 4: HIGH_R or HIGH_S strict-MIC classification (per drug breakpoints)
  Stage 5: with assembly_accession (NCBI-downloadable)
  Stage 6: passing assembly QC (contig_count <= 500, n50 >= 50_000)

Per-drug breakpoints (CLSI 2024 + EUCAST 14.0 for E. coli):

  - Ciprofloxacin: CLSI R>=2, S<=0.5  / EUCAST R>=1, S<=0.25
  - Ceftriaxone:   CLSI R>=4, S<=1    / EUCAST R>=2, S<=1
  - Tetracycline:  CLSI R>=16, S<=4   / EUCAST: no breakpoints for E. coli (ECOFF only)
  - Gentamicin:    CLSI R>=16, S<=4   / EUCAST R>=4, S<=2  (aminoglycoside; 4th-mechanism-class)

HIGH_R = median MIC >= 4 * CLSI-R breakpoint (4x safety margin)
HIGH_S = median MIC <= CLSI-S / 4

Strict-MIC pass = HIGH_R or HIGH_S only (BORDERLINE / AMBIGUOUS / CONFLICT / NO_MIC dropped).

Output: wiki/bvbrc_strict_mic_4drug_census_<date>.md + .json sidecar.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import date as _date
from math import isnan
from pathlib import Path
from typing import Optional

import pandas as pd

from dna_decode.data.ast_data import load_bvbrc_ast
from dna_decode.data.bvbrc_genome import load_bvbrc_genome_metadata
from dna_decode.data.mic_tiers import (
    DRUG_BREAKPOINTS,
    breakpoints_for,
    classify_tier as classify_strict_mic,
)

# Assembly QC thresholds (Stage 2 cohort builder default)
DEFAULT_QC_MAX_CONTIGS = 500
DEFAULT_QC_MIN_N50 = 50_000


# ---------------------------------------------------------------------------
# Per-drug census
# ---------------------------------------------------------------------------


def census_drug(
    ast: pd.DataFrame,
    metadata: dict,
    drug: str,
    qc_max_contigs: int = DEFAULT_QC_MAX_CONTIGS,
    qc_min_n50: int = DEFAULT_QC_MIN_N50,
) -> dict:
    """Run the 6-stage feasibility census for one drug.

    Returns a dict with stage counts + per-strain HIGH_R/HIGH_S breakdown
    + the list of feasible strain_ids at the final stage.
    """
    breakpoints = DRUG_BREAKPOINTS.get(drug.lower())
    if breakpoints is None:
        raise ValueError(f"No breakpoints configured for drug: {drug!r}")

    drug_ast = ast[ast["antibiotic"].str.lower() == drug.lower()].copy()
    stage1_total_ast_rows = len(drug_ast)

    distinct_genomes = drug_ast["strain_id"].dropna().unique()
    stage2_distinct_genomes = len(distinct_genomes)

    # Aggregate MIC + label per strain
    per_strain: dict[str, dict] = defaultdict(lambda: {"mics": [], "labels": set()})
    for _, row in drug_ast.iterrows():
        sid = str(row["strain_id"]).strip()
        if not sid:
            continue
        mic = row.get("mic_value")
        if mic is not None and not (isinstance(mic, float) and isnan(mic)):
            per_strain[sid]["mics"].append(float(mic))
        label = str(row.get("susceptibility_label", "")).strip().upper()
        if label:
            per_strain[sid]["labels"].add(label)

    stage3_with_mic = sum(1 for s in per_strain.values() if s["mics"])

    # Classify each strain. `load_bvbrc_ast` keeps `susceptibility_label` as the
    # uppercased string ("R" / "RESISTANT" / "S" / "SUSCEPTIBLE" / "I" /
    # "INTERMEDIATE"); binarization happens downstream at `get_binary_labels` /
    # `candidates_from_bvbrc_ast`, not here. We pass the raw label set directly.
    classifications: dict[str, str] = {}
    strict_r: list[str] = []      # HIGH_R only
    strict_s: list[str] = []      # HIGH_S only
    decisive_r: list[str] = []    # DECISIVE_R (between 2x and 4x CLSI-R)
    decisive_s: list[str] = []    # DECISIVE_S (between CLSI_S/4 and CLSI_S/2)
    bucket_counts: dict[str, int] = defaultdict(int)
    for sid, agg in per_strain.items():
        distinct_calls = set(agg["labels"])
        tier = classify_strict_mic(agg["mics"], distinct_calls, breakpoints)
        classifications[sid] = tier
        bucket_counts[tier] += 1
        if tier == "HIGH_R":
            strict_r.append(sid)
        elif tier == "HIGH_S":
            strict_s.append(sid)
        elif tier == "DECISIVE_R":
            decisive_r.append(sid)
        elif tier == "DECISIVE_S":
            decisive_s.append(sid)

    stage4_strict_mic = len(strict_r) + len(strict_s)
    stage4_relaxed_mic = stage4_strict_mic + len(decisive_r) + len(decisive_s)

    # Stage 5: with assembly_accession. Note: `load_bvbrc_genome_metadata`
    # already drops rows lacking assembly_accession at parse time, so this
    # check primarily filters strains whose strain_id has no metadata entry
    # at all (returns {} → empty accession).
    def _has_accession(sid: str) -> bool:
        meta = metadata.get(sid, {})
        acc = str(meta.get("assembly_accession", "")).strip()
        return bool(acc) and acc.lower() not in ("nan", "none")

    # Stage 6: assembly QC
    def _passes_qc(sid: str) -> bool:
        meta = metadata.get(sid, {})
        try:
            contigs = int(meta.get("contig_count", 0) or 0)
            n50 = int(meta.get("n50", 0) or 0)
        except (TypeError, ValueError):
            return False
        return contigs > 0 and contigs <= qc_max_contigs and n50 >= qc_min_n50

    def _filter_pipeline(strains: list[str]) -> tuple[list[str], list[str]]:
        with_acc = [s for s in strains if _has_accession(s)]
        with_qc = [s for s in with_acc if _passes_qc(s)]
        return with_acc, with_qc

    strict_r_acc, strict_r_qc = _filter_pipeline(strict_r)
    strict_s_acc, strict_s_qc = _filter_pipeline(strict_s)
    decisive_r_acc, decisive_r_qc = _filter_pipeline(decisive_r)
    decisive_s_acc, decisive_s_qc = _filter_pipeline(decisive_s)

    relaxed_r_qc = strict_r_qc + decisive_r_qc
    relaxed_s_qc = strict_s_qc + decisive_s_qc

    stage5_strict_with_acc = len(strict_r_acc) + len(strict_s_acc)
    stage5_relaxed_with_acc = stage5_strict_with_acc + len(decisive_r_acc) + len(decisive_s_acc)
    stage6_strict_qc = len(strict_r_qc) + len(strict_s_qc)
    stage6_relaxed_qc = stage6_strict_qc + len(decisive_r_qc) + len(decisive_s_qc)

    return {
        "drug": drug,
        "breakpoints": {
            "clsi_r": breakpoints["clsi_r"], "clsi_s": breakpoints["clsi_s"],
            "eucast_r": breakpoints.get("eucast_r"), "eucast_s": breakpoints.get("eucast_s"),
        },
        "stages": {
            "1_total_ast_rows": stage1_total_ast_rows,
            "2_distinct_genomes": stage2_distinct_genomes,
            "3_with_mic_value": stage3_with_mic,
            "4_strict_mic_pass": stage4_strict_mic,
            "4_relaxed_mic_pass": stage4_relaxed_mic,
            "5_strict_with_assembly_accession": stage5_strict_with_acc,
            "5_relaxed_with_assembly_accession": stage5_relaxed_with_acc,
            "6_strict_passing_assembly_qc": stage6_strict_qc,
            "6_relaxed_passing_assembly_qc": stage6_relaxed_qc,
        },
        "tier_buckets": dict(bucket_counts),
        "strict": {
            "high_r": len(strict_r_qc),
            "high_s": len(strict_s_qc),
            "total": stage6_strict_qc,
        },
        "relaxed": {
            "high_or_decisive_r": len(relaxed_r_qc),
            "high_or_decisive_s": len(relaxed_s_qc),
            "total": stage6_relaxed_qc,
        },
        "feasibility_strict": {
            "n150_per_class": (len(strict_r_qc) >= 75 and len(strict_s_qc) >= 75),
            "n100_per_class": (len(strict_r_qc) >= 50 and len(strict_s_qc) >= 50),
            "n60_smoke": (len(strict_r_qc) >= 30 and len(strict_s_qc) >= 30),
        },
        "feasibility_relaxed": {
            "n150_per_class": (len(relaxed_r_qc) >= 75 and len(relaxed_s_qc) >= 75),
            "n100_per_class": (len(relaxed_r_qc) >= 50 and len(relaxed_s_qc) >= 50),
            "n60_smoke": (len(relaxed_r_qc) >= 30 and len(relaxed_s_qc) >= 30),
        },
        "feasible_strain_ids": {
            "strict_r": strict_r_qc,
            "strict_s": strict_s_qc,
            "decisive_r": decisive_r_qc,
            "decisive_s": decisive_s_qc,
        },
    }


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------


def _verdict_string(feasibility: dict, r_count: int, s_count: int, label: str) -> str:
    if feasibility["n150_per_class"]:
        return f"**{label}: N=150 PER-CLASS (75R/75S) FEASIBLE** (R={r_count}, S={s_count})"
    if feasibility["n100_per_class"]:
        return f"**{label}: N=100 PER-CLASS feasible; N=150 short** (R={r_count}, S={s_count})"
    if feasibility["n60_smoke"]:
        return f"**{label}: N=60 smoke feasible; N=100+ short** (R={r_count}, S={s_count})"
    return f"**{label}: INFEASIBLE at any N>=60** (R={r_count}, S={s_count})"


def render_markdown(results: list[dict], output_path: Path) -> None:
    today = _date.today().isoformat()
    lines = [
        f"# BV-BRC MIC 4-drug feasibility census — {today}",
        "",
        "Phase 2 entry: counts per-drug feasibility of building N=150 cohorts from BV-BRC AST data at two label-quality bars.",
        "",
        "## Methodology",
        "",
        "Strains gated through 6 stages, computed at TWO label-quality bars:",
        "",
        "- **Strict-MIC** = HIGH_R + HIGH_S only (median MIC >=4x CLSI-R or <=CLSI-S/4; 4x safety margin)",
        "- **Relaxed-MIC** = strict + DECISIVE_R + DECISIVE_S (median MIC outside borderline gray zone but no 4x margin)",
        "",
        "Strains pass through:",
        "1. E. coli broth-microdilution AST rows for the drug",
        "2. Distinct genome IDs",
        "3. With parseable MIC value (numeric, not NA)",
        "4. Strict-MIC OR relaxed-MIC classification",
        "5. With assembly_accession (NCBI-downloadable)",
        f"6. Passing assembly QC (contig_count <= {DEFAULT_QC_MAX_CONTIGS}, n50 >= {DEFAULT_QC_MIN_N50})",
        "",
        "Breakpoints: CLSI 2024 + EUCAST 14.0 (E. coli).",
        "",
        "**Why two bars:** the audit framework (mechanism x MIC x opacity merge with SUSPEND gate) is the downstream cohort-quality gate. Strict-MIC pre-filters on a 4x safety margin upstream — biology-aware noise reduction, but doing the audit framework's job twice. Relaxed-MIC accepts any MIC clearly outside the gray zone and lets the audit framework flag noisy strains downstream. Both numbers reported so the user can pick the bar that matches the project's failure-tolerance.",
        "",
        "## Per-drug census",
        "",
    ]

    for r in results:
        bp = r["breakpoints"]
        eucast = f"EUCAST R>={bp['eucast_r']}, S<={bp['eucast_s']}" if bp["eucast_r"] else "EUCAST: no breakpoints (ECOFF only)"
        lines.extend([
            f"### {r['drug'].capitalize()}",
            "",
            f"Breakpoints: CLSI R>={bp['clsi_r']}, S<={bp['clsi_s']} / {eucast}",
            "",
            "| Stage | Strict-MIC | Relaxed-MIC |",
            "|---|---:|---:|",
            f"| 1. Total E. coli broth-MIC AST rows | {r['stages']['1_total_ast_rows']:,} | {r['stages']['1_total_ast_rows']:,} |",
            f"| 2. Distinct genome IDs | {r['stages']['2_distinct_genomes']:,} | {r['stages']['2_distinct_genomes']:,} |",
            f"| 3. With parseable MIC value | {r['stages']['3_with_mic_value']:,} | {r['stages']['3_with_mic_value']:,} |",
            f"| 4. MIC classification pass | {r['stages']['4_strict_mic_pass']:,} | {r['stages']['4_relaxed_mic_pass']:,} |",
            f"| 5. With assembly_accession | {r['stages']['5_strict_with_assembly_accession']:,} | {r['stages']['5_relaxed_with_assembly_accession']:,} |",
            f"| 6. Passing assembly QC | {r['stages']['6_strict_passing_assembly_qc']:,} | {r['stages']['6_relaxed_passing_assembly_qc']:,} |",
            "",
            "Tier distribution (all distinct genome IDs, stage 2):",
            "",
            "| Tier | Count |",
            "|---|---:|",
        ])
        for tier, count in sorted(r["tier_buckets"].items(), key=lambda kv: -kv[1]):
            lines.append(f"| {tier} | {count:,} |")
        lines.extend([
            "",
            f"**Strict-MIC final:** R={r['strict']['high_r']}, S={r['strict']['high_s']}, Total={r['strict']['total']}",
            f"**Relaxed-MIC final:** R={r['relaxed']['high_or_decisive_r']}, S={r['relaxed']['high_or_decisive_s']}, Total={r['relaxed']['total']}",
            "",
            f"**Verdict:** {_verdict_string(r['feasibility_strict'], r['strict']['high_r'], r['strict']['high_s'], 'Strict-MIC')}",
            f"**Verdict:** {_verdict_string(r['feasibility_relaxed'], r['relaxed']['high_or_decisive_r'], r['relaxed']['high_or_decisive_s'], 'Relaxed-MIC')}",
            "",
        ])

    lines.extend([
        "## Cross-drug summary",
        "",
        "### Strict-MIC (HIGH_R + HIGH_S only; 4x safety margin)",
        "",
        "| Drug | R | S | Total | N=150 per-class? | N=100 per-class? | N=60 smoke? |",
        "|---|---:|---:|---:|---|---|---|",
    ])
    for r in results:
        fs = r["feasibility_strict"]
        check_150 = "YES" if fs["n150_per_class"] else "no"
        check_100 = "YES" if fs["n100_per_class"] else "no"
        check_60 = "YES" if fs["n60_smoke"] else "no"
        lines.append(
            f"| {r['drug'].capitalize()} | {r['strict']['high_r']:,} | {r['strict']['high_s']:,} "
            f"| {r['strict']['total']:,} | {check_150} | {check_100} | {check_60} |"
        )
    lines.extend([
        "",
        "### Relaxed-MIC (HIGH + DECISIVE; lean on audit framework for noise)",
        "",
        "| Drug | R | S | Total | N=150 per-class? | N=100 per-class? | N=60 smoke? |",
        "|---|---:|---:|---:|---|---|---|",
    ])
    for r in results:
        fr = r["feasibility_relaxed"]
        check_150 = "YES" if fr["n150_per_class"] else "no"
        check_100 = "YES" if fr["n100_per_class"] else "no"
        check_60 = "YES" if fr["n60_smoke"] else "no"
        lines.append(
            f"| {r['drug'].capitalize()} | {r['relaxed']['high_or_decisive_r']:,} | {r['relaxed']['high_or_decisive_s']:,} "
            f"| {r['relaxed']['total']:,} | {check_150} | {check_100} | {check_60} |"
        )
    lines.append("")

    lines.extend([
        "## What this informs",
        "",
        "- **Stage 2 cohort building**: which drugs have enough strict-MIC strains for N=150 per-class cohorts.",
        "- **Drug-parameterize the audit infrastructure** ([[Candidate-1-framing]]): which 4th-mechanism-class drug (gentamicin = aminoglycoside) clears the smoke threshold for the architectural-finding falsifier.",
        "- **Compute budget allocation**: only drugs that pass N=150 per-class warrant Databricks NT cache populates.",
        "",
        "## What this does NOT do",
        "",
        "- Does not count MLST coverage on the feasible subset (separate diagnostic at `scripts/diagnose_bvbrc_mlst_gaps.py`)",
        "- Does not check NCBI Datasets API availability for each accession (downloadable in practice may be lower than `assembly_accession`-present)",
        "- Does not run AMRFinder mechanism audit on the feasible subset (separate step; only worthwhile post-cohort-build)",
        "",
    ])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")


def render_json(results: list[dict], output_path: Path) -> None:
    """Emit JSON sidecar. Strain-id lists omitted to keep file small; the
    in-process `census_drug()` return dict carries them under
    `feasible_strain_ids` if a caller needs them.
    """
    payload = {
        "generated": _date.today().isoformat(),
        "qc_max_contigs": DEFAULT_QC_MAX_CONTIGS,
        "qc_min_n50": DEFAULT_QC_MIN_N50,
        "drugs": [
            {
                "drug": r["drug"],
                "breakpoints": r["breakpoints"],
                "stages": r["stages"],
                "tier_buckets": r["tier_buckets"],
                "strict": r["strict"],
                "relaxed": r["relaxed"],
                "feasibility_strict": r["feasibility_strict"],
                "feasibility_relaxed": r["feasibility_relaxed"],
            }
            for r in results
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_cohort_from_census(
    census_result: dict,
    metadata: dict,
    label_quality: str,
    target_per_class: int,
    balance_slack: int,
    output_path: Path,
) -> dict:
    """Build a cohort parquet from a census_drug result.

    Args:
      census_result: output of census_drug() for ONE drug.
      metadata: same metadata dict passed to census_drug (for MLST + accession lookup).
      label_quality: "strict" (HIGH_R + HIGH_S only) or "relaxed" (HIGH + DECISIVE).
      target_per_class: target R + S count separately.
      balance_slack: acceptable per-class deviation.
      output_path: cohort parquet output.

    Returns dict with selected R/S counts + saved path. Raises ValueError if
    feasibility doesn't permit the target.
    """
    # Local imports to avoid pulling cohort/StrainCohort into the census module
    # by default (keeps the feasibility-only mode lightweight).
    from dna_decode.data.cohort import (
        CandidateStrain,
        StrainCohort,
        _mlst_balanced_selection,
        save_cohort,
    )

    drug = census_result["drug"]
    feasible = census_result["feasible_strain_ids"]

    if label_quality == "strict":
        r_ids = list(feasible["strict_r"])
        s_ids = list(feasible["strict_s"])
    elif label_quality == "relaxed":
        r_ids = list(feasible["strict_r"]) + list(feasible["decisive_r"])
        s_ids = list(feasible["strict_s"]) + list(feasible["decisive_s"])
    else:
        raise ValueError(f"unknown label_quality: {label_quality!r}")

    # Hydrate strain_ids -> CandidateStrain objects via metadata (for MLST + accession)
    def _make_candidate(sid: str, label: int) -> CandidateStrain | None:
        meta = metadata.get(sid, {})
        acc = str(meta.get("assembly_accession", "")).strip()
        mlst = str(meta.get("mlst", "") or "").strip()
        if not acc or not mlst or mlst.lower() in ("none", "nan"):
            return None
        try:
            contigs = int(meta.get("contig_count", 0) or 0)
            n50 = int(meta.get("n50", 0) or 0)
        except (TypeError, ValueError):
            return None
        return CandidateStrain(
            strain_id=sid,
            assembly_accession=acc,
            mlst=mlst,
            country=str(meta.get("country", "") or ""),
            year=str(meta.get("year", "") or ""),
            contig_count=contigs,
            n50=n50,
            ast_labels={drug: label},
        )

    r_candidates = [c for c in (_make_candidate(s, 1) for s in r_ids) if c is not None]
    s_candidates = [c for c in (_make_candidate(s, 0) for s in s_ids) if c is not None]

    if len(r_candidates) < target_per_class - balance_slack:
        raise ValueError(
            f"{drug} R pool size {len(r_candidates)} is below target_per_class - balance_slack "
            f"({target_per_class - balance_slack}); cannot build cohort at this tier"
        )
    if len(s_candidates) < target_per_class - balance_slack:
        raise ValueError(
            f"{drug} S pool size {len(s_candidates)} is below target_per_class - balance_slack "
            f"({target_per_class - balance_slack}); cannot build cohort at this tier"
        )

    # Per-class MLST-balanced selection (matches build_stage2_n150_cohort.py pattern)
    selected_r = _mlst_balanced_selection(r_candidates, target_per_class)
    selected_s = _mlst_balanced_selection(s_candidates, target_per_class)

    selected_strains = selected_r + selected_s
    strain_ids_list = [c.strain_id for c in selected_strains]
    cohort = StrainCohort(
        strains=selected_strains,
        per_drug_strain_ids={drug: strain_ids_list},
        three_drug_intersection=[],
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    save_cohort(cohort, output_path)

    return {
        "drug": drug,
        "label_quality": label_quality,
        "r_count": len(selected_r),
        "s_count": len(selected_s),
        "output_path": str(output_path),
        "pool_sizes": {"r_pool": len(r_candidates), "s_pool": len(s_candidates)},
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--ast-csv",
        type=Path,
        default=Path("C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv"),
        help="BV-BRC AST CSV (default: user's Downloads)",
    )
    p.add_argument(
        "--genome-csv",
        type=Path,
        default=Path("C:/Users/Farshad/Downloads/BVBRC_genome (1).csv"),
        help="BV-BRC Genomes CSV with assembly metadata (default: user's Downloads)",
    )
    p.add_argument(
        "--drugs",
        default="ciprofloxacin,ceftriaxone,tetracycline,gentamicin",
        help="Comma-separated drug names (default: cipro+cef+tet+gentamicin)",
    )
    p.add_argument(
        "--output-md",
        type=Path,
        default=None,
        help="Markdown output path (default: wiki/bvbrc_strict_mic_4drug_census_<date>.md)",
    )
    p.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="JSON sidecar path (default: alongside markdown)",
    )
    p.add_argument(
        "--qc-max-contigs", type=int, default=DEFAULT_QC_MAX_CONTIGS,
        help=f"Assembly QC max contig count (default: {DEFAULT_QC_MAX_CONTIGS})",
    )
    p.add_argument(
        "--qc-min-n50", type=int, default=DEFAULT_QC_MIN_N50,
        help=f"Assembly QC min N50 (default: {DEFAULT_QC_MIN_N50})",
    )
    # EP-2.5 (2026-05-25): cohort-build mode extension. Lets ONE invocation
    # do feasibility + cohort-build in one step. Satisfies the EP-2.5
    # terminal claim from plans/Post_V0_EP_Ladder_Plan.md.
    p.add_argument(
        "--build",
        default=None,
        help="If set, build a cohort parquet for this drug AFTER running feasibility. "
             "Must be one of the drugs in --drugs.",
    )
    p.add_argument(
        "--label-quality",
        choices=["strict", "relaxed"],
        default="strict",
        help="When --build is set: which label-quality tier to pick from. "
             "strict = HIGH_R + HIGH_S only (4x safety margin). "
             "relaxed = HIGH + DECISIVE in both directions. Default: strict.",
    )
    p.add_argument(
        "--target-per-class",
        type=int,
        default=75,
        help="When --build is set: target R-count + S-count separately (default 75).",
    )
    p.add_argument(
        "--balance-slack",
        type=int,
        default=15,
        help="When --build is set: acceptable per-class deviation (default +/-15).",
    )
    p.add_argument(
        "--cohort-output",
        type=Path,
        default=None,
        help="When --build is set: cohort parquet output path "
             "(default: data/processed/<drug>_<label_quality>_cohort_<date>.parquet).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.ast_csv.exists():
        print(f"[census] AST CSV not found: {args.ast_csv}", file=sys.stderr)
        return 2
    if not args.genome_csv.exists():
        print(f"[census] Genome CSV not found: {args.genome_csv}", file=sys.stderr)
        return 2

    drugs = [d.strip().lower() for d in args.drugs.split(",") if d.strip()]
    unknown = [d for d in drugs if d not in DRUG_BREAKPOINTS]
    if unknown:
        print(
            f"[census] unknown drugs (no breakpoints configured): {unknown}",
            file=sys.stderr,
        )
        print(f"  configured: {sorted(DRUG_BREAKPOINTS.keys())}", file=sys.stderr)
        return 2

    today = _date.today().isoformat()
    if args.output_md is None:
        args.output_md = Path(f"wiki/bvbrc_strict_mic_4drug_census_{today}.md")
    if args.output_json is None:
        args.output_json = args.output_md.with_suffix(".json")

    print(f"[census] loading AST CSV: {args.ast_csv}")
    ast = load_bvbrc_ast(args.ast_csv)
    print(f"[census]   AST rows after E. coli + broth-MIC filter: {len(ast):,}")

    print(f"[census] loading Genomes CSV: {args.genome_csv}")
    metadata = load_bvbrc_genome_metadata(args.genome_csv)
    print(f"[census]   metadata entries: {len(metadata):,}")

    results: list[dict] = []
    for drug in drugs:
        print(f"[census] running census for {drug}...")
        r = census_drug(
            ast, metadata, drug,
            qc_max_contigs=args.qc_max_contigs,
            qc_min_n50=args.qc_min_n50,
        )
        strict, relaxed = r["strict"], r["relaxed"]
        print(
            f"[census]   {drug}: strict R={strict['high_r']}/S={strict['high_s']} "
            f"({strict['total']}) | relaxed R={relaxed['high_or_decisive_r']}/"
            f"S={relaxed['high_or_decisive_s']} ({relaxed['total']})"
        )
        results.append(r)

    render_markdown(results, args.output_md)
    render_json(results, args.output_json)

    print(f"\n[census] markdown: {args.output_md}")
    print(f"[census] json:     {args.output_json}")

    # Print final verdicts to stdout, both bars
    print("\n[census] === STRICT-MIC FEASIBILITY (HIGH only, 4x safety margin) ===")
    for r in results:
        v = _verdict_string(r["feasibility_strict"], r["strict"]["high_r"], r["strict"]["high_s"], "")
        print(f"  {r['drug']:>16s}: {v.lstrip('*').strip()}")

    print("\n[census] === RELAXED-MIC FEASIBILITY (HIGH + DECISIVE; audit framework as downstream gate) ===")
    for r in results:
        v = _verdict_string(r["feasibility_relaxed"], r["relaxed"]["high_or_decisive_r"], r["relaxed"]["high_or_decisive_s"], "")
        print(f"  {r['drug']:>16s}: {v.lstrip('*').strip()}")

    # EP-2.5 (2026-05-25): cohort-build mode. Runs AFTER feasibility so the
    # report is always emitted, even if cohort-build fails.
    if args.build:
        build_drug = args.build.lower()
        target_result = next((r for r in results if r["drug"] == build_drug), None)
        if target_result is None:
            print(
                f"\n[build] ERROR: --build drug {args.build!r} was not in --drugs list; "
                f"add it to --drugs and re-run.",
                file=sys.stderr,
            )
            return 2
        if args.cohort_output is None:
            args.cohort_output = Path(
                f"data/processed/{build_drug}_{args.label_quality}_cohort_{today}.parquet"
            )
        print(
            f"\n[build] building cohort for {build_drug} @ label_quality={args.label_quality}, "
            f"target_per_class={args.target_per_class}, slack=+/-{args.balance_slack}"
        )
        try:
            build_result = build_cohort_from_census(
                target_result,
                metadata,
                label_quality=args.label_quality,
                target_per_class=args.target_per_class,
                balance_slack=args.balance_slack,
                output_path=args.cohort_output,
            )
        except ValueError as e:
            print(f"[build] ERROR: {e}", file=sys.stderr)
            return 3
        print(
            f"[build] cohort: {build_result['r_count']}R + {build_result['s_count']}S "
            f"(from pools R={build_result['pool_sizes']['r_pool']} / "
            f"S={build_result['pool_sizes']['s_pool']})"
        )
        print(f"[build] saved: {build_result['output_path']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
