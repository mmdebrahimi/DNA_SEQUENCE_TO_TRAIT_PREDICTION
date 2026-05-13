"""Cohort audit report generator (/brainstorm Option E + ChatGPT consult).

Emits a markdown audit of a cohort parquet (+ optional AST source) covering:
- Cohort overview (size, per-drug R/S balance, class ratios)
- Clade composition (MLST counts, largest-clade fraction, per-clade R fraction)
- Duplicate detection (same MLST + N50 + contig_count signature)
- Metadata completeness (% missing per field)
- Assembly QC quantiles (N50, contig_count) + below-threshold counts
- AST method breakdown (broth microdilution / disk diffusion / etest / other)
- Breakpoint standard breakdown (CLSI vs EUCAST vs unknown)
- Pre-Gate-B verdict (GO / WARN / NO-GO with rule list)

Phase 2.5 hardening — surfaces real-cohort data-quality issues that mock
fixtures never exercise. The /brainstorm specifically called this out as
core scientific infrastructure: prevents misreading Gate B output as model
failure when it's actually cohort-quality failure.

Phase 2.5+ extensions (deferred):
- ANI / Mash clade clustering (needs Mash CLI install)
- Plasmid-content cross-strain similarity (needs AMRFinder integration)
- More-aggressive duplicate detection (full-genome Mash distance)

Usage:
    uv run python scripts/audit_cohort.py \
        --cohort data/processed/gate_b_cohort.parquet \
        --ast Downloads/BVBRC_genome_amr.csv \
        --output reports/cohort_audit.md
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from dna_decode.data.cohort import CandidateStrain, StrainCohort, load_cohort


PHASE1_DRUGS_DEFAULT = ("ciprofloxacin", "ceftriaxone", "tetracycline")
N50_MIN_DEFAULT = 50_000
CONTIG_MAX_DEFAULT = 500


@dataclass
class VerdictRules:
    """Thresholds for the GO / WARN / NO-GO verdict."""

    target_per_drug: int = 150
    min_minority_class: int = 30
    max_pct_missing_metadata: float = 20.0
    min_pct_broth_microdilution: float = 80.0
    n50_min: int = N50_MIN_DEFAULT
    contig_count_max: int = CONTIG_MAX_DEFAULT


PHASE1_DEFAULT_RULES = VerdictRules()
"""Canonical Phase 1 ship thresholds (150 strains/drug, 30 minority, 20% max missing
metadata, 80% min broth-microdilution, N50 ≥ 50K, contig_count ≤ 500). Any caller
deviating from these triggers the non-default warning banner — prevents the
calibration miscarriage of commit b646fc9 where relaxed CLI flags produced a
context-free 'GO' on the 67-strain cohort."""


def non_default_threshold_fields(
    rules: VerdictRules, default: VerdictRules = PHASE1_DEFAULT_RULES
) -> list[tuple[str, object, object, str]]:
    """Return fields where `rules` is more permissive than `default`.

    Returns a list of (field_name, current_value, default_value, direction) tuples.
    `direction` is "lower" or "higher" describing how `rules` deviates. "Permissive"
    means: lower target_per_drug / lower min_minority_class / higher max_pct_missing_metadata
    / lower min_pct_broth_microdilution / lower n50_min / higher contig_count_max.
    """
    deviations: list[tuple[str, object, object, str]] = []
    # Lower-is-permissive fields
    for fname in ("target_per_drug", "min_minority_class", "min_pct_broth_microdilution", "n50_min"):
        cur = getattr(rules, fname)
        dflt = getattr(default, fname)
        if cur < dflt:
            deviations.append((fname, cur, dflt, "lower"))
    # Higher-is-permissive fields
    for fname in ("max_pct_missing_metadata", "contig_count_max"):
        cur = getattr(rules, fname)
        dflt = getattr(default, fname)
        if cur > dflt:
            deviations.append((fname, cur, dflt, "higher"))
    return deviations


@dataclass
class VerdictResult:
    """Verdict + rule-by-rule outcome."""

    verdict: str  # GO / WARN / NO-GO
    rules: list[tuple[str, str, str]] = field(default_factory=list)
    # tuples of (rule_name, outcome PASS/WARN/FAIL, detail)


def _quantiles(values: list[int]) -> dict[str, int]:
    """5/25/50/75/95 quantiles. Empty list returns zeros."""
    if not values:
        return {"p5": 0, "p25": 0, "p50": 0, "p75": 0, "p95": 0}
    sorted_v = sorted(values)
    n = len(sorted_v)

    def at(p: float) -> int:
        idx = min(n - 1, max(0, int(p * n)))
        return sorted_v[idx]

    return {"p5": at(0.05), "p25": at(0.25), "p50": at(0.5), "p75": at(0.75), "p95": at(0.95)}


def _format_int(n: int) -> str:
    return f"{n:,}"


def _format_pct(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "n/a"
    return f"{100.0 * numerator / denominator:.1f}%"


def overview_section(cohort: StrainCohort, drugs: Iterable[str]) -> list[str]:
    """Cohort overview: size + per-drug R/S balance + class ratios."""
    lines = ["## Cohort overview", ""]
    lines.append(f"- Total isolates: **{_format_int(len(cohort))}**")
    drug_list = list(drugs)
    if not drug_list:
        lines.append("- No drugs evaluated (drug list empty)")
        return lines
    lines.append(f"- Drugs evaluated: {', '.join(drug_list)}")
    lines.append("")
    lines.append("| Drug | n | R | S | Minority | R/S ratio |")
    lines.append("|---|---|---|---|---|---|")
    for drug in drug_list:
        d_lower = drug.lower()
        ids = cohort.per_drug_strain_ids.get(d_lower, [])
        strain_objs = [s for s in cohort.strains if s.strain_id in set(ids)]
        r = sum(1 for s in strain_objs if s.ast_labels.get(d_lower) == 1)
        s_count = sum(1 for s in strain_objs if s.ast_labels.get(d_lower) == 0)
        n = r + s_count
        minority = min(r, s_count)
        ratio = (
            f"{r}:{s_count}" if (r > 0 or s_count > 0) else "0:0"
        )
        lines.append(f"| {drug} | {n} | {r} | {s_count} | {minority} | {ratio} |")
    lines.append("")
    return lines


def clade_section(cohort: StrainCohort, drugs: Iterable[str]) -> list[str]:
    """MLST clade composition + per-clade R fraction (lineage-confound signal)."""
    lines = ["## Clade composition (MLST)", ""]
    mlst_counts = Counter(s.mlst or "(blank)" for s in cohort.strains)
    total = sum(mlst_counts.values())
    largest_count = max(mlst_counts.values()) if mlst_counts else 0
    largest_frac = (100.0 * largest_count / total) if total else 0.0

    lines.append(f"- Unique MLSTs: **{len(mlst_counts)}**")
    lines.append(f"- Largest-clade fraction: **{largest_frac:.1f}%** "
                 f"({largest_count} / {total}) — high values risk lineage dominance")
    lines.append("")
    lines.append("### Top 10 clades by size")
    lines.append("")
    lines.append("| MLST | n | fraction |")
    lines.append("|---|---|---|")
    for mlst, count in mlst_counts.most_common(10):
        frac = 100.0 * count / total if total else 0
        lines.append(f"| {mlst} | {count} | {frac:.1f}% |")
    lines.append("")

    # Per-clade resistant fraction per drug (lineage confounding indicator)
    drug_list = list(drugs)
    if drug_list and mlst_counts:
        lines.append("### Resistant fraction by clade (top 8 clades × Phase 1 drugs)")
        lines.append("")
        lines.append("Pure-R or pure-S clades risk teaching the model lineage signature "
                     "rather than mechanism. Watch for cells = 0% or 100% on multi-strain "
                     "clades.")
        lines.append("")
        header = "| MLST | n | " + " | ".join(f"{d} %R" for d in drug_list) + " |"
        sep = "|---|---|" + "|".join("---" for _ in drug_list) + "|"
        lines.append(header)
        lines.append(sep)
        for mlst, count in mlst_counts.most_common(8):
            row_strains = [s for s in cohort.strains if (s.mlst or "(blank)") == mlst]
            cells = [f"{mlst}", str(count)]
            for drug in drug_list:
                d_lower = drug.lower()
                labeled = [s for s in row_strains if d_lower in s.ast_labels]
                if not labeled:
                    cells.append("—")
                else:
                    r = sum(1 for s in labeled if s.ast_labels[d_lower] == 1)
                    cells.append(f"{100.0*r/len(labeled):.0f}% (n={len(labeled)})")
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")

    # Duplicate-signature detection: same MLST + N50 + contig_count
    sig_counts = Counter(
        (s.mlst, s.n50, s.contig_count) for s in cohort.strains if s.mlst
    )
    dupes = [(sig, n) for sig, n in sig_counts.items() if n > 1]
    if dupes:
        lines.append("### Suspected duplicate isolates")
        lines.append("")
        lines.append("Strains sharing (MLST, N50, contig_count) signature — potential "
                     "near-clones that would leak across CV folds.")
        lines.append("")
        lines.append("| MLST | N50 | contig_count | count |")
        lines.append("|---|---|---|---|")
        for (mlst, n50, contig_count), count in dupes[:10]:
            lines.append(
                f"| {mlst} | {_format_int(n50)} | {contig_count} | {count} |"
            )
        lines.append("")
    return lines


def metadata_section(cohort: StrainCohort) -> list[str]:
    """Per-field metadata-completeness percentages."""
    lines = ["## Metadata completeness", ""]
    n = len(cohort)
    if n == 0:
        lines.append("- Empty cohort.")
        return lines

    missing = {
        "assembly_accession": sum(1 for s in cohort.strains if not s.assembly_accession),
        "mlst": sum(1 for s in cohort.strains if not s.mlst),
        "contig_count": sum(1 for s in cohort.strains if not s.contig_count),
        "n50": sum(1 for s in cohort.strains if not s.n50),
        "country": sum(1 for s in cohort.strains if not s.country),
        "year": sum(1 for s in cohort.strains if not s.year),
    }
    lines.append("| Field | Missing | % Missing |")
    lines.append("|---|---|---|")
    for field_name, m in missing.items():
        lines.append(f"| {field_name} | {m} | {_format_pct(m, n)} |")
    lines.append("")
    return lines


def assembly_qc_section(
    cohort: StrainCohort, n50_min: int, contig_max: int
) -> list[str]:
    """Assembly N50 + contig_count quantiles + below-threshold counts."""
    lines = ["## Assembly QC distribution", ""]
    n50s = [s.n50 for s in cohort.strains if s.n50 > 0]
    contigs = [s.contig_count for s in cohort.strains if s.contig_count > 0]

    if n50s:
        q = _quantiles(n50s)
        lines.append(
            f"**N50** (n={len(n50s)}): "
            f"p5={_format_int(q['p5'])}, p25={_format_int(q['p25'])}, "
            f"p50={_format_int(q['p50'])}, p75={_format_int(q['p75'])}, "
            f"p95={_format_int(q['p95'])}"
        )
        below = sum(1 for v in n50s if v < n50_min)
        lines.append(
            f"- Below threshold (N50 < {_format_int(n50_min)}): "
            f"**{below}** ({_format_pct(below, len(n50s))})"
        )

    if contigs:
        q = _quantiles(contigs)
        lines.append("")
        lines.append(
            f"**Contig count** (n={len(contigs)}): "
            f"p5={q['p5']}, p25={q['p25']}, p50={q['p50']}, "
            f"p75={q['p75']}, p95={q['p95']}"
        )
        above = sum(1 for v in contigs if v > contig_max)
        lines.append(
            f"- Above threshold (contig_count > {contig_max}): "
            f"**{above}** ({_format_pct(above, len(contigs))})"
        )

    lines.append("")
    return lines


def ast_breakdown_section(ast_df, drugs: Iterable[str]) -> list[str]:
    """AST method + breakpoint standard breakdown (requires --ast)."""
    import pandas as pd  # type: ignore[import-not-found]

    lines = ["## AST source breakdown (from --ast)", ""]
    if ast_df is None or len(ast_df) == 0:
        lines.append("- No AST table provided (skip section)")
        lines.append("")
        return lines

    total = len(ast_df)
    lines.append(f"Total AST rows post-organism filter: **{_format_int(total)}**")
    lines.append("")
    method_counts = ast_df["measurement_method"].value_counts()
    lines.append("### Method breakdown")
    lines.append("")
    lines.append("| Method | rows | fraction |")
    lines.append("|---|---|---|")
    for method, count in method_counts.items():
        lines.append(f"| {method or '(blank)'} | {count} | {_format_pct(count, total)} |")
    lines.append("")

    if "source" in ast_df.columns:
        std_counts = ast_df["source"].value_counts()
        lines.append("### Breakpoint standard")
        lines.append("")
        lines.append("| Standard | rows | fraction |")
        lines.append("|---|---|---|")
        for std, count in std_counts.head(8).items():
            label = std if std else "(unknown)"
            lines.append(f"| {label} | {count} | {_format_pct(count, total)} |")
        lines.append("")

    # Per-drug method-mix
    drug_list = [d.lower() for d in drugs]
    if drug_list:
        lines.append("### Per-drug method mix")
        lines.append("")
        methods = sorted(ast_df["measurement_method"].dropna().unique())
        header = "| Drug | n | " + " | ".join(methods) + " |"
        sep = "|---|---|" + "|".join("---" for _ in methods) + "|"
        lines.append(header)
        lines.append(sep)
        for drug in drug_list:
            sub = ast_df[ast_df["antibiotic"] == drug]
            if len(sub) == 0:
                lines.append(f"| {drug} | 0 | " + " | ".join("—" for _ in methods) + " |")
                continue
            row = [drug, str(len(sub))]
            for method in methods:
                m_count = (sub["measurement_method"] == method).sum()
                row.append(f"{m_count} ({_format_pct(m_count, len(sub))})")
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
    return lines


def evaluate_verdict(
    cohort: StrainCohort,
    drugs: Iterable[str],
    ast_df,
    rules: VerdictRules,
) -> VerdictResult:
    """Apply rule list, return GO / WARN / NO-GO."""
    drug_list = list(drugs)
    failures = 0
    warnings = 0
    results: list[tuple[str, str, str]] = []

    # Rule 1: ≥ target_per_drug strains per drug
    for drug in drug_list:
        d_lower = drug.lower()
        ids = cohort.per_drug_strain_ids.get(d_lower, [])
        n = len(ids)
        if n >= rules.target_per_drug:
            results.append((f"strain count ({drug})", "PASS", f"{n} ≥ {rules.target_per_drug}"))
        elif n >= rules.target_per_drug // 3:
            warnings += 1
            results.append((f"strain count ({drug})", "WARN",
                            f"{n} < target {rules.target_per_drug} but above floor"))
        else:
            failures += 1
            results.append((f"strain count ({drug})", "FAIL",
                            f"{n} too few — Phase 1 needs ≥{rules.target_per_drug}"))

    # Rule 2: minority class ≥ min_minority_class per drug
    for drug in drug_list:
        d_lower = drug.lower()
        ids = set(cohort.per_drug_strain_ids.get(d_lower, []))
        labeled = [s for s in cohort.strains if s.strain_id in ids]
        r = sum(1 for s in labeled if s.ast_labels.get(d_lower) == 1)
        s_count = sum(1 for s in labeled if s.ast_labels.get(d_lower) == 0)
        minority = min(r, s_count)
        if minority >= rules.min_minority_class:
            results.append((f"minority class ({drug})", "PASS",
                            f"min(R={r}, S={s_count}) = {minority}"))
        elif minority >= 10:
            warnings += 1
            results.append((f"minority class ({drug})", "WARN",
                            f"min(R={r}, S={s_count}) = {minority} — sparse"))
        else:
            failures += 1
            results.append((f"minority class ({drug})", "FAIL",
                            f"min(R={r}, S={s_count}) = {minority} — below classifier guard"))

    # Rule 3: ≤ max_pct_missing_metadata on each critical field
    n = len(cohort)
    if n > 0:
        for field_name in ("assembly_accession", "mlst", "n50"):
            missing = sum(
                1 for s in cohort.strains if not getattr(s, field_name)
            )
            pct = 100.0 * missing / n
            if pct <= rules.max_pct_missing_metadata:
                results.append((f"{field_name} completeness", "PASS",
                                f"{pct:.1f}% missing"))
            elif pct <= 50.0:
                warnings += 1
                results.append((f"{field_name} completeness", "WARN",
                                f"{pct:.1f}% missing"))
            else:
                failures += 1
                results.append((f"{field_name} completeness", "FAIL",
                                f"{pct:.1f}% missing — substantial signal loss"))

    # Rule 4: ≥ min_pct_broth_microdilution if AST provided
    if ast_df is not None and len(ast_df) > 0:
        bm = (ast_df["measurement_method"] == "broth_microdilution").sum()
        bm_pct = 100.0 * bm / len(ast_df)
        if bm_pct >= rules.min_pct_broth_microdilution:
            results.append(("broth microdilution coverage", "PASS",
                            f"{bm_pct:.1f}%"))
        else:
            warnings += 1
            results.append(("broth microdilution coverage", "WARN",
                            f"{bm_pct:.1f}% — mixed method-mix detected"))

    if failures > 0:
        verdict = "NO-GO"
    elif warnings > 0:
        verdict = "WARN"
    else:
        verdict = "GO"
    return VerdictResult(verdict=verdict, rules=results)


def thresholds_block(rules: VerdictRules) -> list[str]:
    """Render the 'Thresholds applied' block for the report header.

    Lists all six VerdictRules knobs with their current values. Pattern mirrors
    the other `*_section()` helpers — returns a list of markdown lines (no trailing
    join). Empty trailing line gives section separation.
    """
    return [
        "**Thresholds applied:**",
        f"- target_per_drug: {rules.target_per_drug}",
        f"- min_minority_class: {rules.min_minority_class}",
        f"- max_pct_missing_metadata: {rules.max_pct_missing_metadata}",
        f"- min_pct_broth_microdilution: {rules.min_pct_broth_microdilution}",
        f"- n50_min: {rules.n50_min}",
        f"- contig_count_max: {rules.contig_count_max}",
        "",
    ]


def relaxed_flags_warning(rules: VerdictRules) -> list[str]:
    """Top-of-report banner when ANY VerdictRules field is more permissive than Phase 1 default.

    Empty list when rules match defaults (no banner needed). When deviations exist,
    returns a multi-line warning that names each deviated field + its current value
    + the canonical default. Goes before the timestamp/header lines so it's the
    first thing a reader sees.

    Per /brainstorm patches D6 + D7: covers all six VerdictRules fields, not just
    target_per_drug and min_minority_class.
    """
    deviations = non_default_threshold_fields(rules)
    if not deviations:
        return []
    lines = [
        "> :warning: **WARNING: NON-DEFAULT THRESHOLDS APPLIED — VERDICT NOT COMPARABLE TO PHASE 1 GATE.**",
        ">",
        "> The following thresholds deviate from canonical Phase 1 defaults:",
        ">",
    ]
    for fname, cur, dflt, direction in deviations:
        lines.append(f"> - `{fname}` = {cur} ({direction} than default {dflt})")
    lines.append("")
    return lines


def stdout_warning_lines(rules: VerdictRules) -> list[str]:
    """Stdout-side equivalent of relaxed_flags_warning (covers the CLI escape hatch).

    Per /brainstorm patch D6: report-side banner alone left a stdout escape hatch
    where CI/pipeline consumers received a context-free `verdict: GO`. This helper
    returns the per-line strings to emit BEFORE the verdict echo on any deviation.
    """
    deviations = non_default_threshold_fields(rules)
    if not deviations:
        return []
    summary = ", ".join(f"{name}={cur}" for name, cur, _, _ in deviations)
    return [
        f"[audit_cohort] WARNING: relaxed thresholds ({summary}) — "
        f"verdict not comparable to Phase 1 gate"
    ]


def verdict_section(result: VerdictResult) -> list[str]:
    lines = [f"## Pre-Gate-B verdict: **{result.verdict}**", ""]
    if not result.rules:
        lines.append("No rules evaluated.")
        return lines
    lines.append("| Rule | Outcome | Detail |")
    lines.append("|---|---|---|")
    for name, outcome, detail in result.rules:
        lines.append(f"| {name} | {outcome} | {detail} |")
    lines.append("")
    return lines


def build_report(
    cohort: StrainCohort,
    ast_df,
    drugs: Iterable[str],
    rules: VerdictRules,
    cohort_path: Path,
) -> str:
    """Compose the full markdown report."""
    import datetime as _dt

    drug_list = list(drugs)
    lines = [
        "# Cohort Audit Report",
        "",
    ]
    # Banner BEFORE source/timestamp so it's the first thing a reader sees on a
    # relaxed-thresholds run (per /brainstorm patch D6 + D7).
    lines += relaxed_flags_warning(rules)
    lines += [
        f"**Source cohort:** `{cohort_path}`",
        f"**Generated:** {_dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
    ]
    lines += thresholds_block(rules)
    lines += overview_section(cohort, drug_list)
    lines += clade_section(cohort, drug_list)
    lines += metadata_section(cohort)
    lines += assembly_qc_section(cohort, rules.n50_min, rules.contig_count_max)
    if ast_df is not None:
        lines += ast_breakdown_section(ast_df, drug_list)
    verdict = evaluate_verdict(cohort, drug_list, ast_df, rules)
    lines += verdict_section(verdict)
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cohort audit report generator (Phase 2.5 hardening)."
    )
    parser.add_argument("--cohort", required=True, type=Path)
    parser.add_argument(
        "--ast",
        default=None,
        type=Path,
        help="Optional BV-BRC AMR CSV/TSV for method + breakpoint breakdown.",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Markdown output path"
    )
    parser.add_argument(
        "--drugs",
        default=",".join(PHASE1_DRUGS_DEFAULT),
        help="Comma-separated drug list for per-drug breakdowns",
    )
    parser.add_argument("--target-per-drug", type=int, default=150)
    parser.add_argument("--min-minority-class", type=int, default=30)
    parser.add_argument("--max-pct-missing-metadata", type=float, default=20.0)
    parser.add_argument("--min-pct-broth-microdilution", type=float, default=80.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.cohort.exists():
        print(f"[audit_cohort] cohort not found: {args.cohort}", file=sys.stderr)
        return 2

    cohort = load_cohort(args.cohort)
    drugs = tuple(d.strip() for d in args.drugs.split(",") if d.strip())

    ast_df = None
    if args.ast:
        if not args.ast.exists():
            print(f"[audit_cohort] --ast file not found: {args.ast}", file=sys.stderr)
            return 2
        from dna_decode.data.ast_data import load_bvbrc_ast

        ast_df = load_bvbrc_ast(args.ast)

    rules = VerdictRules(
        target_per_drug=args.target_per_drug,
        min_minority_class=args.min_minority_class,
        max_pct_missing_metadata=args.max_pct_missing_metadata,
        min_pct_broth_microdilution=args.min_pct_broth_microdilution,
    )

    report = build_report(cohort, ast_df, drugs, rules, args.cohort)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    print(f"[audit_cohort] wrote report: {args.output} ({len(report)} chars)")

    # Echo the verdict line to stdout for CI / pipeline integration.
    # Per /brainstorm patch D6: prepend a WARNING line if any threshold deviates
    # from Phase 1 default — stdout-side coverage of the calibration discipline.
    verdict = evaluate_verdict(cohort, drugs, ast_df, rules)
    for warning_line in stdout_warning_lines(rules):
        print(warning_line)
    print(f"[audit_cohort] verdict: {verdict.verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
