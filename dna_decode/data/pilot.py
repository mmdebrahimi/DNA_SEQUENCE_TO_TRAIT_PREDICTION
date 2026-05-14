"""Step 0.5 — Real-data pilot gate (HARD).

Estimates per-drug labeled-isolate counts after filters. Two paths:

1. **Local-TSV path (runnable today; Wave 1.5 hardening)**: reads a pre-downloaded
   BV-BRC AST TSV from a path supplied via `--ast-tsv <path>`, `config.bvbrc_ast.
   local_tsv_path`, or the `BVBRC_AST_TSV` env var. Reuses `load_bvbrc_ast` from
   Step 5 for parsing.

2. **Live-API path (deferred; NotImplementedError)**: real BV-BRC REST endpoint
   selection deferred until first real-data run. The local-TSV path is the
   preferred Phase 1 mode because it operates on a downloadable artifact.

HARD-gate semantics: returns exit-code 0 only when every Phase 1 drug has
>= target_per_drug strains after all filters AND the 3-drug intersection
target is met. Non-zero exit halts /execute-plan at Step 6 (cohort).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from dna_decode.data.ast_data import DEFAULT_ORGANISM, load_bvbrc_ast

DEFAULT_CONFIG_PATH = Path("config/datasources.yaml")


@dataclass(frozen=True)
class CohortSelectionCriteria:
    """Filter parameters for cohort selection. Mirrors Step 6's CohortSelectionCriteria."""

    target_per_drug: int = 150
    measurement_method_filter: tuple[str, ...] = ("broth_microdilution",)
    assembly_contig_count_max: int = 500
    assembly_n50_min: int = 50_000
    three_drug_intersection_target: int = 75


@dataclass
class PerDrugStageCounts:
    """Strain counts at each filter stage for a single drug."""

    drug: str
    raw: int = 0
    after_method_filter: int = 0
    after_assembly_filter: int = 0


@dataclass
class PilotReport:
    """Result of the pilot gate run."""

    drugs: tuple[str, ...]
    criteria: CohortSelectionCriteria
    per_drug: list[PerDrugStageCounts] = field(default_factory=list)
    three_drug_intersection: int = 0
    go_no_go: str = "PENDING"
    failure_reasons: list[str] = field(default_factory=list)
    estimated_download_gb: float = 0.0
    estimated_embedding_minutes_rtx4090: float = 0.0

    def decide(self) -> str:
        """Compute GO/NO-GO verdict based on per-drug + intersection counts."""
        reasons: list[str] = []
        for d in self.per_drug:
            if d.after_assembly_filter < self.criteria.target_per_drug:
                reasons.append(
                    f"{d.drug}: {d.after_assembly_filter} strains after filters "
                    f"< target {self.criteria.target_per_drug}"
                )
        if self.three_drug_intersection < self.criteria.three_drug_intersection_target:
            reasons.append(
                f"3-drug intersection: {self.three_drug_intersection} "
                f"< target {self.criteria.three_drug_intersection_target}"
            )
        self.failure_reasons = reasons
        self.go_no_go = "GO" if not reasons else "NO-GO"
        return self.go_no_go


class PilotGateError(Exception):
    """Raised when the pilot gate cannot complete (network failure, missing config, etc.)."""


def _load_config(config_path: Path) -> dict[str, Any]:
    """Load the datasources config; raise PilotGateError if missing or malformed."""
    if not config_path.exists():
        raise PilotGateError(f"config not found: {config_path}")
    try:
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise PilotGateError(f"config malformed: {e}") from e
    if "phase1_drugs" not in cfg or "bvbrc_ast" not in cfg:
        raise PilotGateError("config missing required keys: phase1_drugs, bvbrc_ast")
    return cfg


def fetch_bvbrc_drug_counts(
    drugs: tuple[str, ...],
    cfg: dict[str, Any],
    ast_tsv_path: Path | str | None = None,
) -> dict[str, int]:
    """Return raw drug-labeled isolate counts per drug.

    Resolution order:
    1. Explicit `ast_tsv_path` parameter (CLI --ast-tsv flag).
    2. `BVBRC_AST_TSV` environment variable.
    3. `cfg["bvbrc_ast"]["local_tsv_path"]` from datasources.yaml.
    4. Live BV-BRC API → NotImplementedError (deferred until first real-data run).

    Reuses `load_bvbrc_ast` from Step 5 — same broth-microdilution + organism
    filters apply, so the pilot's per-drug counts match what the downstream
    Step 5 ingestion will see.
    """
    resolved_path: Path | None = None

    if ast_tsv_path:
        resolved_path = Path(ast_tsv_path)
    elif env_path := os.environ.get("BVBRC_AST_TSV"):
        resolved_path = Path(env_path)
    elif (bvbrc_cfg := cfg.get("bvbrc_ast", {})) and bvbrc_cfg.get("local_tsv_path"):
        resolved_path = Path(bvbrc_cfg["local_tsv_path"])

    if resolved_path is None:
        raise NotImplementedError(
            "Live BV-BRC API integration deferred. Provide a pre-downloaded "
            "BV-BRC AST TSV via --ast-tsv <path>, BVBRC_AST_TSV env var, or "
            "config.bvbrc_ast.local_tsv_path. Download from ftp.bvbrc.org."
        )

    if not resolved_path.exists():
        raise PilotGateError(f"BV-BRC AST TSV not found at {resolved_path}")

    organism = cfg.get("bvbrc_ast", {}).get("default_filters", {}).get(
        "organism", DEFAULT_ORGANISM
    )
    method_filter_list = cfg.get("bvbrc_ast", {}).get("default_filters", {}).get(
        "measurement_method", ["broth_microdilution"]
    )
    method_filter = tuple(method_filter_list)

    ast = load_bvbrc_ast(resolved_path, organism=organism, method_filter=method_filter)

    counts: dict[str, int] = {}
    for drug in drugs:
        drug_lower = drug.lower()
        counts[drug] = int((ast["antibiotic"] == drug_lower).sum())
    return counts


def fetch_ncbi_assembly_quality(strain_ids: list[str]) -> dict[str, dict[str, int]]:
    """Fetch per-strain assembly contig_count + N50 from NCBI Datasets API.

    Returns mapping strain_id -> {"contig_count": int, "n50": int}.

    Scaffold; tests must monkeypatch with mocked metadata.
    """
    raise NotImplementedError(
        "NCBI assembly quality fetch is scaffolded; implement at first real-data run. "
        "Tests must monkeypatch this function with mocked metadata."
    )


def apply_method_filter(
    raw_counts: dict[str, int], method_filter_factor: float = 1.0
) -> dict[str, int]:
    """Pass-through when fetch_bvbrc_drug_counts already applied the filter.

    Default `method_filter_factor=1.0` means no further reduction. The local-TSV
    path applies broth-microdilution filtering inside `load_bvbrc_ast`, so counts
    coming out of `fetch_bvbrc_drug_counts` are already post-method-filter.

    A smaller factor (<1.0) is meaningful only for the not-yet-implemented live
    API path where the raw count would not yet be filtered.
    """
    return {drug: int(count * method_filter_factor) for drug, count in raw_counts.items()}


def apply_assembly_filter(
    method_filtered: dict[str, int], assembly_survival_factor: float = 0.85
) -> dict[str, int]:
    """Estimate post-assembly-quality-filter counts.

    Real implementation: join strains with NCBI assembly summaries, filter
    contig_count <= max + N50 >= min. For the pilot, empirical survival
    factor (typically ~0.85 for E. coli where most assemblies pass quality).
    """
    return {drug: int(count * assembly_survival_factor) for drug, count in method_filtered.items()}


def estimate_three_drug_intersection(
    per_drug_counts: dict[str, int], overlap_factor: float = 0.45
) -> int:
    """Estimate the count of strains labeled for all 3 drugs.

    Real implementation: set-intersect strain IDs across drugs. Pilot uses an
    overlap factor (typically ~0.45 for cipro+ceftriaxone+tet in E. coli
    based on BV-BRC empirical density).
    """
    if not per_drug_counts:
        return 0
    return int(min(per_drug_counts.values()) * overlap_factor)


def estimate_download_volume_gb(strain_count: int, avg_genome_mb: float = 5.0) -> float:
    """Rough download-volume estimate. E. coli genomes ~5 Mbp; FASTA + GFF3 + GenBank ~5MB/strain."""
    return (strain_count * avg_genome_mb) / 1024.0


def estimate_embedding_minutes(strain_count: int, gene_per_strain: int = 5000, minutes_per_million_seq: float = 30.0) -> float:
    """Rough embedding-compute estimate; original target was RTX 4090 + 4-bit Evo (still the formula's baseline). Actual hardware is GTX 860M (CC=5.0) → 4-bit Evo unreachable; NT v2 100M is the working model and runs ~30× slower than this estimate suggests. Treat output as optimistic upper bound on possible compute, not predicted wallclock."""
    total_sequences = strain_count * gene_per_strain
    return (total_sequences / 1_000_000) * minutes_per_million_seq


def run_pilot_gate(
    drugs: tuple[str, ...] | None = None,
    criteria: CohortSelectionCriteria | None = None,
    config_path: Path | str = DEFAULT_CONFIG_PATH,
    ast_tsv_path: Path | str | None = None,
) -> PilotReport:
    """Run the pilot gate.

    Returns a PilotReport with per-drug counts at each filter stage + GO/NO-GO
    verdict. The local-TSV path (preferred) is triggered when `ast_tsv_path`,
    `BVBRC_AST_TSV` env var, or `config.bvbrc_ast.local_tsv_path` is set.
    Otherwise raises NotImplementedError (live API deferred).
    """
    cfg = _load_config(Path(config_path))
    criteria = criteria or CohortSelectionCriteria()
    if drugs is None:
        drugs = tuple(d["name"] for d in cfg["phase1_drugs"])

    raw_counts = fetch_bvbrc_drug_counts(drugs, cfg, ast_tsv_path=ast_tsv_path)
    method_filtered = apply_method_filter(raw_counts)
    assembly_filtered = apply_assembly_filter(method_filtered)
    intersection = estimate_three_drug_intersection(assembly_filtered)

    per_drug = [
        PerDrugStageCounts(
            drug=d,
            raw=raw_counts.get(d, 0),
            after_method_filter=method_filtered.get(d, 0),
            after_assembly_filter=assembly_filtered.get(d, 0),
        )
        for d in drugs
    ]

    report = PilotReport(
        drugs=drugs,
        criteria=criteria,
        per_drug=per_drug,
        three_drug_intersection=intersection,
        estimated_download_gb=estimate_download_volume_gb(max(assembly_filtered.values()) if assembly_filtered else 0),
        estimated_embedding_minutes_rtx4090=estimate_embedding_minutes(
            max(assembly_filtered.values()) if assembly_filtered else 0
        ),
    )
    report.decide()
    return report


def write_pilot_report(report: PilotReport, output_path: Path | str) -> Path:
    """Write the pilot report as a human-readable markdown table."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("# Pilot Gate Report\n")
    lines.append(f"**Verdict:** `{report.go_no_go}`\n")
    lines.append(f"**Drugs:** {', '.join(report.drugs)}\n")
    lines.append(f"**Target per drug:** {report.criteria.target_per_drug}\n")
    lines.append(f"**Target 3-drug intersection:** {report.criteria.three_drug_intersection_target}\n")
    lines.append("")
    lines.append("## Per-drug counts by filter stage\n")
    lines.append("| Drug | Raw | After method filter | After assembly filter |")
    lines.append("|---|---|---|---|")
    for d in report.per_drug:
        lines.append(f"| {d.drug} | {d.raw} | {d.after_method_filter} | {d.after_assembly_filter} |")
    lines.append("")
    lines.append(f"## 3-drug intersection: {report.three_drug_intersection}\n")
    lines.append(f"## Estimated download volume: {report.estimated_download_gb:.1f} GB\n")
    lines.append(
        f"## Estimated embedding compute (formula baseline: RTX 4090 + 4-bit Evo; "
        f"actual hardware is GTX 860M + NT v2 100M, expect ~30× slower): "
        f"{report.estimated_embedding_minutes_rtx4090:.0f} min\n"
    )
    if report.failure_reasons:
        lines.append("## Failure reasons\n")
        for r in report.failure_reasons:
            lines.append(f"- {r}")
        lines.append("")
        lines.append("**Recommended remediation:** relax one filter / pick alternate drug / drop a drug from Phase 1.")

    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path
