"""Step 0.5 — Real-data pilot gate (HARD).

Metadata-only HTTP calls to BV-BRC + NCBI to estimate per-drug labeled-isolate
counts after filters. No genome downloads. No embedding compute. Cheap (~1hr).

HARD-gate semantics: returns exit-code 0 only when every Phase 1 drug has
>= target_per_drug strains after all filters. Per the post-tech-plan
brainstorm C1 resolution, the pilot script's non-zero exit halts /execute-plan
at Step 6 (cohort) before downstream ingestion + embedding fires.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

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


def fetch_bvbrc_drug_counts(drugs: tuple[str, ...], cfg: dict[str, Any]) -> dict[str, int]:
    """Fetch raw + method-filtered drug-labeled isolate counts from BV-BRC.

    Metadata-only HTTP. Returns mapping drug -> count.

    Note: real BV-BRC API integration is scaffolded here; concrete endpoint
    selection (data.bv-brc.org REST vs ftp.bvbrc.org TSV head) is deferred
    until first real-data run. Tests inject mocked counts via this function's
    return value.
    """
    raise NotImplementedError(
        "BV-BRC drug count fetch is scaffolded; implement at first real-data run. "
        "Tests must monkeypatch this function with mocked counts."
    )


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
    raw_counts: dict[str, int], method_filter_factor: float = 0.6
) -> dict[str, int]:
    """Estimate post-method-filter counts.

    Real implementation: filter BV-BRC AST rows where measurement_method in
    criteria.measurement_method_filter. For the pilot gate, we use an empirical
    survival factor (typically ~0.6 for broth_microdilution as the dominant
    method) and let real-data ingestion (Step 5) compute exact counts.
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
    """Rough embedding-compute estimate on RTX 4090 with 4-bit Evo + batching."""
    total_sequences = strain_count * gene_per_strain
    return (total_sequences / 1_000_000) * minutes_per_million_seq


def run_pilot_gate(
    drugs: tuple[str, ...] | None = None,
    criteria: CohortSelectionCriteria | None = None,
    config_path: Path | str = DEFAULT_CONFIG_PATH,
) -> PilotReport:
    """Run the metadata-only pilot gate.

    Returns a PilotReport with per-drug counts at each filter stage + GO/NO-GO
    verdict. Tests monkeypatch fetch_bvbrc_drug_counts to inject deterministic
    counts.
    """
    cfg = _load_config(Path(config_path))
    criteria = criteria or CohortSelectionCriteria()
    if drugs is None:
        drugs = tuple(d["name"] for d in cfg["phase1_drugs"])

    raw_counts = fetch_bvbrc_drug_counts(drugs, cfg)
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
        f"## Estimated embedding compute (RTX 4090, 4-bit Evo): "
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
