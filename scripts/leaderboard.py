"""Step 17 — Comparative model leaderboard (shell-loop style).

Per ship-path plan D5: NOT a dedicated `dna_decode/eval/leaderboard.py` module.
Loops over the configured foundation models + classical baselines, invokes
`pipeline.py train` per (model × drug), aggregates results into one
markdown report.

Phase 1 default: Evo + DNABERT-2 only (NT + GENA-LM added incrementally).
Classical baselines (AMRFinder + k-mer + gene-presence) run alongside via
their own training paths.

Output: `data/processed/leaderboard.md` with per-drug × per-model rows.
"""
from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class LeaderboardRow:
    """One row of the leaderboard — drug × model combination."""

    drug: str
    model_name: str
    is_classical: bool
    auroc: float | None = None
    auprc: float | None = None
    clade_baseline_gap: float | None = None
    attribution_tier_1_3: float | None = None
    notes: str = ""


@dataclass
class LeaderboardReport:
    """Aggregate leaderboard output."""

    rows: list[LeaderboardRow] = field(default_factory=list)

    def as_markdown(self) -> str:
        lines: list[str] = ["# Phase 1 Leaderboard", ""]
        # Group by drug for readability
        drugs = sorted({r.drug for r in self.rows})
        for drug in drugs:
            lines.append(f"## {drug}")
            lines.append("")
            lines.append("| Model | Type | AUROC | AUPRC | Clade-gap | Tier 1-3 frac | Notes |")
            lines.append("|---|---|---|---|---|---|---|")
            drug_rows = [r for r in self.rows if r.drug == drug]
            drug_rows.sort(key=lambda r: (r.auroc is None, -(r.auroc or 0)))
            for r in drug_rows:
                lines.append(
                    f"| {r.model_name} | {'classical' if r.is_classical else 'foundation'} | "
                    f"{_fmt(r.auroc)} | {_fmt(r.auprc)} | {_fmt(r.clade_baseline_gap, signed=True)} | "
                    f"{_fmt(r.attribution_tier_1_3)} | {r.notes} |"
                )
            lines.append("")
        return "\n".join(lines)


def _fmt(value: float | None, signed: bool = False) -> str:
    if value is None:
        return "—"
    if signed:
        return f"{value:+.3f}"
    return f"{value:.3f}"


def run_foundation_model_for_drug(
    model_name: str,
    drug: str,
    project_root: Path,
    extra_args: list[str] | None = None,
) -> LeaderboardRow:
    """Invoke `pipeline.py train` for one (model, drug) combination.

    Reads the pickled bundle written by pipeline.py train to extract metrics.
    """
    cmd = [
        sys.executable,
        "-m",
        "scripts.pipeline",
        "train",
        "--drug",
        drug,
        "--model",
        model_name,
        "--include-clade-baseline",
    ]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, cwd=str(project_root), capture_output=True, text=True)
    if result.returncode != 0:
        return LeaderboardRow(
            drug=drug,
            model_name=model_name,
            is_classical=False,
            notes=f"train failed: exit {result.returncode}; {result.stderr.strip()[:120]}",
        )

    # Read the trained model bundle for canonical metrics
    bundle_path = (
        project_root / "data" / "processed" / "models" / f"{drug}_{model_name}.pkl"
    )
    if not bundle_path.exists():
        return LeaderboardRow(
            drug=drug,
            model_name=model_name,
            is_classical=False,
            notes=f"bundle not found at {bundle_path}",
        )
    with open(bundle_path, "rb") as f:
        bundle = pickle.load(f)

    return LeaderboardRow(
        drug=drug,
        model_name=model_name,
        is_classical=False,
        auroc=bundle.get("auroc_loso"),
    )


def run_classical_baseline_for_drug(
    baseline_name: str, drug: str, project_root: Path
) -> LeaderboardRow:
    """Placeholder classical-baseline invocation.

    Phase 1 wires the three classical baseline trainers (Step 18) but doesn't
    expose them via the pipeline.py CLI yet. Leaderboard rows for classical
    baselines are filled in when callers run them externally + drop pickled
    bundles in the same `data/processed/models/<drug>_<baseline>.pkl` shape.

    This function reads the bundle if it exists; otherwise marks the row
    as 'pending'.
    """
    bundle_path = (
        project_root / "data" / "processed" / "models" / f"{drug}_{baseline_name}.pkl"
    )
    if not bundle_path.exists():
        return LeaderboardRow(
            drug=drug,
            model_name=baseline_name,
            is_classical=True,
            notes="pending — run classical baseline trainer externally + drop bundle",
        )
    with open(bundle_path, "rb") as f:
        bundle = pickle.load(f)
    return LeaderboardRow(
        drug=drug,
        model_name=baseline_name,
        is_classical=True,
        auroc=bundle.get("auroc_loso"),
    )


def run_leaderboard(
    drugs: Iterable[str],
    foundation_models: Iterable[str],
    classical_baselines: Iterable[str],
    project_root: Path,
    output_path: Path,
) -> LeaderboardReport:
    """Iterate over (drug, model) pairs + build the leaderboard report."""
    report = LeaderboardReport()

    for drug in drugs:
        print(f"\n[leaderboard] === drug: {drug} ===")
        for model in foundation_models:
            print(f"[leaderboard]  foundation: {model}")
            row = run_foundation_model_for_drug(model, drug, project_root)
            report.rows.append(row)
            auroc_str = f"{row.auroc:.3f}" if row.auroc is not None else "—"
            print(f"[leaderboard]    AUROC: {auroc_str} {row.notes}")
        for baseline in classical_baselines:
            print(f"[leaderboard]  classical: {baseline}")
            row = run_classical_baseline_for_drug(baseline, drug, project_root)
            report.rows.append(row)
            auroc_str = f"{row.auroc:.3f}" if row.auroc is not None else "—"
            print(f"[leaderboard]    AUROC: {auroc_str} {row.notes}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report.as_markdown(), encoding="utf-8")
    print(f"\n[leaderboard] report written: {output_path}")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase 1 leaderboard — fan out pipeline.py train over models × drugs"
    )
    parser.add_argument(
        "--drugs",
        default="ciprofloxacin,ceftriaxone,tetracycline",
        help="Comma-separated drugs (default: Phase 1 trio).",
    )
    parser.add_argument(
        "--models",
        default="evo,dnabert2",
        help="Comma-separated foundation models (default: evo + dnabert2; "
        "NT + GENA-LM added incrementally).",
    )
    parser.add_argument(
        "--classical",
        default="classical_amrfinder,classical_kmer_logreg,classical_gene_presence",
        help="Comma-separated classical baselines (matched against bundles in "
        "data/processed/models/<drug>_<baseline>.pkl).",
    )
    parser.add_argument(
        "--output",
        default="data/processed/leaderboard.md",
        help="Where to write the leaderboard markdown report.",
    )
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root path (where scripts/pipeline.py lives).",
    )
    args = parser.parse_args(argv)

    drugs = [d.strip() for d in args.drugs.split(",") if d.strip()]
    foundation_models = [m.strip() for m in args.models.split(",") if m.strip()]
    classical_baselines = [c.strip() for c in args.classical.split(",") if c.strip()]
    project_root = Path(args.project_root).resolve()

    report = run_leaderboard(
        drugs=drugs,
        foundation_models=foundation_models,
        classical_baselines=classical_baselines,
        project_root=project_root,
        output_path=Path(args.output),
    )

    # Exit 1 if every row is missing AUROC (run produced no useful results)
    if not any(r.auroc is not None for r in report.rows):
        print("[leaderboard] no AUROC values populated; check pipeline.py train outputs", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
