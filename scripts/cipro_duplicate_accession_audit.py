"""Audit duplicate assembly_accession rows in the current cipro cohort.

Purpose:
  - detect same-assembly duplicates inside the ciprofloxacin pool
  - flag potential LOSO leakage because training currently folds by strain_id
  - optionally inspect a trained model bundle for provenance/context
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

from dna_decode.data.cohort import find_duplicate_accessions, load_cohort


def _render_markdown(report: dict) -> str:
    lines = [
        "# Cipro duplicate accession audit",
        "",
        f"- Cohort file: `{report['cohort_path']}`",
        f"- Drug: `{report['drug']}`",
        f"- Pool size: `{report['n_pool_strains']}`",
        f"- Duplicate accessions in pool: `{report['n_duplicate_accessions']}`",
        f"- Strains participating in duplicate accessions: `{report['n_duplicate_strains']}`",
        f"- LOSO leakage present: `{report['loso_leakage_present']}`",
        "",
    ]
    if report.get("model_bundle"):
        model = report["model_bundle"]
        lines.extend(
            [
                "## Model bundle",
                "",
                f"- Model path: `{model['model_path']}`",
                f"- Drug: `{model.get('drug')}`",
                f"- Model name: `{model.get('model_name')}`",
                f"- CV strategy: `{model.get('cv_strategy')}`",
                f"- CV grouping: `{model.get('cv_grouping')}`",
                f"- Training cohort: `{model.get('training_cohort')}`",
                f"- Trained on: `{model.get('trained_on')}`",
                f"- Primary CV AUROC: `{model.get('auroc_cv_primary')}`",
                f"- Legacy auroc_loso field: `{model.get('auroc_loso')}`",
                f"- n_strains: `{model.get('n_strains')}`",
                "",
            ]
        )

    lines.extend(["## Duplicate accessions", ""])
    if not report["duplicates"]:
        lines.extend(["(none)", ""])
    else:
        lines.extend(["| Assembly accession | Strain IDs |", "|---|---|"])
        for row in report["duplicates"]:
            lines.append(f"| {row['assembly_accession']} | {', '.join(row['strain_ids'])} |")
        lines.append("")

    lines.extend(
        [
            "## Verdict",
            "",
            report["verdict"],
            "",
        ]
    )
    return "\n".join(lines)


def build_report(cohort_path: Path, drug: str, model_path: Path | None = None) -> dict:
    cohort = load_cohort(cohort_path)
    drug_lower = drug.lower()
    if drug_lower not in cohort.per_drug_strain_ids:
        raise ValueError(
            f"Drug {drug!r} not present in cohort per_drug_strain_ids: "
            f"{sorted(cohort.per_drug_strain_ids)}"
        )

    pool_ids = cohort.per_drug_strain_ids[drug_lower]
    duplicate_map = find_duplicate_accessions(cohort, restrict_to_strain_ids=pool_ids)
    duplicates = [
        {"assembly_accession": accession, "strain_ids": strain_ids}
        for accession, strain_ids in sorted(duplicate_map.items())
    ]

    model_bundle = None
    if model_path is not None:
        with open(model_path, "rb") as f:
            bundle = pickle.load(f)
        model_bundle = {
            "model_path": str(model_path),
            "drug": bundle.get("drug"),
            "model_name": bundle.get("model_name"),
            "cv_strategy": bundle.get("cv_strategy"),
            "cv_grouping": bundle.get("cv_grouping"),
            "training_cohort": bundle.get("training_cohort"),
            "trained_on": bundle.get("trained_on"),
            "auroc_cv_primary": bundle.get("auroc_cv_primary", bundle.get("auroc_loso")),
            "auroc_loso": bundle.get("auroc_loso"),
            "n_strains": bundle.get("n_strains"),
        }

    loso_leakage_present = bool(duplicates)
    model_bundle_is_accession_safe = bool(
        model_bundle
        and (
            model_bundle.get("cv_strategy") == "leave_one_accession_out"
            or model_bundle.get("cv_grouping") == "assembly_accession"
        )
    )
    if loso_leakage_present and model_bundle_is_accession_safe:
        verdict = (
            "PASS_WITH_MITIGATION: duplicate assembly_accession rows exist in the cipro "
            "cohort, but the inspected model bundle uses accession-safe CV "
            "(leave_one_accession_out / assembly_accession grouping)."
        )
    elif loso_leakage_present:
        verdict = (
            "FAIL: same assembly_accession appears under multiple strain_ids inside the "
            "ciprofloxacin pool. LOSO by strain_id can leak the same assembly across "
            "train and held-out folds."
        )
    else:
        verdict = "PASS: no duplicated non-empty assembly_accession values found in the cipro pool."

    return {
        "cohort_path": str(cohort_path),
        "drug": drug_lower,
        "n_pool_strains": len(pool_ids),
        "n_duplicate_accessions": len(duplicates),
        "n_duplicate_strains": sum(len(row["strain_ids"]) for row in duplicates),
        "loso_leakage_present": loso_leakage_present,
        "model_bundle_is_accession_safe": model_bundle_is_accession_safe,
        "duplicates": duplicates,
        "model_bundle": model_bundle,
        "verdict": verdict,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True, help="Path to cohort parquet")
    parser.add_argument("--drug", default="ciprofloxacin", help="Drug pool to audit")
    parser.add_argument("--model-path", help="Optional trained model pickle to inspect")
    parser.add_argument("--output", required=True, help="Markdown output path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cohort_path = Path(args.cohort)
    if not cohort_path.exists():
        print(f"[dup-audit] cohort not found at {cohort_path}", file=sys.stderr)
        return 2

    model_path = Path(args.model_path) if args.model_path else None
    if model_path is not None and not model_path.exists():
        print(f"[dup-audit] model not found at {model_path}", file=sys.stderr)
        return 2

    try:
        report = build_report(cohort_path, args.drug, model_path)
    except Exception as e:  # pragma: no cover - CLI guard
        print(f"[dup-audit] {e}", file=sys.stderr)
        return 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(_render_markdown(report), encoding="utf-8")
    json_path = output_path.with_suffix(".json")
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[dup-audit] wrote {output_path}")
    print(f"[dup-audit] wrote {json_path}")
    if report["loso_leakage_present"] and not report["model_bundle_is_accession_safe"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
