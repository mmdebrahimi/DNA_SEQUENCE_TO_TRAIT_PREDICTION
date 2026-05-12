"""Step 15 — End-to-end smoke pipeline on a synthetic 5-strain corpus.

Exercises the full Wave 1-3 contract: ingest synthetic fixtures → build
cohort → populate embedding cache via MockFoundationModel → train XGBoost
classifier → leave-one-strain-out CV → ISM → tier classification → write
report. Should complete in <60s on CPU; no GPU, no foundation-model
downloads, no network.

Used both as a regression test for the pipeline + as user-facing
documentation of the end-to-end flow.

Synthetic resistance signal: gene `g1` (the "synthetic gyrase A" locus)
carries a 50-bp consensus seeded into RESISTANT strains' embedding via
the MockFoundationModel hash function. SUSCEPTIBLE strains lack the
consensus. The classifier + ISM + tier classification should all surface
g1 / TAG_001 as the top attribution hit.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# Quiet expected sklearn warnings on tiny synthetic data
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="dna_decode")


from dna_decode.data.annotations import parse_gff3
from dna_decode.data.cohort import (
    CandidateStrain,
    CohortSelectionCriteria,
    StrainCohort,
    build_cohort,
)
from dna_decode.data.resistance_db import ResistanceCatalog, ResistanceEntry
from dna_decode.eval.cv import leave_one_strain_out_cv
from dna_decode.eval.metrics import compute_metrics
from dna_decode.interp.mutagenesis import (
    build_attribution_report,
    gene_level_mutagenesis,
)
from dna_decode.models.classifiers import (
    aggregate_strain_features,
    predict_proba,
    train_xgboost_classifier,
)
from dna_decode.models.foundation import MockFoundationModel, ModelMetadata


N_STRAINS = 10  # 5 resistant + 5 susceptible
N_GENES = 10
EMBEDDING_DIM = 16
SIGNAL_GENE = "g1"
SIGNAL_LOCUS_TAG = "TAG_001"


def _build_synthetic_strains(
    fixtures_dir: Path,
) -> tuple[list[CandidateStrain], dict[str, dict[str, np.ndarray]], dict[str, int]]:
    """Generate 10 synthetic strains with seeded resistance signal in g1."""
    annotations = parse_gff3(fixtures_dir / "annotations.gff3")
    gene_ids = annotations[annotations["type"] == "CDS"]["gene_id"].tolist()

    mock = MockFoundationModel(
        ModelMetadata(
            name="mock",
            huggingface_id="mock://smoke",
            embedding_dim=EMBEDDING_DIM,
            max_context=10_000,
        )
    )

    strain_embeddings: dict[str, dict[str, np.ndarray]] = {}
    candidates: list[CandidateStrain] = []
    labels: dict[str, int] = {}

    for i in range(N_STRAINS):
        is_resistant = i < N_STRAINS // 2
        strain_id = f"synth_{i:03d}"

        per_gene_embeddings: dict[str, np.ndarray] = {}
        for gene_id in gene_ids:
            if gene_id == SIGNAL_GENE and is_resistant:
                signal_seq = f"resistant_signal_for_{gene_id}_in_{strain_id}"
            elif gene_id == SIGNAL_GENE:
                signal_seq = f"susceptible_baseline_for_{gene_id}_in_{strain_id}"
            else:
                signal_seq = f"background_for_{gene_id}_in_{strain_id}"
            per_gene_embeddings[gene_id] = mock._embed_window(signal_seq)

        # Reinforce signal: resistant strains' g1 embedding has higher dim-0 value
        if is_resistant:
            per_gene_embeddings[SIGNAL_GENE] = per_gene_embeddings[SIGNAL_GENE].copy()
            per_gene_embeddings[SIGNAL_GENE][0] += 5.0
        strain_embeddings[strain_id] = per_gene_embeddings

        candidates.append(
            CandidateStrain(
                strain_id=strain_id,
                assembly_accession=f"SYNTH_{i:03d}.1",
                mlst=f"ST{i % 3}",  # 3 synthetic MLSTs across 10 strains
                contig_count=10,
                n50=200_000,
                ast_labels={"ciprofloxacin": 1 if is_resistant else 0},
            )
        )
        labels[strain_id] = 1 if is_resistant else 0

    return candidates, strain_embeddings, labels


def _pool_strain_features(
    strain_embeddings: dict[str, dict[str, np.ndarray]],
    strain_ids: list[str],
    gene_ids: list[str],
) -> np.ndarray:
    """Build (n_strains, embedding_dim) mean-pooled feature matrix."""
    out = np.empty((len(strain_ids), EMBEDDING_DIM), dtype=np.float32)
    for i, sid in enumerate(strain_ids):
        gene_matrix = np.stack([strain_embeddings[sid][g] for g in gene_ids])
        out[i] = aggregate_strain_features(gene_matrix, "mean")
    return out


def run_smoke(fixtures_dir: Path, output_dir: Path) -> dict[str, object]:
    """End-to-end smoke run. Returns a results dict."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[smoke] loading fixture annotations from {fixtures_dir}/annotations.gff3")
    annotations = parse_gff3(fixtures_dir / "annotations.gff3")
    gene_ids = annotations[annotations["type"] == "CDS"]["gene_id"].tolist()
    print(f"[smoke] {len(gene_ids)} genes in synthetic annotation")

    print(f"[smoke] generating {N_STRAINS} synthetic strains with seeded signal in {SIGNAL_GENE}")
    candidates, strain_embeddings, labels_dict = _build_synthetic_strains(fixtures_dir)

    print(f"[smoke] building cohort (criteria: target_per_drug=4, intersection=2)")
    cohort = build_cohort(
        candidates,
        ("ciprofloxacin",),
        CohortSelectionCriteria(
            target_per_drug=4,
            three_drug_intersection_target=2,
            assembly_contig_count_max=500,
            assembly_n50_min=50_000,
        ),
    )
    cohort_strain_ids = [s.strain_id for s in cohort.strains]
    print(f"[smoke] cohort size: {len(cohort_strain_ids)} strains")

    # Pool features
    print(f"[smoke] mean-pooling per-gene embeddings → strain features")
    X = _pool_strain_features(strain_embeddings, cohort_strain_ids, gene_ids)
    y = np.array([labels_dict[sid] for sid in cohort_strain_ids], dtype=int)

    # LOSO CV
    print(f"[smoke] running leave-one-strain-out CV with XGBoost classifier")

    def _train(X_train: np.ndarray, y_train: np.ndarray):
        return train_xgboost_classifier(
            X_train, y_train, drug_name="ciprofloxacin", calibrate=False
        )

    def _predict(model, X_test: np.ndarray) -> np.ndarray:
        return predict_proba(model, X_test)

    cv_result = leave_one_strain_out_cv(
        X, y, cohort_strain_ids, _train, _predict, drug="ciprofloxacin"
    )

    metrics = compute_metrics(cv_result.all_y_true, cv_result.all_y_score)
    print(f"[smoke] LOSO CV AUROC: {metrics.auroc:.3f} (n={metrics.n_samples})")

    # Train on full data for attribution
    print(f"[smoke] training final classifier on full synthetic cohort")
    full_clf = train_xgboost_classifier(X, y, drug_name="ciprofloxacin", calibrate=False)

    # ISM on one resistant strain
    target_strain = next(sid for sid in cohort_strain_ids if labels_dict[sid] == 1)
    print(f"[smoke] gene-level ISM on resistant strain {target_strain}")
    gene_effects = gene_level_mutagenesis(
        full_clf, strain_embeddings[target_strain], annotations
    )
    print(f"[smoke] top-3 attributed genes by |delta|:")
    for _, row in gene_effects.head(3).iterrows():
        print(
            f"        {row['gene_id']:>6s}  ({row['locus_tag']:>10s})  "
            f"delta={row['prediction_delta']:+.3f}"
        )

    # Tier classification
    cat = ResistanceCatalog()
    cat.add(
        ResistanceEntry(
            "TAG_001",
            "synthetic-resistance",
            "fluoroquinolone",
            "synthetic",
            "smoke",
            "SYNTH_001",
        )
    )
    report = build_attribution_report(
        gene_effects, "fluoroquinolone", cat, top_k=5, annotations=annotations
    )
    print(f"[smoke] attribution tier counts: {report.tier_counts()}")
    print(f"[smoke] fraction Tier 1-3: {report.fraction_tier_1_to_3():.2f}")

    # Write smoke report
    report_path = output_dir / "smoke_report.md"
    lines = [
        "# Smoke Pipeline Report",
        "",
        f"**Cohort size:** {len(cohort_strain_ids)} strains",
        f"**LOSO AUROC:** {metrics.auroc:.3f} (n={metrics.n_samples})",
        f"**Top-3 attributed genes:**",
    ]
    for _, row in gene_effects.head(3).iterrows():
        lines.append(
            f"- {row['gene_id']} ({row['locus_tag']}) — delta={row['prediction_delta']:+.3f}"
        )
    lines += [
        "",
        f"**Attribution tier counts (top-5):** {report.tier_counts()}",
        f"**Fraction Tier 1-3:** {report.fraction_tier_1_to_3():.2f}",
        f"**Fraction Fail:** {report.fraction_fail():.2f}",
        "",
        f"**Expected:** top-1 attribution gene = `{SIGNAL_GENE}` (locus `{SIGNAL_LOCUS_TAG}`)",
        f"**Actual top-1:** `{gene_effects.iloc[0]['gene_id']}` (locus `{gene_effects.iloc[0]['locus_tag']}`)",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"[smoke] report written: {report_path}")

    return {
        "auroc": float(metrics.auroc) if not np.isnan(metrics.auroc) else None,
        "n_samples": int(metrics.n_samples),
        "top_1_gene": gene_effects.iloc[0]["gene_id"],
        "top_1_locus_tag": gene_effects.iloc[0]["locus_tag"],
        "tier_1_3_fraction": report.fraction_tier_1_to_3(),
        "fail_fraction": report.fraction_fail(),
        "report_path": str(report_path),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase 1 smoke pipeline — end-to-end on synthetic 5-strain corpus."
    )
    parser.add_argument(
        "--fixtures-dir",
        default="tests/fixtures/ecoli_mini",
        help="Path to the synthetic fixture corpus (default: tests/fixtures/ecoli_mini).",
    )
    parser.add_argument(
        "--output-dir",
        default="data/processed",
        help="Where to write the smoke report (default: data/processed).",
    )
    args = parser.parse_args(argv)

    try:
        results = run_smoke(Path(args.fixtures_dir), Path(args.output_dir))
    except Exception as e:
        print(f"[smoke] FAILED: {type(e).__name__}: {e}", file=sys.stderr)
        return 1

    # Sanity: top-1 should be the seeded signal gene
    if results["top_1_gene"] != SIGNAL_GENE:
        print(
            f"[smoke] WARNING: top-1 attribution is {results['top_1_gene']!r}, "
            f"expected {SIGNAL_GENE!r}",
            file=sys.stderr,
        )
        return 2

    auroc = results["auroc"]
    if auroc is None or auroc < 0.85:
        print(
            f"[smoke] WARNING: LOSO AUROC {auroc} below 0.85 threshold on seeded signal",
            file=sys.stderr,
        )
        return 3

    print(f"[smoke] PASS: AUROC={auroc:.3f}, top-1={results['top_1_gene']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
