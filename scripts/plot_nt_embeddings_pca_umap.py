"""Diagnostic plot: PCA + UMAP of the 12-strain cipro NT embeddings.

Loads per-strain mean-pooled NT embeddings from the existing cache + cohort
parquet, projects to 2-D via PCA (top-2 PCs) and UMAP (n_neighbors=5, default
otherwise), and writes a 2-panel PNG with R/S labels + MLST clade overlay.

NOT a go/no-go gate (per LESSONS_LEARNED.md 2026-05-14: 2D projection at N=12
is informational, not decisive). Used to check whether separation tracks
resistance OR lineage — if it tracks lineage, the 12-strain smoke result is
lineage-confounded and Mash-clustering-stratified CV may be needed at Stage 2.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

from dna_decode.data.cohort import load_cohort
from dna_decode.models.cache import EmbeddingCache
from dna_decode.models.classifiers import aggregate_strain_features


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_mini_cohort.parquet"))
    parser.add_argument("--cache", type=Path, default=Path("data/processed/mini_cipro_nt_cache.h5"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("wiki/pca_umap_12strain_cipro_2026-05-14.png"),
    )
    args = parser.parse_args(argv)

    cohort = load_cohort(args.cohort)
    cache = EmbeddingCache(
        args.cache,
        model_name="nucleotide_transformer",
        model_version="InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
        embedding_dim=512,
    )

    drug_lower = args.drug.lower()
    X_rows: list[np.ndarray] = []
    labels: list[int] = []
    strain_ids: list[str] = []
    mlsts: list[str] = []
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        gene_ids = cache.list_genes(s.strain_id)
        if not gene_ids:
            continue
        gene_matrix = cache.bulk_get([(s.strain_id, g) for g in gene_ids])
        X_rows.append(aggregate_strain_features(gene_matrix, "mean"))
        labels.append(int(s.ast_labels[drug_lower]))
        strain_ids.append(s.strain_id)
        mlsts.append(str(getattr(s, "mlst", "unknown")))

    X = np.stack(X_rows)
    y = np.array(labels, dtype=int)
    print(f"[plot] N={len(X_rows)} strains, embedding shape={X.shape}, R/S = {(y==1).sum()}/{(y==0).sum()}")

    # PCA
    pca = PCA(n_components=2, random_state=42)
    Xp = pca.fit_transform(X)
    pca_var = pca.explained_variance_ratio_
    print(f"[plot] PCA explained variance: PC1={pca_var[0]:.3f}, PC2={pca_var[1]:.3f}")

    # UMAP — optional dep; degrade gracefully if missing
    Xu: np.ndarray | None = None
    try:
        import umap

        reducer = umap.UMAP(n_neighbors=5, min_dist=0.1, random_state=42)
        Xu = reducer.fit_transform(X)
        print(f"[plot] UMAP fit complete (n_neighbors=5)")
    except ImportError:
        print("[plot] umap-learn not installed — PCA-only output")

    # Plot
    fig, axes = plt.subplots(1, 2 if Xu is not None else 1, figsize=(12, 5) if Xu is not None else (6, 5))
    if Xu is None:
        axes = [axes]

    def _scatter(ax, X2, title, xlabel, ylabel):
        # Color by R/S, mark MLST with text annotation
        for cls, cls_name, marker in [(1, "Resistant", "o"), (0, "Susceptible", "^")]:
            mask = y == cls
            ax.scatter(X2[mask, 0], X2[mask, 1], s=80, alpha=0.7, marker=marker, label=cls_name)
        for i, (sid, mlst) in enumerate(zip(strain_ids, mlsts)):
            ax.annotate(
                mlst,
                (X2[i, 0], X2[i, 1]),
                fontsize=7,
                alpha=0.6,
                xytext=(5, 5),
                textcoords="offset points",
            )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.3)

    _scatter(
        axes[0],
        Xp,
        f"PCA — 12-strain cipro NT embeddings (PC1 {pca_var[0]:.1%}, PC2 {pca_var[1]:.1%})",
        f"PC1 ({pca_var[0]:.1%} var)",
        f"PC2 ({pca_var[1]:.1%} var)",
    )
    if Xu is not None:
        _scatter(axes[1], Xu, "UMAP (n_neighbors=5)", "UMAP1", "UMAP2")

    fig.suptitle(
        "PCA + UMAP of mean-pooled NT-v2-100M embeddings — 12-strain cipro mini cohort\n"
        "Diagnostic only; not a go/no-go gate at N=12. MLST shown as text label.",
        fontsize=10,
    )
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150)
    print(f"[plot] wrote {args.output}")

    # Quick textual diagnostic: does PC1 sign correlate with R/S?
    pc1_R_mean = float(Xp[y == 1, 0].mean())
    pc1_S_mean = float(Xp[y == 0, 0].mean())
    print(f"[plot] PC1 mean(R) = {pc1_R_mean:+.3f}  vs  mean(S) = {pc1_S_mean:+.3f}  (separation = {abs(pc1_R_mean - pc1_S_mean):.3f})")
    # MLST cardinality
    unique_mlst = len(set(mlsts))
    print(f"[plot] unique MLSTs = {unique_mlst} of {len(mlsts)} strains")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
