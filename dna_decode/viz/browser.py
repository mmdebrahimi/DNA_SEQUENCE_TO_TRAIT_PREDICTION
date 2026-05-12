"""Step 13 — Genome-browser-lite visualization (matplotlib + TSV export).

Phase 1 ships matplotlib + TSV (per ship-path plan D4); pygenometracks adapter
deferred to Phase 2. Outputs:
- TSV exports of `GeneEffectTable` and `PositionEffectTable` for downstream
  tooling (load into IGV, R, Excel).
- matplotlib PNG: per-gene bar of `prediction_delta` (gene-level ISM) +
  per-position line plot of `|prediction_delta|` (saturation mutagenesis).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def export_attribution_tsv(
    gene_effects: pd.DataFrame,
    output_path: Path | str,
    position_effects: pd.DataFrame | None = None,
    position_output_path: Path | str | None = None,
) -> tuple[Path, Path | None]:
    """Dump gene-level + (optional) saturation tables as TSV files.

    Args:
        gene_effects: `GeneEffectTable` from Step 11's `gene_level_mutagenesis`.
        output_path: TSV output for gene_effects.
        position_effects: optional `PositionEffectTable` from saturation_mutagenesis.
        position_output_path: required when position_effects is provided.

    Returns:
        (gene_tsv_path, position_tsv_path_or_None)
    """
    gene_path = Path(output_path)
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    gene_effects.to_csv(gene_path, sep="\t", index=False)

    pos_path: Path | None = None
    if position_effects is not None:
        if position_output_path is None:
            raise ValueError("position_output_path required when position_effects supplied")
        pos_path = Path(position_output_path)
        pos_path.parent.mkdir(parents=True, exist_ok=True)
        position_effects.to_csv(pos_path, sep="\t", index=False)

    return gene_path, pos_path


def render_attribution_plot(
    gene_effects: pd.DataFrame,
    output_path: Path | str,
    top_k: int = 20,
    drug_name: str = "",
    position_effects: pd.DataFrame | None = None,
) -> Path:
    """Render a matplotlib PNG of attribution results.

    Top panel: bar chart of top-K genes by |prediction_delta|.
    Bottom panel (optional): per-position |prediction_delta| line plot for
    the top-1 gene from `gene_effects`, sourced from `position_effects`.

    Args:
        gene_effects: GeneEffectTable from Step 11.
        output_path: PNG output path.
        top_k: number of top genes to visualize.
        drug_name: optional title annotation.
        position_effects: optional PositionEffectTable for the bottom panel.

    Returns:
        Path to the written PNG.
    """
    import matplotlib

    matplotlib.use("Agg")  # headless backend; no display required
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    has_position_panel = position_effects is not None and len(position_effects) > 0
    n_panels = 2 if has_position_panel else 1
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 4 * n_panels), squeeze=False)

    # Top panel — gene-level bar chart
    ax_top = axes[0, 0]
    top_rows = gene_effects.head(top_k).copy()
    labels = [
        f"{row['gene_id']}\n({row.get('locus_tag', '')})" if row.get("locus_tag") else row["gene_id"]
        for _, row in top_rows.iterrows()
    ]
    deltas = top_rows["prediction_delta"].abs().values
    bar_colors = ["#d62728" if d > 0 else "#1f77b4" for d in top_rows["prediction_delta"].values]
    ax_top.bar(range(len(top_rows)), deltas, color=bar_colors)
    ax_top.set_xticks(range(len(top_rows)))
    ax_top.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax_top.set_ylabel("|prediction_delta|")
    title = f"Top-{top_k} gene-level attribution"
    if drug_name:
        title += f" — {drug_name}"
    ax_top.set_title(title)
    ax_top.grid(axis="y", alpha=0.3)

    # Bottom panel — saturation mutagenesis on top-1 gene
    if has_position_panel:
        ax_bot = axes[1, 0]
        top_gene_id = top_rows.iloc[0]["gene_id"]
        gene_positions = position_effects[position_effects["gene_id"] == top_gene_id]
        if len(gene_positions) > 0:
            # Aggregate per position: max-abs delta across alt bases
            agg = (
                gene_positions.groupby("position")["prediction_delta"]
                .apply(lambda s: s.abs().max())
                .reset_index()
                .sort_values("position")
            )
            ax_bot.plot(agg["position"], agg["prediction_delta"], color="#2ca02c", linewidth=1.2)
            ax_bot.fill_between(
                agg["position"], 0, agg["prediction_delta"], color="#2ca02c", alpha=0.3
            )
            ax_bot.set_xlabel(f"Position within {top_gene_id}")
            ax_bot.set_ylabel("|max prediction_delta|")
            ax_bot.set_title(f"Saturation-mutagenesis attribution map — {top_gene_id}")
            ax_bot.grid(alpha=0.3)
        else:
            ax_bot.text(
                0.5,
                0.5,
                f"No saturation data for {top_gene_id}",
                ha="center",
                va="center",
                transform=ax_bot.transAxes,
            )

    plt.tight_layout()
    plt.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)

    return output_path
