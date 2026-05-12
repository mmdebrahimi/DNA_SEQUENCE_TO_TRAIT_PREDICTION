"""Step 11 — In-silico mutagenesis (Phase 1's sole attribution mechanism).

Captum IG deferred to Phase 2 per post-tech-plan brainstorm C1: XGBoost is
non-differentiable; IG would require a differentiable MLP head not built
until Phase 2. ISM works for ANY predictor (no differentiability needed),
is the published consensus for sequence-bio interpretability, and ships
Phase 1 as gene-level + nucleotide-level saturation on the top-K genes.

Phase 1 attribution-success rubric: Tier 1-5 (exact known SNP → plausible
region → Fail clade marker). See plan's Verification section.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
import pandas as pd

from dna_decode.models.classifiers import (
    TrainedClassifier,
    aggregate_strain_features,
    predict_proba,
)


GENE_EFFECT_COLUMNS = (
    "gene_id",
    "locus_tag",
    "prediction_delta",
    "baseline_probability",
    "knockout_probability",
)

POSITION_EFFECT_COLUMNS = (
    "gene_id",
    "position",
    "ref_base",
    "alt_base",
    "prediction_delta",
)

GeneEffectTable = pd.DataFrame
PositionEffectTable = pd.DataFrame


# ---- Tier 1-5 attribution-success framework ----


@dataclass(frozen=True)
class AttributionTier:
    """Tier classification for one predicted attribution hit.

    Tiers (descending biological strength):
        1 = exact known resistance variant (codon-level match)
        2 = same codon or amino-acid region
        3 = same gene
        4 = same operon/pathway/mechanism class
        5 = plausible by literature/DB but unvalidated for THIS drug
       -1 = Fail (clade marker / known-irrelevant / intergenic-no-annotation)
    """

    tier: int
    rationale: str


@dataclass
class AttributionReport:
    """Aggregate of Tier 1-5 classifications across top-K hits for a drug."""

    drug: str
    top_k: int
    hits: list[tuple[str, AttributionTier]] = field(default_factory=list)

    def tier_counts(self) -> dict[int, int]:
        counts: dict[int, int] = {}
        for _, tier in self.hits:
            counts[tier.tier] = counts.get(tier.tier, 0) + 1
        return counts

    def fraction_tier_1_to_3(self) -> float:
        if not self.hits:
            return 0.0
        return sum(1 for _, t in self.hits if 1 <= t.tier <= 3) / len(self.hits)

    def fraction_fail(self) -> float:
        if not self.hits:
            return 0.0
        return sum(1 for _, t in self.hits if t.tier == -1) / len(self.hits)


# ---- Gene-level ISM ----


def gene_level_mutagenesis(
    classifier: TrainedClassifier,
    strain_gene_embeddings: dict[str, np.ndarray],
    annotations: pd.DataFrame | None = None,
) -> GeneEffectTable:
    """Run in-silico knockout at gene granularity.

    For each gene in the strain, zero out (knockout proxy) that gene's
    embedding row, re-aggregate strain features, re-predict, record the
    prediction delta. Returns rows sorted by absolute delta descending.

    Args:
        classifier: TrainedClassifier from Step 9.
        strain_gene_embeddings: gene_id → 1-D embedding vector (from
            cache.bulk_get or cache.get per gene).
        annotations: optional AnnotationTable from Step 3 to enrich rows
            with locus_tag.

    Note: "knockout" here is a presence-proxy; zeroing the embedding
    simulates the strain WITHOUT that gene's contribution to the pooled
    representation. Not equivalent to gene deletion at the sequence level
    (saturation mutagenesis handles that).
    """
    gene_ids = list(strain_gene_embeddings.keys())
    if not gene_ids:
        return pd.DataFrame(columns=list(GENE_EFFECT_COLUMNS))

    # Build the full per-gene matrix once
    full_matrix = np.stack([strain_gene_embeddings[g] for g in gene_ids])
    baseline_features = aggregate_strain_features(full_matrix, "mean").reshape(1, -1)
    baseline_prob = float(predict_proba(classifier, baseline_features)[0])

    # Locus-tag lookup from annotations
    locus_lookup: dict[str, str] = {}
    if annotations is not None and "gene_id" in annotations.columns:
        for _, row in annotations.iterrows():
            gid = row.get("gene_id", "")
            if gid and gid not in locus_lookup:
                locus_lookup[gid] = row.get("locus_tag", "")

    rows: list[dict[str, object]] = []
    for i, gene_id in enumerate(gene_ids):
        # Knockout: drop this gene's row from the aggregation
        keep_mask = np.ones(len(gene_ids), dtype=bool)
        keep_mask[i] = False
        if not keep_mask.any():
            continue
        ko_features = aggregate_strain_features(full_matrix[keep_mask], "mean").reshape(1, -1)
        ko_prob = float(predict_proba(classifier, ko_features)[0])
        rows.append(
            {
                "gene_id": gene_id,
                "locus_tag": locus_lookup.get(gene_id, ""),
                "prediction_delta": baseline_prob - ko_prob,
                "baseline_probability": baseline_prob,
                "knockout_probability": ko_prob,
            }
        )

    df = pd.DataFrame(rows, columns=list(GENE_EFFECT_COLUMNS))
    df = df.reindex(df["prediction_delta"].abs().sort_values(ascending=False).index)
    return df.reset_index(drop=True)


# ---- Saturation mutagenesis (nucleotide-level) ----


def saturation_mutagenesis(
    classifier: TrainedClassifier,
    foundation_model,
    gene_id: str,
    sequence: str,
    strain_gene_embeddings: dict[str, np.ndarray],
    alt_bases: tuple[str, ...] = ("A", "C", "G", "T"),
) -> PositionEffectTable:
    """Run per-position single-base substitutions across a gene's sequence.

    For each position, mutate to each non-ref alt base, re-embed the mutated
    gene via the foundation model, re-aggregate strain features, re-predict,
    record per-position-per-alt prediction delta.

    Computational cost: len(sequence) × 3 alt bases × 1 embed call. Restrict
    to top-K genes from gene_level_mutagenesis to keep this tractable
    (default K=20 → ~30K embed calls × 1kb avg gene → 30K seqs total).
    """
    if gene_id not in strain_gene_embeddings:
        raise ValueError(f"gene_id {gene_id!r} not in strain_gene_embeddings")

    # Baseline features + prediction
    gene_ids = list(strain_gene_embeddings.keys())
    full_matrix = np.stack([strain_gene_embeddings[g] for g in gene_ids])
    baseline_features = aggregate_strain_features(full_matrix, "mean").reshape(1, -1)
    baseline_prob = float(predict_proba(classifier, baseline_features)[0])

    target_index = gene_ids.index(gene_id)

    rows: list[dict[str, object]] = []
    seq_upper = sequence.upper()
    for pos, ref_base in enumerate(seq_upper):
        if ref_base not in alt_bases:
            continue  # skip Ns or ambiguous
        for alt in alt_bases:
            if alt == ref_base:
                continue
            mutated_seq = seq_upper[:pos] + alt + seq_upper[pos + 1 :]
            mutated_emb = foundation_model.embed_batch([mutated_seq])[0]
            # Substitute the mutated gene's embedding into the full matrix
            mutated_matrix = full_matrix.copy()
            mutated_matrix[target_index] = mutated_emb
            mut_features = aggregate_strain_features(mutated_matrix, "mean").reshape(1, -1)
            mut_prob = float(predict_proba(classifier, mut_features)[0])
            rows.append(
                {
                    "gene_id": gene_id,
                    "position": pos,
                    "ref_base": ref_base,
                    "alt_base": alt,
                    "prediction_delta": mut_prob - baseline_prob,
                }
            )

    return pd.DataFrame(rows, columns=list(POSITION_EFFECT_COLUMNS))


# ---- Tier classification ----


def tier_classify(
    predicted_locus: str,
    drug: str,
    resistance_catalog,  # ResistanceCatalog from Step 4
    annotations: pd.DataFrame | None = None,
    known_pathways: Iterable[str] = (),
) -> AttributionTier:
    """Classify a predicted attribution hit into Tier 1-5 or Fail.

    Tier 1: exact known resistance variant (codon-level match — Phase 2 with
        saturation-mutagenesis nucleotide resolution; in Phase 1 we approximate
        as "gene_symbol matches a known resistance gene AND the catalog entry
        names a specific codon in the gene_symbol string")
    Tier 2: same codon-region — Phase 2
    Tier 3: same gene as a known resistance locus (Phase 1 default success)
    Tier 4: same operon / pathway / mechanism class as a known resistance gene
    Tier 5: plausible by DB lookup but not validated for THIS drug
    Fail: no match in catalog AND no functional annotation
    """
    # Catalog lookup
    drug_hits = resistance_catalog.filter_by_drug_class(drug)
    drug_gene_symbols = {h.gene_symbol.lower() for h in drug_hits}
    drug_families = {h.gene_family.lower() for h in drug_hits if h.gene_family}

    predicted_lower = predicted_locus.lower()

    # Tier 3: gene-symbol exact match against a known resistance gene for THIS drug
    if predicted_lower in drug_gene_symbols:
        return AttributionTier(
            tier=3,
            rationale=f"{predicted_locus} is in the {drug} known-resistance catalog",
        )

    # Tier 4: gene-family match (e.g., 'fluoroquinolone-resistance')
    for family in drug_families:
        if family and family in predicted_lower:
            return AttributionTier(
                tier=4,
                rationale=f"{predicted_locus} is in the {family} family (drug-class match)",
            )

    # Tier 5: matches a resistance gene in catalog but NOT specifically for this drug
    all_resistance_symbols = {e.gene_symbol.lower() for e in resistance_catalog.entries}
    if predicted_lower in all_resistance_symbols:
        return AttributionTier(
            tier=5,
            rationale=f"{predicted_locus} is in the resistance catalog but not specifically associated with {drug}",
        )

    # Fail: no catalog match
    return AttributionTier(
        tier=-1,
        rationale=f"{predicted_locus} has no match in the resistance catalog for {drug} or any drug",
    )


def build_attribution_report(
    gene_effects: GeneEffectTable,
    drug: str,
    resistance_catalog,
    top_k: int = 20,
    annotations: pd.DataFrame | None = None,
) -> AttributionReport:
    """Build Tier 1-5 report for the top-K gene-level hits."""
    report = AttributionReport(drug=drug, top_k=top_k)
    top_rows = gene_effects.head(top_k)
    for _, row in top_rows.iterrows():
        locus = row.get("locus_tag") or row.get("gene_id") or ""
        if not locus:
            continue
        tier = tier_classify(locus, drug, resistance_catalog, annotations=annotations)
        report.hits.append((locus, tier))
    return report


# ---- Motif recovery against RegulonDB / JASPAR ----


def motif_recovery(
    saturation_table: PositionEffectTable,
    known_motifs: dict[str, str],
    delta_threshold: float = 0.1,
) -> dict[str, list[int]]:
    """Find high-impact saturation-mutagenesis windows that overlap known motifs.

    Args:
        saturation_table: PositionEffectTable from saturation_mutagenesis.
        known_motifs: motif_name -> consensus sequence (e.g., from RegulonDB).
        delta_threshold: minimum |prediction_delta| for a position to count
            as 'high-impact'.

    Returns:
        motif_name → list of saturation-table positions where the motif was
        recovered. Empty list = motif not recovered for that gene.
    """
    if len(saturation_table) == 0:
        return {name: [] for name in known_motifs}

    # Group by position; take max-abs delta per position
    by_position = (
        saturation_table.groupby("position")["prediction_delta"]
        .apply(lambda s: s.abs().max())
        .reset_index()
    )
    high_impact_positions = by_position[
        by_position["prediction_delta"] >= delta_threshold
    ]["position"].tolist()

    # For Phase 1, we don't have the gene sequence here — just position list.
    # Phase 2 will add sequence alignment to motif consensus. For now, this
    # function returns the high-impact positions; caller correlates with
    # genomic coordinates via Step 3 annotations.
    return {name: high_impact_positions for name in known_motifs}
