"""Tests for Step 11 — In-silico mutagenesis."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

xgboost = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")

from dna_decode.data.resistance_db import (  # noqa: E402
    ResistanceCatalog,
    ResistanceEntry,
)
from dna_decode.interp.mutagenesis import (  # noqa: E402
    AttributionReport,
    AttributionTier,
    GENE_EFFECT_COLUMNS,
    POSITION_EFFECT_COLUMNS,
    build_attribution_report,
    gene_level_mutagenesis,
    motif_recovery,
    saturation_mutagenesis,
    tier_classify,
)
from dna_decode.models.classifiers import train_xgboost_classifier  # noqa: E402
from dna_decode.models.foundation import MockFoundationModel, ModelMetadata  # noqa: E402


# ---- Test fixtures ----


def _train_synthetic_classifier(seed: int = 0, embedding_dim: int = 8):
    """Train an XGBoost classifier on a synthetic feature signal in dim 0."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((40, embedding_dim))
    y = (X[:, 0] > 0).astype(int)
    return train_xgboost_classifier(X, y, drug_name="cipro", calibrate=False)


def _strain_gene_embeddings(
    seed: int = 1, n_genes: int = 5, embedding_dim: int = 8, signal_gene_index: int = 0
) -> dict[str, np.ndarray]:
    """Build a strain's per-gene embeddings; one gene carries the signal."""
    rng = np.random.default_rng(seed)
    embeddings = rng.standard_normal((n_genes, embedding_dim)).astype(np.float32)
    # Pump signal into the chosen gene's embedding (high value on dim 0)
    embeddings[signal_gene_index, 0] = 5.0
    return {f"gene_{i:02d}": embeddings[i] for i in range(n_genes)}


# ---- Gene-level mutagenesis ----


def test_gene_level_mutagenesis_returns_stable_columns():
    classifier = _train_synthetic_classifier()
    gene_embeds = _strain_gene_embeddings()
    df = gene_level_mutagenesis(classifier, gene_embeds)
    assert list(df.columns) == list(GENE_EFFECT_COLUMNS)


def test_gene_level_mutagenesis_one_row_per_gene():
    classifier = _train_synthetic_classifier()
    gene_embeds = _strain_gene_embeddings(n_genes=7)
    df = gene_level_mutagenesis(classifier, gene_embeds)
    assert len(df) == 7


def test_gene_level_mutagenesis_signal_gene_has_largest_delta():
    """The gene carrying the signal (dim 0 > 0) should rank top by |delta|."""
    classifier = _train_synthetic_classifier()
    gene_embeds = _strain_gene_embeddings(signal_gene_index=2)
    df = gene_level_mutagenesis(classifier, gene_embeds)
    # First row by sorted abs delta should be the signal gene
    assert df.iloc[0]["gene_id"] == "gene_02"


def test_gene_level_mutagenesis_empty_input():
    classifier = _train_synthetic_classifier()
    df = gene_level_mutagenesis(classifier, {})
    assert len(df) == 0
    assert list(df.columns) == list(GENE_EFFECT_COLUMNS)


def test_gene_level_mutagenesis_enriches_with_locus_tag():
    classifier = _train_synthetic_classifier()
    gene_embeds = _strain_gene_embeddings(n_genes=3)
    annotations = pd.DataFrame(
        {
            "gene_id": ["gene_00", "gene_01", "gene_02"],
            "locus_tag": ["TAG_001", "TAG_002", "TAG_003"],
        }
    )
    df = gene_level_mutagenesis(classifier, gene_embeds, annotations)
    locus_for_g0 = df[df["gene_id"] == "gene_00"]["locus_tag"].iloc[0]
    assert locus_for_g0 == "TAG_001"


# ---- Saturation mutagenesis ----


def test_saturation_mutagenesis_returns_stable_columns():
    classifier = _train_synthetic_classifier()
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    gene_embeds = _strain_gene_embeddings(embedding_dim=8)
    df = saturation_mutagenesis(classifier, model, "gene_00", "ATGC", gene_embeds)
    assert list(df.columns) == list(POSITION_EFFECT_COLUMNS)


def test_saturation_mutagenesis_row_count_matches_positions_x_3():
    classifier = _train_synthetic_classifier()
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    gene_embeds = _strain_gene_embeddings(embedding_dim=8)
    df = saturation_mutagenesis(classifier, model, "gene_00", "ATGC", gene_embeds)
    # 4 positions × 3 alt-bases each = 12 rows
    assert len(df) == 12


def test_saturation_mutagenesis_skips_ambiguous_bases():
    classifier = _train_synthetic_classifier()
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    gene_embeds = _strain_gene_embeddings(embedding_dim=8)
    # Sequence with N → that position is skipped
    df = saturation_mutagenesis(classifier, model, "gene_00", "ANGC", gene_embeds)
    # 3 valid positions × 3 alt-bases = 9 rows
    assert len(df) == 9


def test_saturation_mutagenesis_missing_gene_raises():
    classifier = _train_synthetic_classifier()
    model = MockFoundationModel(
        ModelMetadata(name="mock", huggingface_id="x", embedding_dim=8, max_context=100)
    )
    gene_embeds = _strain_gene_embeddings(embedding_dim=8)
    with pytest.raises(ValueError, match="not in"):
        saturation_mutagenesis(classifier, model, "nonexistent", "ATGC", gene_embeds)


# ---- Tier classification ----


def test_tier_classify_tier_3_for_known_drug_gene():
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "fluoroquinolone", "target alt", "CARD", "X"))
    tier = tier_classify("gyrA", "fluoroquinolone", cat)
    assert tier.tier == 3


def test_tier_classify_case_insensitive_match():
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "fluoroquinolone", "target alt", "CARD", "X"))
    tier = tier_classify("GYRA", "fluoroquinolone", cat)
    assert tier.tier == 3


def test_tier_classify_tier_4_for_family_match():
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "fluoroquinolone-resistance", "fluoroquinolone", "alt", "CARD", "X"))
    tier = tier_classify("fluoroquinolone-resistance-gene", "fluoroquinolone", cat)
    assert tier.tier == 4


def test_tier_classify_tier_5_for_other_drug_gene():
    """Gene is in catalog but registered for a different drug."""
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("tetA", "tet-efflux", "tetracycline", "efflux", "CARD", "X"))
    tier = tier_classify("tetA", "fluoroquinolone", cat)
    assert tier.tier == 5


def test_tier_classify_fail_for_unknown_gene():
    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "fluoroquinolone", "alt", "CARD", "X"))
    tier = tier_classify("random_gene_xyz", "fluoroquinolone", cat)
    assert tier.tier == -1


# ---- AttributionReport ----


def test_attribution_report_tier_counts():
    report = AttributionReport(drug="cipro", top_k=5)
    report.hits = [
        ("gyrA", AttributionTier(3, "")),
        ("parC", AttributionTier(3, "")),
        ("random", AttributionTier(-1, "")),
        ("nearby_gene", AttributionTier(4, "")),
        ("unrelated", AttributionTier(5, "")),
    ]
    counts = report.tier_counts()
    assert counts[3] == 2
    assert counts[-1] == 1


def test_attribution_report_fraction_tier_1_to_3():
    report = AttributionReport(drug="cipro", top_k=4)
    report.hits = [
        ("a", AttributionTier(1, "")),
        ("b", AttributionTier(3, "")),
        ("c", AttributionTier(4, "")),
        ("d", AttributionTier(-1, "")),
    ]
    # 2 of 4 in Tier 1-3 → 0.5
    assert report.fraction_tier_1_to_3() == 0.5


def test_attribution_report_fraction_fail():
    report = AttributionReport(drug="cipro", top_k=4)
    report.hits = [
        ("a", AttributionTier(1, "")),
        ("b", AttributionTier(-1, "")),
        ("c", AttributionTier(-1, "")),
        ("d", AttributionTier(3, "")),
    ]
    assert report.fraction_fail() == 0.5


def test_attribution_report_empty_hits():
    report = AttributionReport(drug="cipro", top_k=20)
    assert report.fraction_tier_1_to_3() == 0.0
    assert report.fraction_fail() == 0.0


# ---- build_attribution_report integration ----


def test_build_attribution_report_walks_top_k():
    classifier = _train_synthetic_classifier()
    gene_embeds = _strain_gene_embeddings(n_genes=10, signal_gene_index=3)
    gene_effects = gene_level_mutagenesis(classifier, gene_embeds)

    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gene_03", "f", "fluoroquinolone", "alt", "CARD", "X"))

    report = build_attribution_report(gene_effects, "fluoroquinolone", cat, top_k=5)
    assert report.top_k == 5
    assert len(report.hits) > 0
    # gene_03 is the signal gene + the catalog entry → should land in top-K + Tier 3
    matched = [(name, tier.tier) for name, tier in report.hits if "gene_03" in name.lower()]
    assert matched and matched[0][1] == 3


# ---- Wave 3.5 C8: build_attribution_report tries both gene_id and locus_tag ----


def test_attribution_report_matches_via_locus_tag_when_gene_id_is_not_symbol():
    """Bakta-style row: gene_id='g3' + locus_tag='gyrA'. Catalog has 'gyrA'.
    Old behavior preferred locus_tag → would match; the fix preserves this."""
    from dna_decode.interp.mutagenesis import _best_tier_across_candidates

    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "fluoroquinolone", "alt", "CARD", "X"))
    chosen, tier = _best_tier_across_candidates(cat, "fluoroquinolone", "g3", "gyrA")
    assert tier.tier == 3
    assert chosen == "gyrA"


def test_attribution_report_matches_via_gene_id_when_annotation_set_it_to_symbol():
    """RefSeq-style row: gene_id='gyrA' + locus_tag='b2231'. Old behavior tried
    locus_tag first → would fail (b2231 not in catalog). C8 fix tries gene_id too."""
    from dna_decode.interp.mutagenesis import _best_tier_across_candidates

    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "f", "fluoroquinolone", "alt", "CARD", "X"))
    chosen, tier = _best_tier_across_candidates(cat, "fluoroquinolone", "gyrA", "b2231")
    assert tier.tier == 3
    assert chosen == "gyrA"


def test_attribution_report_picks_better_tier_when_one_candidate_better():
    """When both candidates match catalog at different tiers, pick the lower tier."""
    from dna_decode.interp.mutagenesis import _best_tier_across_candidates

    cat = ResistanceCatalog()
    cat.add(ResistanceEntry("gyrA", "fluoroquinolone-resistance", "fluoroquinolone", "alt", "CARD", "X"))
    # 'gyrA' → Tier 3 (direct match); 'fluoroquinolone-resistance-gene' → Tier 4 (family)
    chosen, tier = _best_tier_across_candidates(
        cat, "fluoroquinolone", "gyrA", "fluoroquinolone-resistance-gene"
    )
    assert tier.tier == 3
    assert chosen == "gyrA"


def test_attribution_report_no_candidates_returns_fail():
    """Empty candidates → Fail tier with informative rationale."""
    from dna_decode.interp.mutagenesis import _best_tier_across_candidates

    cat = ResistanceCatalog()
    chosen, tier = _best_tier_across_candidates(cat, "cipro", "", "")
    assert tier.tier == -1
    assert chosen == ""


# ---- Wave 3.5 M4: motif_recovery emits UserWarning ----


def test_motif_recovery_emits_user_warning_on_call():
    """motif_recovery is a placeholder; first call must warn so silent reliance is loud."""
    # Reset the module-level flag to test fresh-call behavior
    import dna_decode.interp.mutagenesis as mod
    mod._MOTIF_RECOVERY_WARNED = False

    empty = pd.DataFrame(columns=list(POSITION_EFFECT_COLUMNS))
    with pytest.warns(UserWarning, match="placeholder"):
        motif_recovery(empty, {"motifA": "ACGT"})


# ---- motif_recovery ----


def test_motif_recovery_returns_dict_keyed_by_motif_name():
    saturation_df = pd.DataFrame(
        {
            "gene_id": ["g1"] * 6,
            "position": [10, 10, 50, 50, 100, 100],
            "ref_base": ["A"] * 6,
            "alt_base": ["C", "G", "A", "T", "A", "C"],
            "prediction_delta": [0.5, 0.6, 0.01, 0.02, 0.4, 0.3],
        }
    )
    motifs = {"motifA": "ACGT", "motifB": "TTGCAA"}
    out = motif_recovery(saturation_df, motifs, delta_threshold=0.1)
    assert set(out.keys()) == {"motifA", "motifB"}
    # Positions with abs delta >= 0.1: 10, 100; 50 has 0.02 (filtered out)
    assert sorted(out["motifA"]) == [10, 100]


def test_motif_recovery_empty_table():
    empty = pd.DataFrame(columns=list(POSITION_EFFECT_COLUMNS))
    out = motif_recovery(empty, {"motifA": "ACGT"})
    assert out == {"motifA": []}
