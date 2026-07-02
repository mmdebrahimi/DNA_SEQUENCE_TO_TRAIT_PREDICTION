"""Lineage/structure-DE-CONFOUNDING primitives + the dataset-candidate scorecard.

The reusable methodology distilled from the 2026-07-02 cross-substrate work (yeast growth + copy-number
attribution; DepMap drug-response multimodal attribution). The session's load-bearing law: an interpretable
decoder recovers the CANONICAL causal gene, de-confounded from population/lineage structure, exactly when the
GENOTYPE FEATURE TYPE matches the MECHANISM TYPE (point mutation / copy number / expression).

These functions are the data-INDEPENDENT core (no dataset download at import). The heavy research runners live
in `scripts/` and import from here. See `wiki/yeast_growth_decoder_result_2026-07-02.md`,
`wiki/yeast_cnv_attribution_result_2026-07-02.md`, `wiki/depmap_multimodal_result_2026-07-02.md`.
"""
from dna_decode.deconfound.deconfound import (  # noqa: F401
    cluster_from_distance,
    cv_r2,
    group_centered_association,
    group_centered_biomarker_t,
    group_centered_spearman,
    permutation_null,
    r2,
    univariate_top,
    within_group_r2,
)
from dna_decode.deconfound.scorecard import (  # noqa: F401
    GATE_KEYS,
    GATES,
    Candidate,
    decoder_type,
    rank,
    score,
)

__all__ = [
    "r2", "cv_r2", "within_group_r2", "cluster_from_distance", "group_centered_spearman",
    "group_centered_biomarker_t", "univariate_top", "permutation_null", "group_centered_association",
    "GATES", "GATE_KEYS", "Candidate", "decoder_type", "score", "rank",
]
