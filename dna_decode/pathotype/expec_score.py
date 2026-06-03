"""Per-gene ExPEC support scoring (EP-4 v0.1 ExPEC-recall hardening, 2026-06-03).

Pure function over a per-gene coverage dict (gene_prefix -> coverage float). Counts distinct ExPEC
support GENES (members of SIDEROPHORES + CAPSULE_SERUM) at or above the confident-coverage bar, so a
multi-gene extraintestinal burden becomes visible to the resolver where the coarse cluster boolean
(present iff ANY member >=0.80) collapses it to one bit.

No I/O, no detection — the per-gene coverage is produced by detection / the committed
data/pathotype_pergene_cache/ artifact (scripts/build_pergene_support_cache.py). Exhaustively
unit-testable on synthetic coverage dicts.
"""
from __future__ import annotations

from dna_decode.pathotype.markers import (
    EXPEC_SUPPORT_GENE_PREFIXES, EXPEC_SUPPORT_GENE_K, CLUSTER_MARKERS,
)

CONFIDENT_COV = 0.80  # same presence bar as detect.CONFIDENT_COV (per-gene count refines the cluster bool)

# Cross-axis support rule (EP-4 v0.1 ExPEC-recall, user round-2 commitment 2026-06-03):
# an ExPEC rescue requires BOTH extraintestinal axes — >=1 iron-acquisition gene AND >=1
# capsule/serum gene (each >=CONFIDENT_COV). Structurally excludes capsule-only genomes (e.g. a lone
# serum-resistance traT hit) that a flat K=1 count over-rescues. No hand-tuned pooled K. The two axes
# ARE the SIDEROPHORES + CAPSULE_SERUM clusters.
IRON_GENES = tuple(CLUSTER_MARKERS["SIDEROPHORES"])
CAPSULE_GENES = tuple(CLUSTER_MARKERS["CAPSULE_SERUM"])


def support_gene_count(pergene_cov: dict[str, float], *, threshold: float = CONFIDENT_COV) -> int:
    """Number of distinct ExPEC support genes >= threshold coverage in `pergene_cov`.

    Keys of `pergene_cov` are lowercase gene-name prefixes (e.g. 'iuta', 'chua', 'kpsmii', 'trat'),
    matching EXPEC_SUPPORT_GENE_PREFIXES. A missing key counts as absent (0.0).
    """
    return sum(1 for g in EXPEC_SUPPORT_GENE_PREFIXES if pergene_cov.get(g, 0.0) >= threshold)


def iron_capsule_counts(pergene_cov: dict[str, float], *, threshold: float = CONFIDENT_COV) -> tuple[int, int]:
    """(#iron genes, #capsule/serum genes) >= threshold coverage."""
    iron = sum(1 for g in IRON_GENES if pergene_cov.get(g, 0.0) >= threshold)
    caps = sum(1 for g in CAPSULE_GENES if pergene_cov.get(g, 0.0) >= threshold)
    return iron, caps


def meets_cross_axis_support(pergene_cov: dict[str, float], *, threshold: float = CONFIDENT_COV) -> bool:
    """True iff BOTH axes present: >=1 iron gene AND >=1 capsule/serum gene (each >= threshold).

    The committed round-2 rule. Structurally excludes capsule-only / iron-only genomes. Replaces the
    earlier flat K=1 burden (which over-rescued a lone-traT genome). Rescue calls stay LOW_CONFIDENCE.
    """
    iron, caps = iron_capsule_counts(pergene_cov, threshold=threshold)
    return iron >= 1 and caps >= 1


def meets_support_burden(pergene_cov: dict[str, float], *, k: int = EXPEC_SUPPORT_GENE_K,
                         threshold: float = CONFIDENT_COV) -> bool:
    """DEPRECATED flat-K burden (over-rescues capsule-only). Retained only so the regression test can
    pin the K=1-vs-cross-axis divergence. Live rescue uses meets_cross_axis_support."""
    return support_gene_count(pergene_cov, threshold=threshold) >= k
