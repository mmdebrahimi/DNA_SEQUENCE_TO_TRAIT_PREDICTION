"""Step 6 tests — G1/G2 gate + cross-genome spike verdict."""
from __future__ import annotations

from dna_decode.genome_map import (
    TIER_CURATED_FUNCTION,
    TIER_DETERMINANT_PHENOTYPE,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.gate import (
    VERDICT_GO,
    VERDICT_NO_GO,
    aggregate_spike_verdict,
    evaluate_gate,
)


def _feature(i, tier, raw_product="", phenotype=None, reason="", raw_gene_symbol=""):
    return {
        "feature_index": i,
        "primary_tier": tier,
        "raw_product": raw_product,
        "raw_gene_symbol": raw_gene_symbol,
        "classification_reason": reason,
        "phenotype": phenotype or [],
    }


def _map(features, *, all_symbol_fallback=False, unknown_rate=0.1, acc="GCA_X"):
    return {
        "genome_accession": acc,
        "features": features,
        "metrics": {
            "all_joins_symbol_fallback": all_symbol_fallback,
            "unknown_under_bakta_db_light": unknown_rate,
        },
    }


def test_g1_pass_three_demotes_one_surface():
    feats = [
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative kinase", reason="low-confidence wording ('putative')"),
        _feature(1, TIER_HOMOLOGY_HYPOTHESIS, "probable transporter"),
        _feature(2, TIER_DETERMINANT_PHENOTYPE, "DNA gyrase subunit A",
                 phenotype=[{"drug": "ciprofloxacin", "phenotype": "R"}]),
        _feature(3, TIER_UNKNOWN, "hypothetical protein"),
    ]
    res = evaluate_gate(_map(feats))
    # 2 demotes + 1 surface = 3 G1 features
    assert len(res["g1_features"]) == 3
    assert res["g1_demote_count"] == 2
    assert res["g1_surface_count"] == 1
    assert res["g1_pass"] is True
    assert res["verdict"] == VERDICT_GO


def test_unknown_rate_alone_is_not_g1():
    # only unknown features -> 0 G1 features -> NO-GO (busywork guard)
    feats = [_feature(i, TIER_UNKNOWN, "hypothetical protein") for i in range(5)]
    res = evaluate_gate(_map(feats))
    assert len(res["g1_features"]) == 0
    assert res["g1_pass"] is False
    assert res["verdict"] == VERDICT_NO_GO


def test_g1_needs_three():
    feats = [
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative kinase"),
        _feature(1, TIER_DETERMINANT_PHENOTYPE, "gyrA", phenotype=[{"drug": "cipro"}]),
    ]
    res = evaluate_gate(_map(feats))
    assert len(res["g1_features"]) == 2
    assert res["g1_pass"] is False


def test_all_symbol_fallback_forces_no_go():
    # even with strong G1, an all-symbol-fallback genome is NO-GO
    feats = [
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative a"),
        _feature(1, TIER_HOMOLOGY_HYPOTHESIS, "putative b"),
        _feature(2, TIER_HOMOLOGY_HYPOTHESIS, "putative c"),
    ]
    res = evaluate_gate(_map(feats, all_symbol_fallback=True))
    assert res["g1_pass"] is True
    assert res["all_joins_symbol_fallback"] is True
    assert res["verdict"] == VERDICT_NO_GO


def test_g2_violation_forces_no_go():
    # a curated-function feature carrying a phenotype = wall breach
    feats = [
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative a"),
        _feature(1, TIER_HOMOLOGY_HYPOTHESIS, "putative b"),
        _feature(2, TIER_HOMOLOGY_HYPOTHESIS, "putative c"),
        _feature(3, TIER_CURATED_FUNCTION, "kinase", phenotype=[{"drug": "x"}]),  # breach
    ]
    res = evaluate_gate(_map(feats))
    assert res["g2_spotcheck"]["pass"] is False
    assert res["verdict"] == VERDICT_NO_GO


def test_demote_only_homology_heavy_genome_passes_g1():
    # the homology-heavy (Gemmata-like) genome: many demotes, no determinants -> G1 via (a)
    feats = [_feature(i, TIER_HOMOLOGY_HYPOTHESIS, "putative protein") for i in range(5)]
    res = evaluate_gate(_map(feats))
    assert res["g1_demote_count"] == 5
    assert res["g1_pass"] is True
    assert res["verdict"] == VERDICT_GO


# ---- aggregate spike verdict ----


def test_aggregate_go_when_one_genome_passes_g1():
    g_pass = evaluate_gate(_map([
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative a"),
        _feature(1, TIER_HOMOLOGY_HYPOTHESIS, "putative b"),
        _feature(2, TIER_HOMOLOGY_HYPOTHESIS, "putative c"),
    ], acc="G1"))
    g_weak = evaluate_gate(_map([_feature(0, TIER_UNKNOWN, "hypothetical protein")], acc="G2"))
    agg = aggregate_spike_verdict([g_pass, g_weak])
    assert agg["spike_g1_pass"] is True
    assert agg["verdict"] == VERDICT_GO


def test_aggregate_no_go_when_any_all_symbol_fallback():
    g_pass = evaluate_gate(_map([
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative a"),
        _feature(1, TIER_HOMOLOGY_HYPOTHESIS, "putative b"),
        _feature(2, TIER_HOMOLOGY_HYPOTHESIS, "putative c"),
    ], acc="G1"))
    g_trap = evaluate_gate(_map([
        _feature(0, TIER_HOMOLOGY_HYPOTHESIS, "putative a"),
        _feature(1, TIER_HOMOLOGY_HYPOTHESIS, "putative b"),
        _feature(2, TIER_HOMOLOGY_HYPOTHESIS, "putative c"),
    ], all_symbol_fallback=True, acc="G2"))
    agg = aggregate_spike_verdict([g_pass, g_trap])
    assert agg["any_all_symbol_fallback"] is True
    assert agg["verdict"] == VERDICT_NO_GO


def test_aggregate_no_go_when_no_genome_passes_g1():
    g1 = evaluate_gate(_map([_feature(0, TIER_UNKNOWN, "hypothetical protein")], acc="G1"))
    g2 = evaluate_gate(_map([_feature(0, TIER_UNKNOWN, "hypothetical protein")], acc="G2"))
    agg = aggregate_spike_verdict([g1, g2])
    assert agg["verdict"] == VERDICT_NO_GO
    assert any("G1 fail" in r for r in agg["reasons"])
