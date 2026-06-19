"""Step 5 tests — map assembler (tiers, phenotype wall, unknown rate, raw fields)."""
from __future__ import annotations

import pandas as pd

from dna_decode.genome_map import (
    TIER_CURATED_FUNCTION,
    TIER_DETERMINANT_PHENOTYPE,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.build_map import build_feature_table, build_genome_map
from dna_decode.genome_map.phenotype_overlay import DeterminantHit, JoinedHit


def _features() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"seqid": "contig_1", "start": 100, "end": 900, "strand": "+", "type": "CDS",
             "gene_id": "G1", "gene_symbol": "gyrA", "locus_tag": "T1", "product": "DNA gyrase subunit A"},
            {"seqid": "contig_1", "start": 2000, "end": 3000, "strand": "-", "type": "CDS",
             "gene_id": "G2", "gene_symbol": "", "locus_tag": "T2", "product": "blaCTX-M-15"},
            {"seqid": "contig_1", "start": 4000, "end": 4500, "strand": "+", "type": "CDS",
             "gene_id": "G3", "gene_symbol": "", "locus_tag": "T3", "product": "putative transporter"},
            {"seqid": "contig_1", "start": 5000, "end": 5400, "strand": "+", "type": "CDS",
             "gene_id": "G4", "gene_symbol": "", "locus_tag": "T4", "product": "hypothetical protein"},
        ]
    )


def _high_join(feature_index, symbol, cls, sub=""):
    return JoinedHit(
        DeterminantHit(symbol, "", cls, sub, "POINTX", f"P{feature_index}", "contig_1", 0, 0),
        feature_index=feature_index, join_confidence="coord",
    )


def _fallback_join(feature_index, symbol, cls):
    return JoinedHit(
        DeterminantHit(symbol, "", cls, "", "POINTX", None, None, None, None),
        feature_index=feature_index, join_confidence="symbol_fallback",
    )


def test_primary_tiers_assigned():
    joined = [_high_join(0, "gyrA_S83L", "QUINOLONE")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ciprofloxacin": {"prediction": "R", "rule": "r"}},
                         drugs=["ciprofloxacin"])
    tiers = [f["primary_tier"] for f in m["features"]]
    assert tiers[0] == TIER_DETERMINANT_PHENOTYPE   # high-confidence determinant beats function
    assert tiers[1] == TIER_CURATED_FUNCTION         # blaCTX-M-15 specific product
    assert tiers[2] == TIER_HOMOLOGY_HYPOTHESIS       # putative
    assert tiers[3] == TIER_UNKNOWN                   # hypothetical


def test_phenotype_wall_only_on_determinant_tier():
    joined = [_high_join(0, "gyrA_S83L", "QUINOLONE")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ciprofloxacin": {"prediction": "R", "rule": "r"}},
                         drugs=["ciprofloxacin"])
    for f in m["features"]:
        if f["primary_tier"] == TIER_DETERMINANT_PHENOTYPE:
            assert f["phenotype"], "determinant feature must carry phenotype"
        else:
            assert f["phenotype"] == [], f"non-determinant feature {f['feature_index']} has phenotype"


def test_symbol_fallback_not_phenotype_but_visible():
    # a symbol-fallback determinant on feature 1 must NOT promote it to determinant-phenotype
    joined = [_fallback_join(1, "blaCTX-M-15", "CEPHALOSPORIN")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 0, "n_symbol_fallback": 1, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ceftriaxone": {"prediction": "R", "rule": "r"}},
                         drugs=["ceftriaxone"])
    f1 = m["features"][1]
    assert f1["primary_tier"] != TIER_DETERMINANT_PHENOTYPE
    assert f1["phenotype"] == []
    # but it IS visible as secondary evidence
    assert any(s["type"] == "determinant_symbol_fallback" for s in f1["secondary_evidence"])
    assert m["metrics"]["all_joins_symbol_fallback"] is True


def test_unknown_rate_db_labelled():
    joined = []
    counts = {"n_main_rows": 0, "n_high_confidence_join": 0, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts)
    # 1 of 4 features is hypothetical -> unknown
    assert m["metrics"]["unknown_under_bakta_db_light"] == 0.25
    assert "unknown_under_bakta_db_light" in m["metrics"]


def test_raw_fields_and_reason_retained():
    joined = []
    counts = {"n_main_rows": 0, "n_high_confidence_join": 0, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts)
    f2 = m["features"][2]
    assert f2["raw_product"] == "putative transporter"
    assert "putative" in f2["classification_reason"]
    assert f2["raw_locus_tag"] == "T3"
    assert f2["raw_feature_type"] == "CDS"


def test_abstain_propagates_to_feature():
    joined = [_high_join(1, "blaX", "CEPHALOSPORIN", "CEPHALOSPORIN")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Klebsiella_pneumoniae", _features(), joined, counts,
                         drug_verdicts={"ceftriaxone": {"prediction": "ABSTAIN", "rule": "cal"}},
                         drugs=["ceftriaxone"])
    f1 = m["features"][1]
    assert f1["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
    assert any(p["phenotype"] == "ABSTAIN" for p in f1["phenotype"])


def test_genome_level_calls_separate():
    joined = [_high_join(0, "gyrA_S83L", "QUINOLONE")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ciprofloxacin": {"prediction": "R", "confidence": "HIGH",
                                                           "n_determinants": 2, "rule": "r"}},
                         drugs=["ciprofloxacin"])
    assert m["metrics"]["genome_level_calls"]["ciprofloxacin"]["prediction"] == "R"


def test_offline_degraded_flag():
    joined = []
    counts = {"n_main_rows": 0, "n_high_confidence_join": 0, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", None, _features(), joined, counts, degraded=True)
    assert m["degraded_coverage"] is True
    # still emits a full map
    assert m["metrics"]["total_features"] == 4


def _acquired_join(feature_index, symbol, cls, sub):
    # an acquired gene (Method EXACTX), NOT a POINT mutation
    return JoinedHit(
        DeterminantHit(symbol, "", cls, sub, "EXACTX", f"P{feature_index}", "contig_1", 0, 0),
        feature_index=feature_index, join_confidence="coord",
    )


def test_refined_excludes_narrow_betalactam_from_ceftriaxone():
    # blaTEM-1 (Subclass BETA-LACTAM, ampicillin-only) must NOT be tagged ceftriaxone
    # (broad class match would have; the refined CEPHALOSPORIN/CARBAPENEM subclass rule excludes it).
    joined = [_acquired_join(1, "blaTEM-1", "BETA-LACTAM", "BETA-LACTAM")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ceftriaxone": {"prediction": "S", "rule": "r"}},
                         drugs=["ceftriaxone"])
    f1 = m["features"][1]
    assert f1["primary_tier"] == TIER_DETERMINANT_PHENOTYPE  # it IS a curated determinant
    # but it carries NO drug — surfaced as DETERMINANT_PRESENT, not ceftriaxone
    assert all(p.get("drug") is None for p in f1["phenotype"])
    assert any(p.get("phenotype") == "DETERMINANT_PRESENT" for p in f1["phenotype"])


def test_refined_keeps_esbl_ceftriaxone():
    # blaCTX-M-15 (Subclass CEPHALOSPORIN) is a real ESBL -> still counts toward ceftriaxone
    joined = [_acquired_join(1, "blaCTX-M-15", "CEPHALOSPORIN", "CEPHALOSPORIN")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ceftriaxone": {"prediction": "R", "rule": "r"}},
                         drugs=["ceftriaxone"])
    pheno = m["features"][1]["phenotype"]
    cef = [p for p in pheno if p.get("drug") == "ceftriaxone"]
    assert cef and cef[0]["drug_rule_counted"] is True
    assert cef[0]["genome_prediction"] == "R"


def test_refined_excludes_efflux_from_ciprofloxacin():
    # a QUINOLONE-class efflux/acquired gene that is NOT a QRDR POINT must NOT be tagged ciprofloxacin
    joined = [_acquired_join(0, "qnrB19", "QUINOLONE", "QUINOLONE")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ciprofloxacin": {"prediction": "S", "rule": "r"}},
                         drugs=["ciprofloxacin"])
    f0 = m["features"][0]
    assert all(p.get("drug") is None for p in f0["phenotype"])  # not ciprofloxacin


def test_refined_keeps_qrdr_point_ciprofloxacin():
    joined = [_high_join(0, "gyrA_S83L", "QUINOLONE")]  # _high_join uses Method POINTX
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ciprofloxacin": {"prediction": "R", "rule": "r"}},
                         drugs=["ciprofloxacin"])
    cip = [p for p in m["features"][0]["phenotype"] if p.get("drug") == "ciprofloxacin"]
    assert cip and cip[0]["drug_rule_counted"] is True


def test_feature_table_flat():
    joined = [_high_join(0, "gyrA_S83L", "QUINOLONE")]
    counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    m = build_genome_map("GCA_1", "Escherichia", _features(), joined, counts,
                         drug_verdicts={"ciprofloxacin": {"prediction": "R", "rule": "r"}},
                         drugs=["ciprofloxacin"])
    table = build_feature_table(m)
    assert len(table) == 4
    assert table[0]["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
    assert table[0]["phenotype"]  # cipro present
    assert table[3]["phenotype"] == ""  # hypothetical -> empty
