"""Virulence-determinant overlay tests (Steps 2-4).

Covers the tier constant (Step 2), the overlay adapter + join + full-contract genome
pathotype call + organism scope (Step 3), and the build_map virulence tier + metric
isolation from the AMR spike gate (Step 4).
"""
from __future__ import annotations

import pandas as pd

from dna_decode.genome_map import (
    PHENOTYPE_TIER,
    TIER_DETERMINANT_PHENOTYPE,
    TIER_PRECEDENCE,
    TIER_VIRULENCE_DETERMINANT,
)
from dna_decode.genome_map.build_map import build_genome_map
from dna_decode.genome_map.gate import aggregate_spike_verdict, evaluate_gate
from dna_decode.genome_map.phenotype_overlay import DeterminantHit, JoinedHit
from dna_decode.genome_map.virulence_overlay import (
    _pergene_cov_from_per_gene,
    all_virulence_symbol_fallback,
    cluster_pathotype_context,
    genome_pathotype_call,
    join_virulence,
    parse_virulence_hits,
    virulence_organism_in_scope,
)
from dna_decode.pathotype.markers import CLUSTER_MARKERS
from dna_decode.pathotype.vf_runner import NON_INDEPENDENCE_CAVEAT

DB_SHA = "deadbeefdeadbeef"


# ---------- fixtures ----------

def _features() -> pd.DataFrame:
    return pd.DataFrame([
        {"seqid": "contig_1", "start": 100, "end": 900, "strand": "+", "type": "CDS",
         "gene_id": "G1", "gene_symbol": "", "locus_tag": "T1", "product": "Shiga toxin 2 subunit A"},
        {"seqid": "contig_1", "start": 2000, "end": 2860, "strand": "+", "type": "CDS",
         "gene_id": "G2", "gene_symbol": "", "locus_tag": "T2", "product": "intimin"},
        {"seqid": "contig_1", "start": 4000, "end": 4860, "strand": "+", "type": "CDS",
         "gene_id": "G3", "gene_symbol": "fimH", "locus_tag": "T3", "product": "type 1 fimbrial adhesin"},
        {"seqid": "contig_1", "start": 6000, "end": 6500, "strand": "+", "type": "CDS",
         "gene_id": "G4", "gene_symbol": "", "locus_tag": "T4", "product": "hypothetical protein"},
    ])


def _per_hit(allele_id, cluster, sseqid, start, stop, *, strand="+", cov=100.0, ident=99.0):
    return {"allele_id": allele_id, "vf_gene": allele_id.split(":")[0], "cluster": cluster,
            "sseqid": sseqid, "start": start, "stop": stop, "strand": strand,
            "percent_identity": ident, "percent_coverage": cov, "called": True}


def _per_cluster(called: dict[str, bool]) -> dict:
    return {c: {"called": bool(called.get(c, False)),
                "best_gene": f"{c.lower()}:acc" if called.get(c) else None,
                "percent_identity": 99.0 if called.get(c) else None,
                "percent_coverage": 100.0 if called.get(c) else None}
            for c in CLUSTER_MARKERS}


def _big_contigs(total_bp: int):
    """One contig sized to land in a chosen assembly_qc band."""
    return [("contig_1", "A" * total_bp)]


# ---------- Step 2: tier constant ----------

def test_tier_constant_precedence():
    assert TIER_VIRULENCE_DETERMINANT == "virulence-determinant"
    assert TIER_PRECEDENCE.index(TIER_VIRULENCE_DETERMINANT) == 1  # immediately after AMR
    assert TIER_PRECEDENCE.index(TIER_DETERMINANT_PHENOTYPE) == 0
    # the AMR phenotype wall is unchanged — virulence has its own presence wall
    assert PHENOTYPE_TIER == TIER_DETERMINANT_PHENOTYPE


# ---------- Step 3: parse_virulence_hits ----------

def test_parse_virulence_hits_includes_tandem_and_excludes_uncalled():
    res = {"status": "ok", "db_sha": DB_SHA, "per_hit": [
        _per_hit("stx2A:a", "STX2", "contig_1", 100, 900),
        _per_hit("blaTEM-1:b", None, "contig_1", 2000, 2860),   # unclustered, copy 1
        _per_hit("blaTEM-1:b", None, "contig_1", 4000, 4860),   # unclustered, copy 2 (tandem)
        {**_per_hit("xyz:c", None, "contig_1", 9, 9), "called": False},  # not called -> dropped
    ]}
    hits = parse_virulence_hits(res)
    assert len(hits) == 3  # uncalled excluded; both tandem copies retained
    stx = next(h for h in hits if h.symbol == "stx2A")
    assert stx.cls == "VIRULENCE" and stx.subclass == "STX2"
    assert stx.protein_id is None and stx.contig == "contig_1"
    assert (stx.start, stx.stop) == (100, 900)
    tem = [h for h in hits if h.symbol == "blaTEM-1"]
    assert len(tem) == 2 and {(h.start, h.stop) for h in tem} == {(2000, 2860), (4000, 4860)}
    assert all(h.subclass == "" for h in tem)  # unclustered -> empty subclass


# ---------- Step 3: join_virulence ----------

def test_join_virulence_coord_high_confidence_tandem():
    res = {"status": "ok", "db_sha": DB_SHA, "per_hit": [
        _per_hit("stx2A:a", "STX2", "contig_1", 100, 900),
        _per_hit("eae:b", "LEE", "contig_1", 2000, 2860),
    ]}
    hits = parse_virulence_hits(res)
    joined, counts = join_virulence(_features(), hits, contig_name_map=None,
                                    contig_names=["contig_1"])
    assert counts["n_vf_rows"] == 2
    assert counts["n_high_confidence_join"] == 2
    assert counts["all_virulence_joins_symbol_fallback"] is False
    assert counts["n_ambiguous_contig"] == 0
    assert all(jh.is_high_confidence for jh in joined)


def test_join_virulence_symbol_fallback_excluded_and_flagged():
    # a hit with no coords whose symbol matches a feature gene_symbol -> symbol-fallback only
    hits = parse_virulence_hits({"status": "ok", "per_hit": [
        {**_per_hit("fimH:a", None, None, None, None), "sseqid": None, "start": None, "stop": None},
    ]})
    joined, counts = join_virulence(_features(), hits, contig_names=["contig_1"])
    assert counts["n_vf_rows"] == 1
    assert counts["n_high_confidence_join"] == 0
    assert counts["n_symbol_fallback"] == 1
    assert counts["all_virulence_joins_symbol_fallback"] is True
    assert not any(jh.is_high_confidence for jh in joined)


def test_join_virulence_ambiguous_contig_counted():
    hits = parse_virulence_hits({"status": "ok", "per_hit": [
        _per_hit("stx2A:a", "STX2", "dupcontig", 100, 900),
    ]})
    # the genome has TWO contigs sharing the first-token "dupcontig" -> ambiguous join
    _, counts = join_virulence(_features(), hits, contig_names=["dupcontig", "dupcontig", "x"])
    assert counts["n_ambiguous_contig"] == 1


# ---------- Step 3: genome_pathotype_call (C1 full contract) ----------

def test_genome_pathotype_call_low_qc_never_commensal():
    res = {"status": "ok", "db_sha": DB_SHA,
           "per_cluster": _per_cluster({}), "per_gene": {}}
    call = genome_pathotype_call(res, _big_contigs(1_000))  # tiny -> QC FAIL
    assert call["status"] == "ok"
    assert call["derived_call"]["primary"] == "AMBIGUOUS_LOW_QC"
    assert call["derived_call"]["confidence_tier"] != "CONFIDENT"
    assert call["caveat"] == NON_INDEPENDENCE_CAVEAT and call["db_sha"] == DB_SHA


def test_genome_pathotype_call_stx_lee_is_ehec_with_caveat():
    res = {"status": "ok", "db_sha": DB_SHA,
           "per_cluster": _per_cluster({"STX2": True, "LEE": True}), "per_gene": {}}
    call = genome_pathotype_call(res, _big_contigs(4_500_000))  # PASS band
    dc = call["derived_call"]
    assert dc["primary"] == "EHEC_COMPATIBLE"
    assert call["caveat"] == NON_INDEPENDENCE_CAVEAT
    assert call["db_sha"] == DB_SHA


def test_genome_pathotype_call_cross_axis_expec_rescue():
    # 0 strong adhesins, but >=1 iron gene AND >=1 capsule gene at confident coverage
    res = {"status": "ok", "db_sha": DB_SHA,
           "per_cluster": _per_cluster({}),
           "per_gene": {"iutA:a": {"percent_coverage": 96.0},
                        "kpsMII:b": {"percent_coverage": 95.0}}}
    call = genome_pathotype_call(res, _big_contigs(4_500_000))
    dc = call["derived_call"]
    assert dc["primary"] == "ExPEC_COMPATIBLE"
    assert dc["confidence_tier"] == "LOW_CONFIDENCE"
    assert dc["rule_id"] == "RULE_EXPEC_002"


def test_genome_pathotype_call_insufficient_context():
    # VF unavailable -> insufficient_context (NOT a confident commensal)
    assert genome_pathotype_call({"status": "unavailable", "db_sha": DB_SHA},
                                 _big_contigs(4_500_000))["status"] == "insufficient_context"
    # no contigs -> insufficient_context even with an ok VF result
    ok = {"status": "ok", "db_sha": DB_SHA, "per_cluster": _per_cluster({}), "per_gene": {}}
    assert genome_pathotype_call(ok, [])["status"] == "insufficient_context"


def test_genome_pathotype_call_no_support_is_confident_commensal_not_insufficient():
    """Absent markers + absent support at good QC is a REAL confident COMMENSAL call —
    NEVER insufficient_context (that status is reserved for unavailable VF / no contigs)."""
    res = {"status": "ok", "db_sha": DB_SHA, "per_cluster": _per_cluster({}), "per_gene": {}}
    call = genome_pathotype_call(res, _big_contigs(4_500_000))
    assert call["status"] == "ok"
    assert call["derived_call"]["primary"] == "COMMENSAL_LOW_MARKER_BURDEN"
    assert call["derived_call"]["confidence_tier"] == "CONFIDENT"


# ---------- Step 3: pure helpers ----------

def test_cluster_pathotype_context_maps_unmaps_and_empty():
    assert cluster_pathotype_context("STX2") == ["STEC/EHEC"]
    assert cluster_pathotype_context("LT") == ["ETEC"]
    # the returned list is a copy (mutating it must not corrupt the catalog)
    ctx = cluster_pathotype_context("STX2")
    ctx.append("X")
    assert cluster_pathotype_context("STX2") == ["STEC/EHEC"]
    # unmapped / falsy -> []
    for c in ("NOT_A_CLUSTER", None, ""):
        assert cluster_pathotype_context(c) == []


def test_all_virulence_symbol_fallback_unit():
    assert all_virulence_symbol_fallback({"n_vf_rows": 0, "n_high_confidence_join": 0}) is False
    assert all_virulence_symbol_fallback({"n_vf_rows": 3, "n_high_confidence_join": 2}) is False
    assert all_virulence_symbol_fallback({"n_vf_rows": 3, "n_high_confidence_join": 0}) is True
    assert all_virulence_symbol_fallback({}) is False  # missing keys default to 0


def test_pergene_cov_from_per_gene_max_across_alleles_and_prefix_filter():
    pg = {
        "iutA:a": {"percent_coverage": 80.0},   # support prefix iuta
        "iutA:b": {"percent_coverage": 96.0},   # higher copy of the same prefix -> max wins
        "kpsMII:c": {"percent_coverage": 91.0},  # capsule support prefix
        "nomatch:d": {"percent_coverage": 99.0},  # no support prefix -> excluded
    }
    cov = _pergene_cov_from_per_gene(pg)
    assert cov["iuta"] == 0.96           # coverage normalized to [0,1]; max across the two alleles
    assert cov["kpsmii"] == 0.91
    assert "nomatch" not in cov
    assert _pergene_cov_from_per_gene({}) == {}
    # a missing/None coverage degrades to 0.0 (no crash)
    assert _pergene_cov_from_per_gene({"iutA:x": {}})["iuta"] == 0.0


# ---------- Step 3: organism scope ----------

def test_virulence_organism_in_scope():
    for o in ("Escherichia", "Escherichia_coli", "Escherichia coli",
              "Escherichia_coli_Shigella", "Shigella", "Shigella_flexneri"):
        assert virulence_organism_in_scope(o) is True
    for o in ("Klebsiella_pneumoniae", "Salmonella", "", None, "Acinetobacter"):
        assert virulence_organism_in_scope(o) is False


# ---------- Step 4: build_map virulence tier + the wall ----------

def _amr_high_join(feature_index, symbol, cls):
    return JoinedHit(
        DeterminantHit(symbol, "", cls, cls, "EXACTX", f"P{feature_index}", "contig_1", 0, 0),
        feature_index=feature_index, join_confidence="coord",
    )


def _vir_overlay(per_hit_list, *, contig_names=("contig_1",)):
    res = {"status": "ok", "db_sha": DB_SHA, "per_hit": per_hit_list}
    hits = parse_virulence_hits(res)
    return join_virulence(_features(), hits, contig_names=list(contig_names))


def test_build_map_virulence_tier_and_walls():
    # AMR determinant at feature 1; a high-confidence VF determinant at feature 0;
    # a symbol-fallback VF (fimH, no coords) at feature 2.
    amr = [_amr_high_join(1, "blaCTX-M-15", "CEPHALOSPORIN")]
    amr_counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    vir_joined, vir_counts = _vir_overlay([
        _per_hit("stx2A:a", "STX2", "contig_1", 100, 900),
        {**_per_hit("fimH:b", None, None, None, None), "sseqid": None, "start": None, "stop": None},
    ])
    patho = {"status": "ok", "derived_call": {"primary": "EHEC_COMPATIBLE"}, "db_sha": DB_SHA}
    gm = build_genome_map(
        "GCA_V", "Escherichia", _features(), amr, amr_counts,
        drug_verdicts={"ceftriaxone": {"prediction": "R", "determinants": [{"symbol": "blaCTX-M-15"}], "rule": "r"}},
        drugs=["ceftriaxone"],
        virulence_joined_hits=vir_joined, virulence_join_counts=vir_counts,
        pathotype_call=patho, virulence_db_sha=DB_SHA,
    )
    feats = gm["features"]
    # feature 0 -> virulence-determinant with the virulence wall populated
    assert feats[0]["primary_tier"] == TIER_VIRULENCE_DETERMINANT
    assert feats[0]["virulence"] and feats[0]["virulence"][0]["vf_gene"] == "stx2A"
    assert feats[0]["virulence"][0]["cluster"] == "STX2"
    assert feats[0]["virulence"][0]["pathotype_context"] == ["STEC/EHEC"]
    assert feats[0]["virulence"][0]["db_sha"] == DB_SHA
    assert feats[0]["phenotype"] == []          # virulence feature carries NO AMR phenotype
    # feature 1 -> AMR determinant-phenotype (AMR wins); phenotype set, virulence empty
    assert feats[1]["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
    assert feats[1]["phenotype"] and feats[1]["virulence"] == []
    # feature 2 -> NOT the virulence tier (symbol-fallback), visible as secondary evidence
    assert feats[2]["primary_tier"] != TIER_VIRULENCE_DETERMINANT
    assert feats[2]["virulence"] == []
    assert any(s["type"] == "virulence_symbol_fallback" for s in feats[2]["secondary_evidence"])
    # raw fields retained
    assert feats[0]["raw_product"] == "Shiga toxin 2 subunit A"

    # WALLS: phenotype only on AMR tier; virulence only on virulence tier
    for f in feats:
        if f["phenotype"]:
            assert f["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
        if f["virulence"]:
            assert f["primary_tier"] == TIER_VIRULENCE_DETERMINANT

    m = gm["metrics"]
    assert m["virulence_determinant_feature_count"] == 1
    assert m["all_virulence_joins_symbol_fallback"] is False
    assert m["genome_pathotype_call"] is patho           # surfaced separately
    assert m["virulence_join_quality"]["n_ambiguous_contig"] == 0


def test_m2_virulence_metrics_isolated_from_amr_gate():
    """An all-symbol-fallback VF set + a clean AMR overlay: the AMR all_joins_symbol_fallback
    + the spike gate are IDENTICAL to a map built with NO virulence at all (M2)."""
    amr = [_amr_high_join(1, "blaCTX-M-15", "CEPHALOSPORIN")]
    amr_counts = {"n_main_rows": 1, "n_high_confidence_join": 1, "n_symbol_fallback": 0, "n_unjoined": 0}
    verdicts = {"ceftriaxone": {"prediction": "R", "determinants": [{"symbol": "blaCTX-M-15"}], "rule": "r"}}

    # all-symbol-fallback virulence set (fimH, no coords)
    vir_joined, vir_counts = _vir_overlay([
        {**_per_hit("fimH:b", None, None, None, None), "sseqid": None, "start": None, "stop": None},
    ])
    gm_vir = build_genome_map("GCA_V", "Escherichia", _features(), amr, amr_counts,
                              drug_verdicts=verdicts, drugs=["ceftriaxone"],
                              virulence_joined_hits=vir_joined, virulence_join_counts=vir_counts)
    gm_plain = build_genome_map("GCA_V", "Escherichia", _features(), amr, amr_counts,
                                drug_verdicts=verdicts, drugs=["ceftriaxone"])

    # virulence metric is set, AND it did NOT leak into the AMR key
    assert gm_vir["metrics"]["all_virulence_joins_symbol_fallback"] is True
    assert gm_vir["metrics"]["all_joins_symbol_fallback"] is False
    # the AMR spike gate is byte-identical with vs without virulence
    g_vir, g_plain = evaluate_gate(gm_vir), evaluate_gate(gm_plain)
    for k in ("verdict", "all_joins_symbol_fallback", "overlay_go", "g1_pass"):
        assert g_vir[k] == g_plain[k]
    assert (aggregate_spike_verdict([g_vir])["verdict"]
            == aggregate_spike_verdict([g_plain])["verdict"])
