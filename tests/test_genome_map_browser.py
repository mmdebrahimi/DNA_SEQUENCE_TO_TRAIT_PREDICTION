"""Tests for the genome-map graphical browser (pure HTML builder; no I/O, no network)."""
from __future__ import annotations

from dna_decode.genome_map.browser import build_genome_map_html, contig_lengths


def _feat(seqid, start, end, tier, ftype="CDS", **kw):
    return {"seqid": seqid, "start": start, "end": end, "strand": "+",
            "primary_tier": tier, "raw_feature_type": ftype, "raw_gene_symbol": kw.get("gene", ""),
            "raw_locus_tag": kw.get("locus", ""), "raw_product": kw.get("product", ""),
            "phenotype": kw.get("phenotype", [])}


def _map(features, **kw):
    m = {"metrics": {"total_features": len(features), "per_tier_counts": kw.get("per_tier", {}),
                     "determinant_phenotype_feature_count": kw.get("dp_count", 0),
                     "genome_level_calls": kw.get("calls", {}), "join_quality": kw.get("jq", {}),
                     "unknown_under_bakta_db_light": kw.get("unk", 0.05),
                     "all_joins_symbol_fallback": kw.get("all_fb", False)}}
    return {"genome_accession": kw.get("acc", "GCA_TEST.1"), "amrfinder_organism": "Escherichia",
            "overlay_status": kw.get("overlay", "FULL"), "features": features, **{k: v for k, v in kw.items()
                                                                                  if k == "virulence_status"},
            "metrics": {**m["metrics"], **kw.get("extra_metrics", {})}}


# ---------------------------------------------------------------- contig lengths

def test_contig_length_from_region_backbone():
    feats = [_feat("c1", 1, 5000, "unknown", ftype="region"),
             _feat("c1", 100, 400, "curated-molecular-function")]
    assert contig_lengths(feats) == {"c1": 5000}


def test_contig_length_falls_back_to_max_end_without_region():
    feats = [_feat("c1", 100, 400, "unknown"), _feat("c1", 900, 1200, "unknown")]
    assert contig_lengths(feats) == {"c1": 1200}


def test_contig_length_handles_multiple_contigs():
    feats = [_feat("c1", 1, 5000, "unknown", ftype="region"),
             _feat("c2", 1, 800, "unknown", ftype="region")]
    assert contig_lengths(feats) == {"c1": 5000, "c2": 800}


# ---------------------------------------------------------------- structure

def test_html_is_self_contained_and_well_formed():
    gm = _map([_feat("c1", 1, 5000, "unknown", ftype="region"),
               _feat("c1", 100, 400, "curated-molecular-function")])
    out = build_genome_map_html(gm, generated="2026-07-11")
    assert out.startswith("<!DOCTYPE html>") and out.rstrip().endswith("</html>")
    assert "GCA_TEST.1" in out
    assert "<script" in out and "http://" not in out and "https://" not in out   # no external assets
    assert "c1 —" in out   # a contig track


def test_backbone_region_is_not_drawn_as_a_block():
    """The whole-contig region feature must NOT become a full-width block covering the track."""
    gm = _map([_feat("c1", 1, 5000, "unknown", ftype="region"),
               _feat("c1", 2500, 2600, "curated-molecular-function")])
    out = build_genome_map_html(gm)
    # exactly one feature block (the CDS), not two; the region is excluded
    assert out.count('class="feat ') == 1


def test_feature_block_count_matches_non_backbone_features():
    feats = [_feat("c1", 1, 9000, "unknown", ftype="region")]
    feats += [_feat("c1", i * 100, i * 100 + 50, "curated-molecular-function") for i in range(1, 6)]
    gm = _map(feats)
    out = build_genome_map_html(gm)
    assert out.count('class="feat ') == 5


# ---------------------------------------------------------------- the honesty wall carries into the visual

def test_determinant_phenotype_feature_is_labelled_others_are_not():
    gm = _map([
        _feat("c1", 1, 9000, "unknown", ftype="region"),
        _feat("c1", 100, 400, "determinant-phenotype", gene="gyrA",
              phenotype=[{"determinant_symbol": "gyrA_S83L", "drug": "ciprofloxacin",
                          "genome_prediction": "R"}]),
        _feat("c1", 500, 800, "curated-molecular-function", gene="dnaA"),
    ])
    out = build_genome_map_html(gm)
    assert '<span class="flabel">gyrA</span>' in out          # DP labelled
    assert '<span class="flabel">dnaA</span>' not in out      # curated NOT labelled on the track


def test_determinant_feature_labels_with_determinant_symbol_when_no_gene_symbol():
    """A DP CDS often has an empty gene symbol; the curated determinant (parC_E84V) is the useful label."""
    gm = _map([
        _feat("c1", 1, 9000, "unknown", ftype="region"),
        _feat("c1", 100, 400, "determinant-phenotype", gene="", locus="LOCUS_00249",
              phenotype=[{"determinant_symbol": "parC_E84V", "drug": "ciprofloxacin"},
                         {"determinant_symbol": "parC_S80I", "drug": "ciprofloxacin"}]),
    ])
    out = build_genome_map_html(gm)
    assert '<span class="flabel">parC_E84V, parC_S80I</span>' in out
    assert "LOCUS_00249" not in out.split('class="flabel"')[1][:40]   # locus tag not used as the label


def test_phenotype_claim_appears_only_in_the_determinant_tooltip():
    gm = _map([
        _feat("c1", 1, 9000, "unknown", ftype="region"),
        _feat("c1", 100, 400, "determinant-phenotype",
              phenotype=[{"determinant_symbol": "parC_S80I", "drug": "ciprofloxacin",
                          "genome_prediction": "R"}]),
        _feat("c1", 500, 800, "unknown", product="hypothetical protein"),
    ])
    out = build_genome_map_html(gm)
    assert "parC_S80I" in out and "ciprofloxacin" in out
    # the phenotype string sits inside the DP block's title, not on the unknown feature
    dp_idx = out.index("t-determinant-phenotype")
    unk_idx = out.index("hypothetical protein")
    assert out.index("parC_S80I") < unk_idx or dp_idx < out.index("parC_S80I")


def test_html_escapes_hostile_raw_text():
    gm = _map([_feat("c1", 1, 5000, "unknown", ftype="region"),
               _feat("c1", 100, 400, "curated-molecular-function",
                     product="<script>alert(1)</script>")])
    out = build_genome_map_html(gm)
    assert "<script>alert(1)</script>" not in out
    assert "&lt;script&gt;" in out


def test_honesty_banner_surfaces_the_wall_and_caveats():
    gm = _map([_feat("c1", 1, 5000, "unknown", ftype="region")],
              overlay="OFFLINE_NO_AMRFINDER", all_fb=True, unk=0.68)
    out = build_genome_map_html(gm)
    assert "Honesty wall" in out
    assert "unknown_under_bakta_db_light" in out
    assert "overlay_status = OFFLINE_NO_AMRFINDER" in out
    assert "all_joins_symbol_fallback = True" in out


def test_empty_features_still_produce_valid_html():
    out = build_genome_map_html(_map([]))
    assert out.startswith("<!DOCTYPE html>") and "</html>" in out


def test_gate_result_summary_is_included_when_supplied():
    gm = _map([_feat("c1", 1, 5000, "unknown", ftype="region")])
    gate = {"g1_features": [1, 2, 3], "g1_demote_count": 2, "g1_surface_count": 1,
            "g2_spotcheck": {"pass": True, "violations": []}, "all_joins_symbol_fallback": False}
    out = build_genome_map_html(gm, gate_result=gate)
    assert "G1 features: 3" in out and "G2 phenotype-wall pass=True" in out
