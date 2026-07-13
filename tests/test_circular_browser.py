"""Tests for the circular genome-ring browser (pure, offline, CI-safe)."""
from __future__ import annotations

import math

import pytest

from dna_decode.genome_map import circular_browser as C


def _gm():
    return {
        "genome_accession": "GCA_TEST.1", "amrfinder_organism": "Escherichia",
        "features": [
            {"seqid": "c1", "start": 1, "end": 1000, "raw_feature_type": "region", "primary_tier": "unknown"},
            {"seqid": "c1", "start": 10, "end": 90, "strand": "+", "raw_feature_type": "CDS",
             "primary_tier": "curated-molecular-function", "raw_gene_symbol": "aaa"},
            {"seqid": "c1", "start": 200, "end": 260, "strand": "-", "raw_feature_type": "CDS",
             "primary_tier": "determinant-phenotype", "raw_gene_symbol": "gyrA",
             "phenotype": [{"determinant_symbol": "gyrA_S83L", "drug": "ciprofloxacin",
                            "genome_prediction": "R"}]},
            {"seqid": "c2", "start": 1, "end": 500, "raw_feature_type": "region", "primary_tier": "unknown"},
            {"seqid": "c2", "start": 5, "end": 60, "strand": "+", "raw_feature_type": "CDS",
             "primary_tier": "homology-only-hypothesis"},
        ],
        "metrics": {"determinant_phenotype_feature_count": 1, "per_tier_counts": {}, "total_features": 5,
                    "join_quality": {}},
    }


def test_polar_top_is_frac_zero():
    x, y = C._pt(100, 100, 50, 0.0)
    assert abs(x - 100) < 1e-6 and abs(y - 50) < 1e-6         # frac 0 -> top (12 o'clock)


def test_layout_spans_sum_within_circle():
    gm = _gm()
    lengths = C.contig_lengths(gm["features"])
    spans, total, order = C._layout(gm["features"], lengths)
    assert total == 1500 and order == ["c1", "c2"]            # c1 longer -> first, largest span
    assert spans["c1"][1] > spans["c2"][1]
    # spans + gaps stay within [0,1)
    last = spans[order[-1]]
    assert last[0] + last[1] <= 1.0


def test_feature_frac_maps_into_contig_span():
    gm = _gm()
    lengths = C.contig_lengths(gm["features"])
    spans, _t, _o = C._layout(gm["features"], lengths)
    start = C._feature_frac(spans, "c1", 1)
    mid = C._feature_frac(spans, "c1", 500)
    assert spans["c1"][0] <= start < mid <= spans["c1"][0] + spans["c1"][1]


def test_arc_path_is_valid_svg():
    p = C._arc(100, 100, 40, 50, 0.1, 0.3)
    assert p.startswith("M") and p.endswith("Z") and " A" in p and " L" in p


def test_html_self_contained_and_honest():
    h = C.build_circular_html(_gm())
    assert h.startswith("<!DOCTYPE html>") and h.rstrip().endswith("</html>")
    for bad in ("src=", "<link", "href=http", "https://", "cdn", "@import"):
        assert bad not in h.lower()
    assert "<svg" in h and "determinant-phenotype" in h
    assert "honesty wall" in h.lower()                        # reused linear banner
    assert 'class="dl"' in h and "gyrA" in h                  # determinant callout label rendered (raw_gene_symbol wins)


def test_backbone_not_drawn_as_feature():
    # the whole-contig `region` backbone defines the arc but is not a feature arc
    h = C.build_circular_html(_gm())
    # 3 non-backbone features -> 3 feat arcs; region rows excluded
    assert h.count('class="feat ') == 3


def test_xss_escape():
    gm = _gm(); gm["genome_accession"] = "<script>x</script>"
    h = C.build_circular_html(gm)
    assert "<script>x</script>" not in h and "&lt;script&gt;x" in h
