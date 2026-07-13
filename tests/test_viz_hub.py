"""Tests for the unified visual-decoder hub (pure, offline, CI-safe)."""
from __future__ import annotations

import pytest

from dna_decode.viz import hub as H


def test_classify_view_file_by_prefix():
    assert H.classify_view_file("network_cooccurrence_ecoli_2026-07-13.html")["scale"] == "network"
    assert H.classify_view_file("circular_genome_GCA_002180195.1.html")["scale"] == "genome"
    assert H.classify_view_file("genome_map_GCA_x.html")["scale"] == "genome"
    assert H.classify_view_file("heatmap_hiv_rt_2026-07-13.html")["scale"] == "protein"
    assert H.classify_view_file("random_file.html") is None
    assert H.classify_view_file("network_cooccurrence_x.json") is None       # non-html ignored


def test_classify_view_href_is_basename():
    v = H.classify_view_file("heatmap_hiv_rt_2026-07-13.html")
    assert v["href"] == "heatmap_hiv_rt_2026-07-13.html"                     # basename, same-dir link
    assert "hiv rt" in v["title"].lower()


def test_build_hub_groups_by_scale_in_order():
    views = [
        {"scale": "protein", "title": "P", "href": "p.html", "blurb": "b", "honesty": "h"},
        {"scale": "genome", "title": "G", "href": "g.html", "blurb": "b", "honesty": "h"},
        {"scale": "network", "title": "N", "href": "n.html", "blurb": "b", "honesty": "h"},
    ]
    h = H.build_hub_html(views)
    # genome (scale 1) section must appear before network (2) before protein (3)
    assert h.index("which parts do what") < h.index("which parts relate") < h.index("what a specific edit does")
    for href in ("g.html", "n.html", "p.html"):
        assert f'href="{href}"' in h


def test_hub_self_contained_and_honest():
    h = H.build_hub_html([{"scale": "network", "title": "N", "href": "n.html", "blurb": "b", "honesty": "h"}])
    assert h.startswith("<!DOCTYPE html>") and h.rstrip().endswith("</html>")
    for bad in ("src=", "<link", "href=http", "https://", "cdn", "@import"):
        assert bad not in h.lower()
    assert "associational" in h.lower() and "not a resistance call" in h.lower()
    assert "frozen decoder surface" in h.lower()
    assert "semantic zoom" in h.lower()


def test_hub_empty_scale_renders_placeholder():
    h = H.build_hub_html([])            # no views -> all three scale sections show the placeholder
    assert h.lower().count("no view generated") == 3


def test_hub_xss_escape():
    h = H.build_hub_html([{"scale": "network", "title": "<script>x</script>", "href": "n.html",
                           "blurb": "b", "honesty": "h"}])
    assert "<script>x</script>" not in h and "&lt;script&gt;x" in h
