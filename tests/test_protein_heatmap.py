"""Tests for the DMS-style protein damage heatmap (pure, offline, CI-safe)."""
from __future__ import annotations

import pytest

from dna_decode.protein_effect.predictor import AA
from dna_decode.viz import protein_heatmap as H


def _logp(seq):
    # tiny synthetic masked-marginals: wt favoured everywhere; one damaging + one favoured off-diagonal
    lp = {}
    for i, wt in enumerate(seq, start=1):
        col = {a: -5.0 for a in AA}
        col[wt] = -0.2                 # wt favoured
        lp[i] = col
    return lp


def test_build_matrix_wt_is_zero_and_sign():
    seq = "MKA"
    lp = _logp(seq)
    lp[2]["N"] = -6.0   # N strongly disfavoured at pos2 (K) -> damaging, positive llr
    lp[2]["R"] = -0.1   # R slightly favoured over K -> negative llr
    m = H.build_matrix(seq, lp)
    assert m["cells"][2]["K"] == 0.0                       # wt is 0 by construction
    assert m["cells"][2]["N"] > 0                          # damaging = positive (logP(wt)-logP(sub))
    assert m["cells"][2]["R"] < 0                          # favoured = negative
    assert m["n_scored_positions"] == 3


def test_missing_position_is_none():
    seq = "MKA"
    lp = _logp(seq)
    del lp[3]                                              # pos 3 unscored
    m = H.build_matrix(seq, lp)
    assert all(m["cells"][3][a] is None for a in AA)
    assert m["n_scored_positions"] == 2


def test_color_diverging_and_none():
    assert H._color(None, 1.0) == "#dddddd"
    assert H._color(0.0, 2.0) == "#f5f5f5"                 # mid = white
    dmg = H._color(2.0, 2.0); fav = H._color(-2.0, 2.0)
    assert dmg == "#b02a37" and fav == "#2b6cb0"           # clamp to red / blue anchors


def test_parse_markers():
    seq = "MKAQ"
    pins = H._parse_markers(["K2N", "Q4L", "X9Z"], seq)
    labels = {p["label"]: p for p in pins}
    assert labels["K2N"]["wt_matches_seq"] is True and labels["K2N"]["mut"] == "N"
    assert labels["Q4L"]["wt_matches_seq"] is True
    assert labels["X9Z"]["wt_matches_seq"] is False       # out of range / wt mismatch, still parsed


def test_html_self_contained_and_honest():
    seq = "MKAQ"
    h = H.build_heatmap_html(seq, _logp(seq), title="test", markers=["K2N"])
    assert h.startswith("<!DOCTYPE html>") and h.rstrip().endswith("</html>")
    for bad in ("src=", "<link", "href=http", "https://", "cdn", "@import"):
        assert bad not in h.lower()
    # honesty rails present: molecular-not-phenotype + frozen surface + separate catalog axis
    assert "not a resistance call" in h.lower()
    assert "molecular" in h.lower() and "phenotype" in h.lower()
    assert "frozen decoder surface" in h.lower()
    assert "catalog markers" in h.lower()                  # K2N pin rendered on its own axis


def test_html_xss_escape():
    seq = "MK"
    h = H.build_heatmap_html(seq, _logp(seq), title="<script>x</script>")
    assert "<script>x</script>" not in h and "&lt;script&gt;x" in h


def test_max_positions_truncates():
    seq = "M" * 50
    h = H.build_heatmap_html(seq, _logp(seq), max_positions=10)
    assert "showing first 10 of 50 positions" in h
