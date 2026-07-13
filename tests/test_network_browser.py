"""Tests for the co-occurrence network adapter + browser (pure, offline, CI-safe)."""
from __future__ import annotations

import json

import pytest

from dna_decode.viz import network_adapter as A
from dna_decode.viz import network_browser as B


def _write_cooc(tmp_path):
    doc = {
        "verdict": "PASS_LINKAGE_STRUCTURE", "honest_caveats": ["not deduped for clonality"],
        "per_organism": {"escherichia_coli_shigella": {
            "organism": "escherichia_coli_shigella", "n_genomes": 240, "n_determinants": 12,
            "per_determinant": {
                "gyrA_S83L": {"auc": 0.86, "ci_lo": 0.76, "ci_hi": 0.95, "n_present": 200, "linked": True},
                "parC_S80I": {"auc": 0.80, "ci_lo": 0.70, "ci_hi": 0.90, "n_present": 140, "linked": True},
                "sul1": {"auc": 0.55, "ci_lo": 0.40, "ci_hi": 0.70, "n_present": 90, "linked": False},
            },
            "lift_table": {
                "gyrA_S83L": [{"det": "parC_S80I", "cooc": 140, "lift": 1.5, "p_given_target": 0.88},
                              {"det": "sul1", "cooc": 5, "lift": 1.1, "p_given_target": 0.10}],  # below min_cooc
                "parC_S80I": [{"det": "gyrA_S83L", "cooc": 140, "lift": 1.5, "p_given_target": 0.88},
                              {"det": "sul1", "cooc": 40, "lift": 1.3, "p_given_target": 0.28}],
            },
        }},
    }
    p = tmp_path / "cooc.json"; p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def _write_crossaxis(tmp_path):
    doc = {"per_gene": {
        "gyrA_S83L": {"n_present": 200, "generalizes_beyond_lineage": True, "clade_concentrated": False},
        "sul1": {"n_present": 90, "generalizes_beyond_lineage": False, "clade_concentrated": True},
        # parC_S80I intentionally absent -> untested
    }}
    p = tmp_path / "cx.json"; p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def test_build_graph_nodes_edges(tmp_path):
    g = A.build_graph(_write_cooc(tmp_path), "escherichia_coli_shigella",
                      _write_crossaxis(tmp_path), min_cooc=8)
    ids = {n["id"] for n in g["nodes"]}
    assert ids == {"gyrA_S83L", "parC_S80I", "sul1"}
    # the cooc=5 gyrA-sul1 edge is pruned; the cooc=140 and cooc=40 survive; deduped to 2 undirected edges
    assert g["meta"]["n_edges"] == 2
    pairs = {frozenset((e["source"], e["target"])) for e in g["edges"]}
    assert frozenset(("gyrA_S83L", "parC_S80I")) in pairs
    assert frozenset(("parC_S80I", "sul1")) in pairs
    assert frozenset(("gyrA_S83L", "sul1")) not in pairs  # pruned (cooc 5 < 8)


def test_lineage_status_maps_deconfound(tmp_path):
    g = A.build_graph(_write_cooc(tmp_path), "escherichia_coli_shigella",
                      _write_crossaxis(tmp_path), min_cooc=8)
    st = {n["id"]: n["lineage_status"] for n in g["nodes"]}
    assert st["gyrA_S83L"] == A.GENERALIZES        # survives clade-grouping
    assert st["sul1"] == A.LINEAGE_MEDIATED         # clade-concentrated
    assert st["parC_S80I"] == A.UNTESTED            # no cross-axis entry


def test_edge_lineage_mediated_if_either_endpoint(tmp_path):
    g = A.build_graph(_write_cooc(tmp_path), "escherichia_coli_shigella",
                      _write_crossaxis(tmp_path), min_cooc=8)
    em = {frozenset((e["source"], e["target"])): e["lineage_mediated"] for e in g["edges"]}
    # parC_S80I<->sul1: sul1 is lineage-mediated -> edge dashed
    assert em[frozenset(("parC_S80I", "sul1"))] is True
    # gyrA<->parC: neither is lineage-mediated (gyrA generalizes, parC untested) -> solid
    assert em[frozenset(("gyrA_S83L", "parC_S80I"))] is False


def test_no_crossaxis_all_untested(tmp_path):
    g = A.build_graph(_write_cooc(tmp_path), "escherichia_coli_shigella", None, min_cooc=8)
    assert all(n["lineage_status"] == A.UNTESTED for n in g["nodes"])
    assert all(e["lineage_mediated"] is False for e in g["edges"])


def test_organism_missing_raises(tmp_path):
    with pytest.raises(KeyError):
        A.build_graph(_write_cooc(tmp_path), "nonexistent_organism", None)


def test_html_is_self_contained_and_escaped(tmp_path):
    g = A.build_graph(_write_cooc(tmp_path), "escherichia_coli_shigella",
                      _write_crossaxis(tmp_path), min_cooc=8)
    h = B.build_network_html(g)
    assert h.startswith("<!DOCTYPE html>") and h.rstrip().endswith("</html>")
    # no external resource FETCHES (fully offline). The SVG/XML namespace URI is not a fetch, so we
    # target actual fetch patterns rather than any http string.
    for bad in ("src=", "cdn", "<link", "href=http", "@import", "https://"):
        assert bad not in h.lower(), f"external ref {bad!r} leaked into the HTML"
    # the only http occurrence allowed is the w3.org SVG namespace constant
    assert h.lower().count("http") == h.count("http://www.w3.org/2000/svg")
    # de-confound is rendered (dashed dasharray present for the lineage-mediated edge/border)
    assert "dasharray" in h and "lineage-mediated" in h
    # honesty banner present
    assert "associational" in h.lower() and "NOT causal".lower() in h.lower()
    assert "frozen decoder surface" in h.lower()


def test_html_xss_escape(tmp_path):
    g = {"meta": {"organism": "<script>evil</script>", "verdict": "v", "n_genomes": 1,
                  "cooc_artifact": "a.json", "crossaxis_artifact": None, "min_cooc": 8,
                  "n_nodes": 1, "n_edges": 0, "honest_caveats": []},
         "nodes": [{"id": "gyrA_S83L", "prevalence": 10, "auc": 0.8, "linked": True,
                    "lineage_status": A.GENERALIZES}],
         "edges": []}
    h = B.build_network_html(g)
    assert "<script>evil" not in h and "&lt;script&gt;evil" in h


def test_node_render_attrs_deterministic(tmp_path):
    g = A.build_graph(_write_cooc(tmp_path), "escherichia_coli_shigella",
                      _write_crossaxis(tmp_path), min_cooc=8)
    a = B._prep_nodes(g["nodes"]); b = B._prep_nodes(g["nodes"])
    assert a == b  # pure/deterministic
    biggest = max(a, key=lambda n: n["r"])
    assert biggest["id"] == "gyrA_S83L"  # highest prevalence -> largest radius
