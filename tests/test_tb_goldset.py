"""Step 7 — independent gold-set ingestion + scoring (deliverable b)."""
from __future__ import annotations

import json

from dna_decode.organism_rules import tb_goldset
from scripts import score_tb_independent_goldset as gs


def test_load_goldset_absent_is_empty(tmp_path):
    assert tb_goldset.load_goldset(tmp_path / "nope.json") == []
    assert tb_goldset.goldset_available(tmp_path / "nope.json") is False


def test_load_goldset_parses_manifest(tmp_path):
    p = tmp_path / "gold.json"
    p.write_text(json.dumps([
        {"strain_id": "X1", "masked_vcf": "a.vcf", "regeno_vcf": "a.regeno.vcf", "label": "r"},
        {"strain_id": "X2", "masked_vcf": "b.vcf", "regeno_vcf": None, "label": "S"},
    ]), encoding="utf-8")
    iso = tb_goldset.load_goldset(p)
    assert [i.strain_id for i in iso] == ["X1", "X2"]
    assert iso[0].label == "R" and iso[0].regeno_vcf == "a.regeno.vcf"
    assert iso[1].regeno_vcf is None


def test_score_independent_no_goldset_is_blocked():
    out = gs.score_independent({}, {}, {}, drug="rifampicin")
    assert out["status"] == gs.BLOCKED_NO_GOLDSET


def test_score_independent_labels_separately_from_baseline():
    clusters = {"s1": 0, "s2": 0, "s3": 1, "s4": 1}
    labels = {"s1": "R", "s2": "R", "s3": "S", "s4": "S"}
    preds = {"s1": "R", "s2": "R", "s3": "S", "s4": "S"}
    out = gs.score_independent(preds, labels, clusters, drug="rifampicin")
    assert out["status"] == gs.INDEPENDENT_LABEL
    assert out["status"] != "WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE"
    assert "independence_note" in out
    assert out["lineage_collapsed"]["sens"] == 1.0
