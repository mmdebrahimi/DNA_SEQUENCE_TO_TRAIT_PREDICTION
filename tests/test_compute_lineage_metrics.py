"""Unit tests for the per-cohort lineage-metrics recompute (scripts/compute_lineage_metrics.py).

No network, no Docker: the Mash call lives only in process_cohort/main; every test here
exercises the pure helpers (parse / gate / reconcile / assembly / upsert) against on-disk
fixtures + monkeypatched rule path.
"""
from __future__ import annotations

import json

import pytest

import scripts.compute_lineage_metrics as m
from dna_decode.data.cell_key import canonical_cell_key, cell_key_str


# --------------------------------------------------------------------------- #
# cohort <-> artifact resolution
# --------------------------------------------------------------------------- #
def test_parse_cohort_dir():
    assert m.parse_cohort_dir("klebsiella_provdisjoint_ciprofloxacin") == ("klebsiella", "ciprofloxacin")
    assert m.parse_cohort_dir("escherichia_coli_shigella_provdisjoint_ciprofloxacin") == (
        "escherichia_coli_shigella", "ciprofloxacin")


def test_parse_cohort_dir_rejects_non_provdisjoint():
    with pytest.raises(ValueError):
        m.parse_cohort_dir("klebsiella_indep_ciprofloxacin")


def test_find_artifact_picks_latest(tmp_path):
    wiki = tmp_path / "wiki"
    wiki.mkdir()
    (wiki / "provenance_disjoint_validation_klebsiella_cipro_2026-06-01.json").write_text("{}")
    (wiki / "provenance_disjoint_validation_klebsiella_cipro_2026-06-10.json").write_text("{}")
    got = m.find_artifact("klebsiella", "ciprofloxacin", wiki=wiki)
    assert got is not None and got.name.endswith("2026-06-10.json")


def test_find_artifact_absent(tmp_path):
    assert m.find_artifact("nope", "ciprofloxacin", wiki=tmp_path) is None


def test_read_selected(tmp_path):
    sel = tmp_path / "selected.tsv"
    sel.write_text("GCA_1\tR\nGCA_2\tS\n")
    assert m.read_selected(sel) == {"GCA_1": 1, "GCA_2": 0}


def test_canonical_key_shared():
    assert cell_key_str("Klebsiella", "Ciprofloxacin") == "klebsiella|ciprofloxacin"
    assert canonical_cell_key("  Klebsiella ", "CIPRO") == ("klebsiella", "cipro")


# --------------------------------------------------------------------------- #
# M1 genome-completeness gate
# --------------------------------------------------------------------------- #
def test_fasta_ok(tmp_path):
    good = tmp_path / "g.fna"
    good.write_text(">contig1\nACGT\n")
    assert m._fasta_ok(good)
    empty = tmp_path / "e.fna"
    empty.write_text("")
    assert not m._fasta_ok(empty)
    notfasta = tmp_path / "n.fna"
    notfasta.write_text("ACGT no header\n")
    assert not m._fasta_ok(notfasta)
    assert not m._fasta_ok(tmp_path / "missing.fna")


def test_ensure_cohort_genomes_partial(tmp_path):
    refseq = tmp_path / "refseq"
    (refseq / "GCA_1").mkdir(parents=True)
    (refseq / "GCA_1" / "genome.fna").write_text(">c\nACGT\n")
    # GCA_2 absent; fetch=False so the gate just reports it missing (no network)
    present, missing = m.ensure_cohort_genomes({"GCA_1": 1, "GCA_2": 0}, refseq, fetch=False)
    assert list(present) == ["GCA_1"] and missing == ["GCA_2"]


# --------------------------------------------------------------------------- #
# M4 reconciliation (raises on mismatch)
# --------------------------------------------------------------------------- #
def _patch_rule_path(monkeypatch, preds_by_acc):
    """Stub _run_dir (always 'found') + call_resistance (acc -> canned prediction)."""
    monkeypatch.setattr(m, "_run_dir", lambda acc, own, glob: own / acc)

    def fake_call(main_tsv, drug, organism=None, **kw):
        acc = main_tsv.parent.name
        return {"prediction": preds_by_acc[acc]}

    monkeypatch.setattr(m, "call_resistance", fake_call)


def test_reconcile_passes_on_match(tmp_path, monkeypatch):
    selected = {"a": 1, "b": 1, "c": 0, "d": 0}
    preds = {"a": "R", "b": "S", "c": "S", "d": "R"}  # tp1 fn1 tn1 fp1
    _patch_rule_path(monkeypatch, preds)
    artifact_metrics = {"tp": 1, "fp": 1, "tn": 1, "fn": 1, "sens": 0.5, "spec": 0.5, "n_scored": 4}
    raw, got_preds = m.reconcile_raw_metrics(selected, tmp_path, "glob", "cipro", "Org", artifact_metrics)
    assert raw["tp"] == 1 and raw["sens"] == 0.5 and got_preds == preds


def test_reconcile_raises_on_mismatch(tmp_path, monkeypatch):
    selected = {"a": 1, "b": 0}
    _patch_rule_path(monkeypatch, {"a": "R", "b": "S"})  # recomputes tp1 tn1 sens1.0 spec1.0
    bad_artifact = {"tp": 2, "fp": 0, "tn": 0, "fn": 0, "sens": 1.0, "spec": None, "n_scored": 2}
    with pytest.raises(m.ReconcileMismatch):
        m.reconcile_raw_metrics(selected, tmp_path, "glob", "cipro", "Org", bad_artifact)


def test_reconcile_skips_accession_with_no_run_dir(tmp_path, monkeypatch):
    # An accession whose _run_dir resolves to None is skipped (not scored), not crashed.
    selected = {"a": 1, "b": 1}
    monkeypatch.setattr(m, "_run_dir", lambda acc, own, glob: None if acc == "b" else own / acc)
    monkeypatch.setattr(m, "call_resistance", lambda mt, drug, organism=None, **kw: {"prediction": "R"})
    artifact_metrics = {"tp": 1, "fp": 0, "tn": 0, "fn": 0, "sens": 1.0, "spec": None, "n_scored": 1}
    raw, preds = m.reconcile_raw_metrics(selected, tmp_path, "glob", "cipro", "Org", artifact_metrics)
    assert raw["n_scored"] == 1 and preds == {"a": "R"}  # "b" skipped entirely


# --------------------------------------------------------------------------- #
# pure metric assembly
# --------------------------------------------------------------------------- #
def test_graded_lineage_bucket_boundaries():
    assert m.graded_lineage_bucket(20).startswith("moderate")
    assert m.graded_lineage_bucket(15).startswith("moderate")
    assert m.graded_lineage_bucket(14).startswith("limited")
    assert m.graded_lineage_bucket(8).startswith("limited")
    assert m.graded_lineage_bucket(7).startswith("scarce")
    assert m.graded_lineage_bucket(3).startswith("scarce")
    assert m.graded_lineage_bucket(2).startswith("clonal")
    assert m.graded_lineage_bucket(0).startswith("clonal")


def test_build_threshold_results_divergence_and_ci():
    # clone of 10 R (pred R) + distinct R lineage missed -> weighted sens 0.5 with a CI.
    preds = {f"r{i}": "R" for i in range(10)} | {"rx": "S", "s0": "S"}
    labels = {f"r{i}": 1 for i in range(10)} | {"rx": 1, "s0": 0}
    clusters = {f"r{i}": 0 for i in range(10)} | {"rx": 1, "s0": 2}
    tr = m.build_threshold_results(preds, labels, {0.001: clusters})
    b = tr["0.001"]
    assert b["effective_lineage_N_R"] == 2 and b["effective_lineage_N_S"] == 1
    cw = b["cluster_weighted"]
    assert cw["sens"] == 0.5 and cw["sens_eff_n"] == 2
    assert isinstance(cw["sens_ci"], list) and len(cw["sens_ci"]) == 2  # C3: CI always present


def test_build_cell_scored_emits_grade_and_tier():
    tr = {"0.005": {"effective_lineage_N_R": 2, "effective_lineage_N_S": 8,
                    "cluster_weighted": {"sens": 0.5, "sens_ci": [0.1, 0.9], "sens_eff_n": 2,
                                         "spec": 1.0, "spec_ci": [0.6, 1.0], "spec_eff_n": 8,
                                         "n_discordant": 1}}}
    cell = m.build_cell(organism="Klebsiella", drug="ciprofloxacin", cohort="c",
                        raw={"n_scored": 60, "sens": 0.97, "spec": 0.97, "tp": 29, "fp": 1, "tn": 29, "fn": 1},
                        raw_reconciled=True, partial=False, n_genomes_missing=0, threshold_results=tr)
    assert cell["lineage_tier_emitted"] is True
    assert cell["lineage_grade"].startswith("clonal")  # 2 R lineages @0.005
    assert cell["raw_N"] == 60


def test_build_cell_partial_emits_no_tier():
    cell = m.build_cell(organism="Klebsiella", drug="tetracycline", cohort="c",
                        raw={"n_scored": 33, "sens": 0.8, "spec": 0.8, "tp": 1, "fp": 0, "tn": 1, "fn": 0},
                        raw_reconciled=True, partial=True, n_genomes_missing=27, threshold_results=None)
    assert cell["lineage_tier_emitted"] is False
    assert cell["thresholds"] == {} and cell["lineage_grade"] is None
    assert cell["n_genomes_missing"] == 27


def test_build_cell_unreconciled_emits_no_tier():
    cell = m.build_cell(organism="X", drug="cipro", cohort="c", raw={"n_scored": 10},
                        raw_reconciled=False, partial=False, n_genomes_missing=0, threshold_results={"0.005": {}})
    assert cell["lineage_tier_emitted"] is False and cell["thresholds"] == {}


# --------------------------------------------------------------------------- #
# idempotent sidecar upsert
# --------------------------------------------------------------------------- #
def test_upsert_replaces_same_canonical_key():
    sc = {"_schema": m.SCHEMA, "cells": []}
    sc = m.upsert_cell(sc, {"organism": "Klebsiella", "drug": "ciprofloxacin", "raw_N": 1})
    sc = m.upsert_cell(sc, {"organism": "klebsiella", "drug": "Ciprofloxacin", "raw_N": 2})  # same key
    assert len(sc["cells"]) == 1 and sc["cells"][0]["raw_N"] == 2


def test_upsert_appends_distinct_keys_sorted():
    sc = {"_schema": m.SCHEMA, "cells": []}
    sc = m.upsert_cell(sc, {"organism": "Klebsiella", "drug": "tetracycline", "raw_N": 1})
    sc = m.upsert_cell(sc, {"organism": "Campylobacter", "drug": "ciprofloxacin", "raw_N": 1})
    keys = [(c["organism"], c["drug"]) for c in sc["cells"]]
    assert keys == [("Campylobacter", "ciprofloxacin"), ("Klebsiella", "tetracycline")]


def test_sidecar_roundtrip(tmp_path):
    p = tmp_path / "lineage.json"
    sc = m.load_sidecar(p)  # absent -> fresh
    assert sc["_schema"] == m.SCHEMA and sc["cells"] == []
    sc = m.upsert_cell(sc, {"organism": "Klebsiella", "drug": "ciprofloxacin", "raw_N": 60})
    m.write_sidecar(sc, p)
    reloaded = m.load_sidecar(p)
    assert reloaded["cells"][0]["raw_N"] == 60


def test_load_sidecar_malformed_returns_fresh(tmp_path):
    # A corrupt sidecar must not crash a re-run — it falls back to a fresh skeleton.
    p = tmp_path / "lineage.json"
    p.write_text("{ this is not json", encoding="utf-8")
    sc = m.load_sidecar(p)
    assert sc["_schema"] == m.SCHEMA and sc["cells"] == []


def test_load_sidecar_existing_without_cells_key(tmp_path):
    # A valid JSON missing the "cells" key gets it defaulted (setdefault path).
    p = tmp_path / "lineage.json"
    p.write_text(json.dumps({"_schema": m.SCHEMA, "date": "2026-06-11"}), encoding="utf-8")
    sc = m.load_sidecar(p)
    assert sc["cells"] == []


# --------------------------------------------------------------------------- #
# process_cohort guard paths (no Docker / network: only the early-exit branches)
# --------------------------------------------------------------------------- #
def test_process_cohort_no_artifact_raises(tmp_path, monkeypatch):
    """A provdisjoint cohort dir with no committed validation JSON raises FileNotFoundError
    BEFORE any Mash/genome work (the artifact is the reconciliation source of truth — M4)."""
    cd = tmp_path / "klebsiella_provdisjoint_ciprofloxacin"
    cd.mkdir()
    monkeypatch.setattr(m, "find_artifact", lambda slug, drug: None)
    with pytest.raises(FileNotFoundError):
        m.process_cohort(cd, use_docker=False)
