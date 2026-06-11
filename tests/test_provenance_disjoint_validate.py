"""Pin the provenance-disjoint validator's main() control flow (scripts/provenance_disjoint_validate.py).

The load-bearing behavior added in this plan is the FAIL-CLOSED manifest gate (C1): an incomplete accession
manifest must REFUSE to score (exit 2) unless explicitly overridden, because an incomplete manifest is a
possible false-independence claim, not a harmless offline fallback. These tests drive main() with the
manifest + selection stages monkeypatched so no network / Docker is touched — the gate fires before any I/O.
"""
from __future__ import annotations

import importlib
import sys

import pytest

import scripts.provenance_disjoint_validate as pdv
from dna_decode.eval.cohort_manifest import Manifest


def _run_main(argv, monkeypatch):
    """Drive main() with a fixed argv (no network: select stage is stubbed by the caller)."""
    monkeypatch.setattr(sys, "argv", ["provenance_disjoint_validate.py", *argv])
    return pdv.main()


def test_main_fail_closed_on_incomplete_manifest(monkeypatch, capsys):
    """C1: incomplete manifest + no override -> exit 2, refuses to score, never reaches selection/network."""
    monkeypatch.setattr(pdv, "build_manifest",
                        lambda *a, **k: Manifest(cohorts=[], incomplete=True,
                                                 warnings=["parquet load failed: x.parquet: OSError: boom"]))

    def _boom(*a, **k):  # selection must NOT run when the gate fails closed
        raise AssertionError("select_disjoint reached despite fail-closed manifest")

    monkeypatch.setattr(pdv, "select_disjoint", _boom)
    rc = _run_main(["--group", "Klebsiella", "--drug", "ciprofloxacin"], monkeypatch)
    assert rc == 2
    out = capsys.readouterr().out
    assert "INCOMPLETE_MANIFEST" in out
    assert "x.parquet" in out  # the warning is surfaced to the operator


def test_main_override_proceeds_past_gate(monkeypatch):
    """--allow-incomplete-manifest lets an incomplete manifest proceed (DEGRADED) — gate does not exit 2.
    Stubbed selection + --select-only returns 0 before any AMRFinder/network step."""
    monkeypatch.setattr(pdv, "build_manifest",
                        lambda *a, **k: Manifest(cohorts=[], incomplete=True, warnings=["degraded"]))
    captured = {}

    def _sel(group, drug, per_class, reuse_glob, selected, exclude_prior):
        captured["exclude_prior"] = exclude_prior
        return {"GCA_1.1": 1, "GCA_2.1": 0}

    monkeypatch.setattr(pdv, "select_disjoint", _sel)
    rc = _run_main(["--group", "Klebsiella", "--drug", "ciprofloxacin",
                    "--allow-incomplete-manifest", "--select-only"], monkeypatch)
    assert rc == 0
    assert captured["exclude_prior"] == set()  # empty manifest -> nothing to exclude


def test_main_complete_manifest_passes_prior_accessions(monkeypatch):
    """Complete manifest -> no exit 2; the exclusion set from prior_accessions reaches select_disjoint, and
    the current output cohort (exact-self) is NOT in it."""
    from dna_decode.eval.cohort_manifest import Cohort

    base_name = "klebsiella_provdisjoint_ciprofloxacin"
    manifest = Manifest(cohorts=[
        Cohort("klebsiella_cipro", "p1", "calibration", "selected_tsv", "klebsiella", "ciprofloxacin",
               {"GCA_prior.1", "GCA_prior.2"}),
        Cohort(base_name, "p2", "validation", "selected_tsv", "klebsiella", "ciprofloxacin", {"GCA_self.1"}),
    ], incomplete=False)
    monkeypatch.setattr(pdv, "build_manifest", lambda *a, **k: manifest)
    seen = {}

    def _sel(group, drug, per_class, reuse_glob, selected, exclude_prior):
        seen["exclude"] = exclude_prior
        return {"GCA_x.1": 1}

    monkeypatch.setattr(pdv, "select_disjoint", _sel)
    rc = _run_main(["--group", "Klebsiella", "--drug", "ciprofloxacin", "--select-only"], monkeypatch)
    assert rc == 0
    assert seen["exclude"] == {"GCA_prior.1", "GCA_prior.2"}  # prior cohort accessions excluded
    assert "GCA_self.1" not in seen["exclude"]                # exact-self cohort NOT excluded
