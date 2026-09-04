"""The sweep's negative must come from looking, not from failing to look.

A sweep that examines nothing returns the same verdict as a sweep that examines everything and finds
nothing. These pin the difference: it must actually compare multiple organisms, and it must still flag
a genuine organism-dependent gap when one is planted.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "organism_scope_sweep.py"
ARTIFACT = ROOT / "wiki" / "organism_scope_sweep_2026-09-04.json"


def _cell(org, drug, spec, share, bp=10):
    return {"organism": org, "drug": drug, "state": "SCORED", "sens": 0.9, "spec": spec, "n": 60,
            "source_concentration": {"largest_share": share, "distinct_bioprojects": bp}}


def _run(cells, tmp_path):
    card = tmp_path / "card.json"
    card.write_text(json.dumps({"cells": cells}), encoding="utf-8")
    out = tmp_path / "o.json"
    r = subprocess.run([sys.executable, str(SCRIPT), "--card", str(card), "--out", str(out)],
                       capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 0, r.stdout + r.stderr
    return json.loads(out.read_text(encoding="utf-8"))


def test_it_flags_a_planted_organism_dependent_gap(tmp_path):
    """A big spec spread in DIVERSE cohorts is a real gap and must not be explained away."""
    d = _run([_cell("escherichia_coli_shigella", "tetracycline", 0.98, 0.20),
              _cell("klebsiella", "tetracycline", 0.40, 0.20)], tmp_path)
    assert d["verdict"] == "ORGANISM_SCOPE_GAP_FOUND"
    f = d["findings"][0]
    assert f["flagged"] and not f["attributable_to_concentration"]


def test_a_spread_carried_by_a_concentrated_cohort_is_attributed_not_alarmed(tmp_path):
    """The same spread, but one cell fails the diversity bar -> evidence about the cohort, not the rule."""
    d = _run([_cell("escherichia_coli_shigella", "tetracycline", 0.98, 0.97),
              _cell("klebsiella", "tetracycline", 0.40, 0.20)], tmp_path)
    assert d["findings"][0]["attributable_to_concentration"] is True
    assert d["verdict"] == "ISOLATED_TO_GENTAMICIN"


def test_a_single_organism_drug_cannot_produce_a_finding(tmp_path):
    """No comparison is possible; the sweep must not manufacture one."""
    d = _run([_cell("klebsiella", "meropenem", 0.40, 0.20)], tmp_path)
    assert d["n_multi_organism_drugs"] == 0 and d["findings"] == []


@pytest.mark.skipif(not ARTIFACT.exists(), reason="committed sweep artifact absent")
def test_the_committed_negative_is_non_vacuous():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["verdict"] == "ISOLATED_TO_GENTAMICIN"
    # It must have actually compared something -- a zero from an empty sweep is not a finding.
    assert d["n_multi_organism_drugs"] >= 3
    assert len(d["findings"]) >= 3
    assert any(f["flagged"] for f in d["findings"]), (
        "if nothing was ever flagged the threshold may be unreachable, which would make the "
        "negative vacuous")
    # And the honest limit that a single-organism drug proves nothing must be stated.
    assert any("ONE organism" in s for s in d["what_this_does_not_show"])
