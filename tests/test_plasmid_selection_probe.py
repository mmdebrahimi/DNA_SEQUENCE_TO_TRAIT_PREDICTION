"""The plasmid null must come from the rule being exercised, not from nothing being called.

A "no difference" result is only meaningful if both orderings actually ran on real calls. These pin the
non-vacuity, the structural claim, and the caveat that secondary fields DO move.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "plasmid_selection_rule_probe.py"
ART = ROOT / "wiki" / "plasmid_selection_rule_probe_2026-09-04.json"


def test_best_per_family_orders_by_the_requested_axis():
    """The two orderings must genuinely differ when identity and coverage disagree."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from plasmid_selection_rule_probe import best_per_family
    # same family, two alleles: one wins on coverage, the other on identity
    per_allele = {
        "IncFIA_1__x": {"called": True, "percent_identity": 99.9, "percent_coverage": 70.0},
        "IncFIA_2__y": {"called": True, "percent_identity": 90.0, "percent_coverage": 100.0},
    }
    cov = best_per_family(per_allele, identity_primary=False)
    ident = best_per_family(per_allele, identity_primary=True)
    fam = next(iter(cov))
    assert cov[fam]["best_allele"] != ident[fam]["best_allele"], (
        "if the two orderings never disagree the probe would be vacuous")
    assert set(cov) == set(ident), "the FAMILY set is invariant -- that is the structural claim"


def test_uncalled_alleles_never_enter_either_ordering():
    sys.path.insert(0, str(ROOT / "scripts"))
    from plasmid_selection_rule_probe import best_per_family
    per_allele = {"IncX1_1__z": {"called": False, "percent_identity": 100.0,
                                 "percent_coverage": 100.0}}
    assert best_per_family(per_allele, identity_primary=True) == {}


def test_probe_refuses_when_nothing_was_called(tmp_path):
    """A null from an empty run is a plumbing result, not a finding."""
    empty = tmp_path / "asm"
    (empty / "GCA_000000000.0").mkdir(parents=True)
    out = tmp_path / "o.json"
    r = subprocess.run([sys.executable, str(SCRIPT), "--asm-dirs", str(empty),
                        "--out", str(out)], capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 2, r.stdout + r.stderr   # no assemblies found at all
    assert not out.exists()


@pytest.mark.skipif(not ART.exists(), reason="probe artifact absent")
def test_the_null_is_non_vacuous():
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["n_assemblies_with_different_replicon_set"] == 0
    assert d["total_replicon_calls"] > 100, "the null must rest on real replicon calls"
    assert d["n_assemblies_with_moved_best_allele"] > 0, (
        "if the ordering never moved ANYTHING the probe never exercised the rule")


@pytest.mark.skipif(not ART.exists(), reason="probe artifact absent")
def test_the_secondary_exposure_is_stated_not_buried():
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["verdict"] == "STRUCTURALLY_INERT_FOR_THE_REPORTED_SET", (
        "the verdict must be scoped to the reported SET, not claim blanket inertness")
    assert any("Secondary fields DO move" in s or "secondary" in s.lower()
               for s in d["honest_limits"])


@pytest.mark.skipif(not ART.exists(), reason="probe artifact absent")
def test_the_prediction_was_recorded_and_marked_survived():
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert "structural_prediction" in d and d["prediction_survived"] is True
