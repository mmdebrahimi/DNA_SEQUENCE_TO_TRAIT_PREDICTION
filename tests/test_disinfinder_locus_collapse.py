"""The ResFinder locus-collapse fix is INERT for DisinFinder — measured, not assumed.

DisinFinder shares the defective pattern literally (it imports `resfinder.gene_of` and keys on the
allele name), so propagating the fix looked obvious. This project's standing rule is that a shared code
pattern is a LEAD, not a diagnosis. These tests pin the measurement and, importantly, that the probe
REFUSES to call inertness when nothing was ever called.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from disinfinder_locus_collapse_probe import new_rule, old_rule  # noqa: E402

ARTIFACT = ROOT / "wiki" / "disinfinder_locus_collapse_probe_2026-09-05.json"


def hit(pid: float, cov: float, contig: str = "c1", s: int = 100, e: int = 900) -> dict:
    return {"called": True, "percent_identity": pid, "percent_coverage": cov,
            "contig": contig, "sstart": s, "send": e}


def test_the_two_rules_agree_when_there_is_no_dense_variant_family():
    """One allele per locus is the DisinFinder shape; both rules must return the same gene."""
    called = [("qacE_1_X", hit(100.0, 100.0)), ("sitABCD_1_Y", hit(99.0, 100.0, contig="c2"))]
    assert old_rule(called) == new_rule(called)


def test_the_rules_DIVERGE_on_a_dense_family_so_the_probe_can_detect_one():
    """The negative result is only meaningful if the probe COULD have found a difference."""
    called = [("qacE_1_X", hit(100.0, 100.0)), ("qacEdelta1_1_Y", hit(99.2, 100.0))]
    assert len(old_rule(called)) == 2       # allele-name keying reports both
    assert len(new_rule(called)) == 1       # one locus -> one gene
    assert old_rule(called) != new_rule(called)


def test_new_rule_keeps_genuinely_separate_loci_apart():
    called = [("qacE_1_X", hit(100.0, 100.0, contig="c1")),
              ("sitABCD_1_Y", hit(100.0, 100.0, contig="c2"))]
    assert len(new_rule(called)) == 2


@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("disinfinder probe artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_artifact_verdict_is_inert(art):
    assert art["verdict"] == "LOCUS_COLLAPSE_IS_INERT_FOR_DISINFINDER"
    assert art["n_genomes_with_different_gene_set"] == 0


def test_artifact_is_non_vacuous(art):
    """Both rules agree on the empty set if nothing was ever called; the probe must have done work."""
    assert art["n_genomes_scored"] > 0
    assert art["total_gene_calls_old_rule"] > 0


def test_the_inertness_is_scoped_to_this_db_build(art):
    """The DB size is the mechanism; a bigger catalogue could behave like the beta-lactamases."""
    assert art["n_db_alleles"] < 50
    assert any("db build" in lim.lower() or "re-measur" in lim.lower()
               for lim in art["honest_limits"])
