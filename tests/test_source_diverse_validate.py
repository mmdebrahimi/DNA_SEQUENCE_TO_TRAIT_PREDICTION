"""Guards for the source-diverse validation arm.

The arm's argument is "your cohort was too concentrated to see this". An arm making that argument must
apply the standard to itself, or it is just a second opinion from a differently-biased sample.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))


def test_the_arm_refuses_a_cohort_as_concentrated_as_the_ones_it_critiques():
    """A 1-BioProject cohort must be refused outright, not reported with a caveat."""
    from source_diverse_validate import diversity_verdict
    ok, why = diversity_verdict({"n": 60, "n_known": 60, "distinct": 1,
                                 "largest_share": 1.0, "dominant": "P"})
    assert not ok and "1 BioProject" in why


def test_the_bar_catches_share_even_when_source_COUNT_looks_fine():
    """5 sources sounds diverse until one holds 68% -- the case that retracted a published number."""
    from source_diverse_validate import diversity_verdict
    ok, why = diversity_verdict({"n": 66, "n_known": 66, "distinct": 5,
                                 "largest_share": 0.682, "dominant": "P"})
    assert not ok and "68%" in why


def test_mostly_unknown_provenance_is_refused_not_treated_as_diverse():
    """Missing metadata makes a cohort LOOK diverse. The arm must not read absence as spread."""
    from source_diverse_validate import diversity_verdict
    ok, why = diversity_verdict({"n": 100, "n_known": 20, "distinct": 20,
                                 "largest_share": 0.05, "dominant": "P"})
    assert not ok and "unknown" in why


def test_a_genuinely_diverse_cohort_passes():
    from source_diverse_validate import diversity_verdict
    ok, why = diversity_verdict({"n": 131, "n_known": 131, "distinct": 8,
                                 "largest_share": 0.305, "dominant": "P"})
    assert ok and why == "source-diverse"


def test_results_use_a_SEPARATE_namespace_from_the_frozen_cells():
    """`provenance_disjoint_validation_*` is the glob the report card's load_scored() reads. Writing
    there would silently overwrite a frozen cell -- the shared-key trap."""
    src = (ROOT / "scripts" / "source_diverse_validate.py").read_text(encoding="utf-8")
    assert "source_diverse_validation_" in src
    assert 'f"provenance_disjoint_validation_' not in src
    for f in (ROOT / "wiki").glob("source_diverse_validation_*.json"):
        d = json.loads(f.read_text(encoding="utf-8"))
        assert d["schema"] == "source-diverse-validation-v1"


def test_refused_cells_carry_no_metrics():
    """A refused cell must not leak an acc/sens/spec a reader could quote anyway."""
    files = list((ROOT / "wiki").glob("source_diverse_validation_*.json"))
    if not files:
        import pytest
        pytest.skip("no artifacts yet")
    refused = [json.loads(f.read_text(encoding="utf-8")) for f in files]
    refused = [d for d in refused if d.get("status") == "source_concentrated"]
    assert refused, "expected at least one refusal -- the bar is not exercised otherwise"
    for d in refused:
        assert "acc" not in d and "sens" not in d and "spec" not in d
