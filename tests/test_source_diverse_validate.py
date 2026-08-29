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


def test_effective_sources_catches_a_cohort_that_passes_count_and_share():
    """THE LOOPHOLE, found on real data. A cohort of 2 real BioProjects plus 3 token ones passes both
    the count bar (5) and the share bar (45% <= 60%) while being 2 sources wearing 5 badges.
    Inverse-Simpson scores it 2.45 vs 4.44 for a genuinely spread cohort.

    Not hypothetical: E. coli x meropenem in the live pool has 5 nominal sources at 48% share and 2.42
    effective, and this rule is what moved it from `underpowered` to `source_concentrated`.
    """
    from source_diverse_validate import effective_sources, diversity_verdict
    assert round(effective_sources([18, 18, 2, 1, 1]), 2) == 2.45
    assert round(effective_sources([12, 10, 8, 6, 4]), 2) == 4.44
    gameable = {"n": 40, "n_known": 40, "distinct": 5, "largest_share": 0.45,
                "dominant": "P", "effective_sources": 2.45}
    ok, why = diversity_verdict(gameable)
    assert not ok and "effective sources" in why
    honest = {"n": 40, "n_known": 40, "distinct": 5, "largest_share": 0.30,
              "dominant": "P", "effective_sources": 4.44}
    assert diversity_verdict(honest)[0]


def test_effective_sources_is_degenerate_safe():
    from source_diverse_validate import effective_sources
    assert effective_sources([]) == 0.0
    assert effective_sources([0, 0]) == 0.0
    assert round(effective_sources([10]), 2) == 1.0     # one source is one source
