"""The multi-hit probe must REFUSE when its parse returns nothing, not report 'inconclusive'.

This guards a defect that shipped in the sibling anatomy probe and, per a peer session, recurred four
times in 48 hours across two projects: an empty PARSE RESULT is indistinguishable from an empty INPUT,
and the empty case is always the reassuring one. A tool that reports "nothing to report" when its
parser failed is worse than one that crashes, because the output is SHAPED like success.

Every test here drives the real `score_cell` with its collaborators patched, so the guard is exercised
rather than asserted.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import cabbage_qrdr_multihit as M  # noqa: E402


class _FakePath:
    """Stands in for the on-disk cache dir. `exists()` decides whether a TSV is 'cached'."""

    def __init__(self, p: str, cached: bool):
        self._p, self._cached = p, cached

    def exists(self) -> bool:
        return self._cached

    def __truediv__(self, other):
        return _FakePath(f"{self._p}/{other}", self._cached)

    def __str__(self) -> str:
        return self._p


def _patch(monkeypatch, *, labels, determinants, cached=True):
    """Point score_cell's collaborators at synthetic data.

    score_cell imports build_labels and CACHE INSIDE the function, so patching the module attributes
    here takes effect at call time.
    """
    import dna_decode.data.amr_portal as AP
    import dna_decode.eval.amr_rules as AR
    import scripts.cabbage_score_cell as CS

    monkeypatch.setattr(CS, "build_labels",
                        lambda species, drug: (labels, {"tier_census": {}}, {}, {}, {}))
    monkeypatch.setattr(CS, "CACHE", _FakePath("D:/fake", cached))
    monkeypatch.setattr(AP, "prepare", lambda raw, out: "fake.tsv")
    monkeypatch.setattr(AR, "call_resistance",
                        lambda path, drug, organism: {"prediction": "R",
                                                      "determinants": list(determinants),
                                                      "n_determinants": len(determinants)})


def test_zero_labels_REFUSES_rather_than_reporting_an_empty_strata_table(monkeypatch):
    """A mis-typed species or drug string yields no labels; without this the verdict reads
    'inconclusive' over a table that was never populated."""
    _patch(monkeypatch, labels={}, determinants=[{"symbol": "gyrA_T86I"}])
    with pytest.raises(M.ParseFailure) as e:
        M.score_cell("Nonexistent species", "Salmonella", "ciprofloxacin")
    assert "ZERO labels" in str(e.value)


def test_zero_determinant_carriers_REFUSES_because_these_cells_have_hundreds(monkeypatch):
    """The exact shape that shipped in the anatomy probe: the record parses, every bucket comes back
    determinant-free, and the verdict is confident and wrong."""
    _patch(monkeypatch, labels={f"GCA_{i}": "R" for i in range(50)}, determinants=[])
    with pytest.raises(M.ParseFailure) as e:
        M.score_cell("Salmonella enterica", "Salmonella", "ciprofloxacin")
    msg = str(e.value)
    assert "NOT ONE carries a counted determinant" in msg
    assert "not a finding" in msg


def test_nothing_scoreable_REFUSES(monkeypatch):
    """Labels exist but no cached TSV does -- nothing was measured, so there is nothing to report."""
    _patch(monkeypatch, labels={f"GCA_{i}": "R" for i in range(50)},
           determinants=[{"symbol": "gyrA_T86I"}], cached=False)
    with pytest.raises(M.ParseFailure) as e:
        M.score_cell("Salmonella enterica", "Salmonella", "ciprofloxacin")
    assert "ZERO isolates scored" in str(e.value)


def test_a_populated_parse_does_NOT_fire_the_guard(monkeypatch):
    """Non-vacuity the other way: the guard must not fire on a real measurement, or it is just a
    crash. A single-determinant carrier set is a legitimate result (it is what Campylobacter is)."""
    _patch(monkeypatch, labels={f"GCA_{i}": ("R" if i % 2 else "S") for i in range(40)},
           determinants=[{"symbol": "gyrA_T86I"}])
    out = M.score_cell("Campylobacter jejuni", "Campylobacter", "ciprofloxacin")
    assert out["by_determinant_count"], "a populated parse produced no strata"
    assert any(r["n_determinants"] == 1 for r in out["by_determinant_count"])


def test_the_guard_is_a_domain_invariant_and_says_so():
    """The refusal rests on these cells being known to contain hundreds of carriers -- a property of
    the data, not a heuristic threshold. If that reasoning is ever dropped the guard becomes
    arbitrary, so it is pinned here."""
    doc = M.ParseFailure.__doc__ or ""
    assert "DOMAIN INVARIANT" in doc
    assert "reassuring" in doc
