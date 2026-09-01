"""Pin the NNRTI curation attempt, and the defect that would have flattered its headline.

F-B measured whether curating the shipped NNRTI catalog beats the free position-novelty doubt layer.
It does not. These tests pin the parser defect that surfaced en route, plus the method contract, so a
future rerun cannot quietly regress into the flattering version.

Offline and pure for the parser tests; the data-dependent ones skip when the gitignored dataset is absent.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.hiv_nnrti_mutant_catalog import (  # noqa: E402
    DEFAULT_DATA, MIN_CARRIERS, RESIST_COEF_MIN, parse_comp_mut_list,
)


# --- the defect: self-to-self entries are not mutations ---

def test_self_to_self_entries_are_not_mutations():
    """THE DEFECT. CompMutList carries tokens like `L234L` / `K238K` / `M230M` / `R72R` where WT == MUT.
    They encode a mixture containing wild-type or an ambiguity, not a substitution. Admitting them let
    the OLS assign real coefficients to what is effectively a sequencing/mixture-status marker -- which
    correlates with treatment experience, so it is a pure confound. Four landed in the first derived
    catalog, and excluding them changed the EFV result (8 additions -> 5, balacc 0.929 -> 0.932). The
    headline would have been wrong in the FLATTERING direction."""
    assert parse_comp_mut_list("L234L, K238K, M230M, R72R") == set()


def test_real_substitutions_still_parse():
    """NON-VACUITY: the filter must not swallow genuine mutations."""
    assert parse_comp_mut_list("D67N, K103N, V118I") == {"D67N", "K103N", "V118I"}


def test_a_mixture_keeps_the_non_wildtype_members_only():
    """`K103NS` = a mixture of N and S at 103. Both are real substitutions; a WT member is dropped."""
    assert parse_comp_mut_list("K103NS") == {"K103N", "K103S"}
    assert parse_comp_mut_list("K103KN") == {"K103N"}      # the K member is WT -> excluded


def test_malformed_tokens_are_ignored_not_guessed():
    assert parse_comp_mut_list("") == set()
    assert parse_comp_mut_list("garbage, 103N, K103") == set()


# --- the method contract ---

def test_the_threshold_is_derived_not_inherited():
    """3x, NOT the 1.5x the NRTI/PI/INSTI builders use. Those restrict candidates to a class's
    catalogued positions; this one spans the whole RT, so the multiple-comparisons burden is far larger.
    At 1.5x it admits 43 mutations and collapses EFV specificity 0.904 -> 0.691."""
    assert RESIST_COEF_MIN == pytest.approx(math.log10(3.0))
    assert MIN_CARRIERS == 5


def test_candidates_span_the_whole_rt_not_just_catalogued_positions():
    """The blind spot's drivers sit OUTSIDE the 8 catalogued positions, so a position-restricted
    candidate set could never reach them. CompMutList reaches the connection domain (>318)."""
    subs = parse_comp_mut_list("K103N, V179D, N348I, G335D")
    assert {int(s[1:-1]) for s in subs} == {103, 179, 348, 335}


# --- the verdict, pinned against the committed artifact ---

def test_the_curation_loses_to_the_free_doubt_layer():
    """THE VERDICT. If a rerun ever shows blind-spot recovery ABOVE the position-novelty flag's 0.604,
    that is a real finding and this must fail loudly rather than let the negative stand unexamined."""
    import json
    hits = sorted((ROOT / "wiki").glob("hiv_nnrti_mutant_catalog_*.json"))
    if not hits:
        pytest.skip("NNRTI catalog artifact not generated on this checkout")
    d = json.loads(hits[-1].read_text(encoding="utf-8"))
    scored = [(k, m) for k, m in d["per_drug"].items() if "blind_spot" in m]
    assert scored, "no drug was scored — the artifact is degenerate"
    for drug, m in scored:
        rr = m["blind_spot"]["recovery_rate"]
        if rr is None:
            continue
        assert rr <= 0.604, (
            f"{drug} blind-spot recovery {rr} EXCEEDS the free position-novelty incumbent 0.604 — "
            "re-examine the verdict rather than keeping the negative")


def test_doravirine_is_a_reported_wall_not_a_guessed_cutoff():
    """DOR postdates Stanford DRMcv.R, so no clinical cutoff is sourced. Scoring it at a guessed
    boundary would manufacture a number."""
    import json
    hits = sorted((ROOT / "wiki").glob("hiv_nnrti_mutant_catalog_*.json"))
    if not hits:
        pytest.skip("artifact absent")
    dor = json.loads(hits[-1].read_text(encoding="utf-8"))["per_drug"].get("doravirine", {})
    assert dor.get("status") == "CUTOFF_UNAVAILABLE"
    assert "blind_spot" not in dor


def test_the_shipped_catalog_was_not_modified():
    """F-B's outcome is a NO-SHIP. `hiv_amr.py` must be untouched: 16 entries, 8 positions."""
    from dna_decode.data.hiv_amr import NNRTI_RT_MAJOR_DRMS
    assert len(NNRTI_RT_MAJOR_DRMS) == 16
    assert {int(s[1:-1]) for s in NNRTI_RT_MAJOR_DRMS} == {100, 101, 103, 106, 181, 188, 190, 230}
    assert "V179D" not in NNRTI_RT_MAJOR_DRMS, "a curation was applied without re-running the verdict"
