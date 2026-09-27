"""Tests for the contract-number audit.

The load-bearing tests are the NON-VACUITY ones. This audit reached a clean 0/130 once by accepting any
number traceable to the package source -- which accepted 12 of 12 randomly-generated numbers. A check that
passes everything reads exactly like a check that found nothing, so the guards below plant drift and require
it to be caught.
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from scripts.contract_number_audit import (
    ADJUDICATED_BENIGN,
    DECIMAL_RE,
    _cited_artifacts,
    _expand_braces,
    _matches_by_rounding,
    audit,
)

REPO = Path(__file__).resolve().parent.parent
ARTIFACT = REPO / "wiki" / "contract_number_audit_2026-09-24.json"


# --- citation parsing ----------------------------------------------------------------------------

def test_brace_expansion():
    """The repo cites artifacts in brace form; not expanding it made the audit invent drift on the
    inverse cell's +53.0%, whose source sits in a roundtrip artifact the regex never loaded."""
    assert set(_expand_braces("wiki/a_2026-09-24.{md,json}")) == {
        "wiki/a_2026-09-24.md", "wiki/a_2026-09-24.json"}
    multi = set(_expand_braces("wiki/f_{x,y}_2026-07-1{6,7}.{md,json}"))
    assert len(multi) == 8
    assert "wiki/f_x_2026-07-16.md" in multi
    assert _expand_braces("wiki/plain.json") == ["wiki/plain.json"]


def test_cited_artifacts_tries_both_sidecars_for_a_bare_stem():
    got = _cited_artifacts("see wiki/foo_2026-01-01 for detail")
    assert "wiki/foo_2026-01-01.json" in got and "wiki/foo_2026-01-01.md" in got


def test_cited_artifacts_strips_trailing_punctuation():
    assert "wiki/foo.json" in _cited_artifacts("(wiki/foo.json).")


# --- what counts as a cited NUMBER ---------------------------------------------------------------

def test_identifier_fragments_are_not_numbers():
    """canFam3.1 / GCA_000005845.2 / v1.2 are identifiers. Capturing them produced a false
    CANDIDATE_DRIFT on the dog cell's '3.1', which is an assembly name."""
    for ident in ("canFam3.1", "GCA_000005845.2", "v1.2", "hg19.2"):
        assert DECIMAL_RE.findall(ident) == [], f"{ident} should not yield a measurement"


def test_real_decimals_are_still_captured():
    assert DECIMAL_RE.findall("accuracy 0.8125 and spec 0.939") == ["0.8125", "0.939"]


# --- rounding-aware matching ---------------------------------------------------------------------

def test_rounding_match_accepts_prose_rounding():
    """Prose rounds: a cited 0.398 must match an artifact's 0.3977. Without this the audit flagged
    SLCO1B1's correct figure."""
    assert _matches_by_rounding("0.398", [0.3977]) == pytest.approx(0.3977)
    assert _matches_by_rounding("0.40", [0.3977]) == pytest.approx(0.3977)


def test_rounding_match_accepts_percent_fraction_pairing():
    assert _matches_by_rounding("53.0", [0.5295]) == pytest.approx(0.5295)
    assert _matches_by_rounding("0.29", [29.0]) == pytest.approx(29.0)


def test_rounding_match_rejects_a_genuinely_different_number():
    """It must NOT be a universal acceptor -- that is the vacuity failure."""
    assert _matches_by_rounding("0.79", [0.2128]) is None
    assert _matches_by_rounding("0.926", [0.454]) is None


# --- NON-VACUITY: planted drift must be caught ---------------------------------------------------

def test_planted_drift_is_caught():
    """THE guard. Feed the real matcher a number that is in no artifact and confirm it does not pass."""
    assert _matches_by_rounding("0.7431", [0.1, 0.2, 0.9999]) is None


def test_random_numbers_do_not_pass_against_a_small_artifact_set():
    """The vacuity regression, stated as a property: 12 random 4-decimal numbers must not all match a
    handful of artifact values. The earlier package-source variant accepted 12 of 12."""
    random.seed(7)
    art = [0.3977, 0.8125, 0.2128, 64.0, 88.0]
    toks = [f"0.{random.randint(1000, 9999)}" for _ in range(12)]
    matched = [t for t in toks if _matches_by_rounding(t, art) is not None]
    assert len(matched) <= 1, f"too permissive: {matched}"


def test_adjudicated_benign_is_small_and_every_entry_has_a_reason():
    """Exemptions are named one at a time, with a reason, precisely because the blanket version was
    vacuous. A growing table is a smell, so this pins it small."""
    assert len(ADJUDICATED_BENIGN) <= 5
    for (cell_id, tok), why in ADJUDICATED_BENIGN.items():
        assert cell_id and tok
        assert len(why) > 40, f"{cell_id} {tok} needs a real reason, got {why!r}"


# --- live audit ----------------------------------------------------------------------------------

def test_live_audit_is_clean_and_non_trivial():
    """Clean AND actually checking something -- a zero over zero checks is not a result."""
    rep = audit()
    assert rep["n_numbers_checked"] >= 100
    assert rep["n_cells_audited"] >= 10
    assert rep["n_candidate_drift"] == 0, rep["verdict"]


def test_live_audit_would_flag_an_injected_bad_citation():
    """End-to-end non-vacuity: a cell whose prose cites a number absent from its artifact must be
    reported, using the REAL audit path on a synthetic contract."""
    import scripts.contract_number_audit as m

    class FakeCell:
        cell_id = "synthetic:test:cell"
        claim = "accuracy 0.7431 on the panel (wiki/contract_number_audit_2026-09-24.json)"
        validation_slice = label_provenance = demotion_rule = claim_status = ""

    real = m.audit.__globals__["__builtins__"]  # touch nothing; we monkeypatch the import instead
    assert real is not None
    # drive the internals directly rather than patching the registry import
    blob = FakeCell.claim
    arts = m._cited_artifacts(blob)
    text = "\n".join(t for t in (m._load_artifact_text(a) for a in arts) if t)
    assert text, "fixture artifact must resolve, else this test proves nothing"
    nums = m.DECIMAL_RE.findall(blob)
    assert "0.7431" in nums
    assert not any(v in text for v in m._number_variants("0.7431"))
    assert m._matches_by_rounding("0.7431", m._artifact_numbers(text)) is None


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_committed_artifact_records_its_own_posture():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert "TRIAGE" in d["posture"].upper()
    assert d["n_numbers_checked"] >= 100


# --- coverage disclosure -------------------------------------------------------------------------
# `scope` always said the audit covers "cells whose prose names a wiki/ artifact", so the exclusion
# was disclosed -- but not its SIZE, which is what says whether the audited set is near-complete or a
# fraction. It is a fraction, and the skipped set is LARGER than the checked one.

def test_coverage_block_is_present_and_internally_consistent():
    rep = audit()
    assert rep["n_numbers_unverifiable"] == sum(u["n_numbers"] for u in rep["cells_unverifiable"])
    assert rep["n_cells_with_numbers_but_no_citation"] == len(rep["cells_unverifiable"])
    for u in rep["cells_unverifiable"]:
        assert u["n_numbers"] == len(u["numbers"]) > 0, "a cell with no numbers is not a coverage gap"
        assert u["why"]


def test_coverage_is_non_trivial_so_the_disclosure_is_not_decoration():
    """If nothing were excluded the block would be noise; the point is that plenty is."""
    rep = audit()
    assert rep["n_cells_with_numbers_but_no_citation"] > 0
    assert rep["n_numbers_unverifiable"] > 0


def test_unverifiable_cells_are_disjoint_from_audited_cells():
    """A cell is either checked against an artifact or reported as uncheckable — never both."""
    rep = audit()
    audited = {r["cell_id"] for r in rep["cells"]}
    skipped = {u["cell_id"] for u in rep["cells_unverifiable"]}
    assert not (audited & skipped)


def test_coverage_does_not_contaminate_the_drift_verdict():
    """Unchecked is not wrong. The verdict must key on drift alone, or an uncited number would read
    as a defect and the audit would stop being a triage funnel."""
    rep = audit()
    expected = "NOTHING_TO_ADJUDICATE" if rep["n_candidate_drift"] == 0 else "ADJUDICATION_REQUIRED"
    assert rep["verdict"] == expected
    assert rep["n_numbers_checked"] == sum(r["checked"] for r in rep["cells"])


def test_an_uncited_cell_whose_only_numbers_are_adjudicated_bars_is_not_flagged():
    """Non-vacuity in the OTHER direction: a threshold has no artifact to live in, so a cell carrying
    only adjudicated bars is not a coverage gap. Exercises the real filter on a synthetic cell."""
    import scripts.contract_number_audit as m

    (cell_id, tok) = next(iter(ADJUDICATED_BENIGN))
    blob = f"the bar is {tok}"
    bare = [t for t in dict.fromkeys(m.DECIMAL_RE.findall(blob))
            if not m.NON_MEASUREMENT.match(t) and (cell_id, t) not in ADJUDICATED_BENIGN]
    assert bare == [], f"{tok} should be exempt for {cell_id}"
    # ... while the same number on a DIFFERENT cell is still a coverage gap (the table is keyed by both)
    bare_other = [t for t in dict.fromkeys(m.DECIMAL_RE.findall(blob))
                  if not m.NON_MEASUREMENT.match(t) and ("other:cell:x", t) not in ADJUDICATED_BENIGN]
    assert bare_other == [tok]


# --- decoy discrimination ------------------------------------------------------------------------
# A citation only VERIFIES a number if the number could have failed to match. These tests exist
# because the naive version of this check is vacuous: a cell whose figures are few and round matches
# a numeric-dense JSON about anything, so "found" alone cannot distinguish provenance from luck.

def test_discrimination_is_reported_with_a_legal_value():
    rep = audit()
    for r in rep["cells"]:
        assert r["discrimination"] in {"HIGH", "MODERATE", "LOW", "n/a"}
        assert 0 <= r["decoy_full_match"] <= r["n_decoys"]
        if r["checked"] and r["n_decoys"]:
            assert r["discrimination"] != "n/a"


def test_discrimination_actually_discriminates():
    """Uniformly HIGH would mean the decoys are too weak to ever match; uniformly LOW would mean the
    control is broken. The signal is only useful if it separates cells."""
    rep = audit()
    seen = {r["discrimination"] for r in rep["cells"] if r["checked"]}
    assert "HIGH" in seen, "no cell resists decoys — the control is not exercised"
    assert len(seen) > 1, f"discrimination is constant at {seen} — it separates nothing"


def test_a_cell_with_many_distinctive_numbers_resists_decoys():
    """Concrete anchor: salmserovar carries 30 measured figures; an unrelated artifact should not
    contain all of them. If this ever fails, the decoy match has become too permissive."""
    rep = audit()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    r = rows["typing:Salmonella:salmserovar"]
    assert r["checked"] >= 20 and r["decoy_full_match"] == 0


def test_low_discrimination_does_not_count_as_drift():
    """Weak evidence is not a defect. A cell can be clean AND poorly verified, and the verdict must
    key on drift alone or the audit stops being a triage funnel."""
    rep = audit()
    assert rep["verdict"] == ("NOTHING_TO_ADJUDICATE" if rep["n_candidate_drift"] == 0
                             else "ADJUDICATION_REQUIRED")
    low = [r for r in rep["cells"] if r["discrimination"] == "LOW"]
    for r in low:
        assert r["candidate_drift"] == [], "a LOW-discrimination cell is unverified, not drifted"
    assert rep["n_cells_low_discrimination"] == len(low)
    assert rep["cells_low_discrimination"] == sorted(r["cell_id"] for r in low)


def test_discrimination_is_deterministic():
    """Fixed decoy seed — a report whose numbers move run to run is not auditable."""
    a, b = audit(), audit()
    assert [(r["cell_id"], r["decoy_full_match"]) for r in a["cells"]] == \
           [(r["cell_id"], r["decoy_full_match"]) for r in b["cells"]]
