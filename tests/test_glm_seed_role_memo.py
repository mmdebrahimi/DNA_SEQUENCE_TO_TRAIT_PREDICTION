"""The memo's numbers must be the ones the committed artifacts hold, and its retractions must survive edits.

Two jobs. First, the prose cannot drift from the computation: every headline figure is re-read from the
artifact and compared. Second -- and this is the part a later edit is most likely to undo -- the memo must
keep RETRACTING the figures it corrects. A tidy-up that deletes the withdrawn 0.0543 band would leave the
record looking like the arm always knew the kmer margin was negative, and the whole value of the correction
is that it was not.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MEMO = ROOT / "wiki" / "glm_seed_role_decomposition_2026-11-21.md"
CEILINGS = ROOT / "wiki" / "glm_seed_role_ceilings_2026-11-21.json"
VARIANCE = ROOT / "wiki" / "glm_seed_role_variance_2026-11-21.json"

pytestmark = pytest.mark.skipif(
    not MEMO.exists() or not CEILINGS.exists() or not VARIANCE.exists(),
    reason="the seed-role arm has not been run on this host")


@pytest.fixture(scope="module")
def memo():
    return MEMO.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def flat(memo):
    """Whitespace-collapsed memo. Multi-word phrases must not break on a line wrap."""
    return re.sub(r"\s+", " ", memo)


@pytest.fixture(scope="module")
def ceilings():
    return json.loads(CEILINGS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def variance():
    return json.loads(VARIANCE.read_text(encoding="utf-8"))


# ------------------------------------------------------------------ ordering (reader-action first)

def test_decisions_needed_comes_before_corrections_and_before_the_result(memo):
    """Required reader action first. A buried decision is the documented failure of this memo's shape."""
    d = memo.index("## DECISIONS NEEDED")
    c = memo.index("## CORRECTIONS")
    r = memo.index("## The result")
    assert d < c < r, (d, c, r)


# ------------------------------------------------------------------ the numbers are the artifact's

def test_the_per_mode_margins_match_the_committed_artifact(memo, ceilings):
    for mode in ("kmer", "both"):
        s = ceilings["summary"][mode]
        for value in (f"{s['null_mean']:.4f}", f"{s['ceiling_mean']:.4f}"):
            assert value in memo, f"{mode} {value} is in the artifact and not in the memo"
    assert "-0.0121" in memo or "−0.0121" in memo, "the kmer margin is missing"
    assert "+0.0964" in memo or "0.0964" in memo, "the both-mode margin is missing"


def test_the_mode_verdicts_match_the_artifact_and_differ_from_each_other(memo, ceilings):
    v = ceilings["verdict_by_mode"]
    assert v["kmer"] in memo and v["both"] in memo
    assert v["kmer"] != v["both"], "the per-mode split is the memo's central finding"


def test_the_variance_arms_and_the_additivity_ratio_match(memo, variance):
    d = variance["decomposition"]["both"]
    for value in (d["sd_all_coupled"], d["sd_generation_only"], d["sd_corpus_only"]):
        assert f"{value:.4f}" in memo, f"{value} is in the artifact and not in the memo"
    ratio = variance["mechanism_variance_additivity"]["both"]["sum_over_coupled_ratio"]
    assert f"{ratio}" in memo, f"the additivity ratio {ratio} is the mechanism and must be stated"


def test_the_anchor_sweep_refutation_is_reported_with_its_range(memo, variance):
    sens = variance["mechanism_anchor_sensitivity"]
    assert sens["rival_status"].startswith("REFUTED")
    assert f"{sens['sd_min']:.4f}" in memo and f"{sens['sd_max']:.4f}" in memo
    assert "refuted" in memo.lower()


# ------------------------------------------------------------------ the retractions must persist

def test_the_withdrawn_kmer_band_is_still_named_as_withdrawn(flat):
    """0.0543 must appear AS AN ERROR. Deleting it would hide that the arm once believed it."""
    assert "0.0543" in flat, "the withdrawn band was deleted instead of retracted"
    assert "0.5990" in flat, "the retraction must name the maximum-of-ten null that produced it"
    assert "wrong twice over" in flat


def test_the_zeroshot_qualifier_on_the_gap_is_stated(flat):
    assert "0.1185" in flat and "0.1355" in flat, (
        "both the zeroshot gap and the all-model mean must appear, or the qualifier is not a correction")


def test_the_plans_own_hypothesis_is_recorded_as_falsified(flat):
    assert "measured false" in flat.lower()
    assert "SUPPRESSES" in flat or "suppresses" in flat


def test_seed_7_is_named_as_an_outlier_in_BOTH_modes(flat):
    assert "0.5865" in flat and "0.5419" in flat
    assert "unrepresentative in both modes" in flat


# ------------------------------------------------------------------ claims it must NOT make

def test_the_memo_does_NOT_claim_monotone_finetuning_degradation(flat):
    """The epoch series looks monotone and is the SAME delta that returned DELTA_NOT_ESTABLISHED."""
    assert "DELTA_NOT_ESTABLISHED" in flat
    assert "0.2072" in flat, "the paired spread that makes the series unpowered must be named"
    assert "not claiming it" in flat.lower()


def test_the_memo_scopes_the_variance_result_to_the_markov_probe(flat):
    assert "Markov null" in flat
    assert "says nothing" in flat.lower() or "different quantity" in flat.lower()


def test_the_memo_keeps_the_necessary_not_sufficient_rail(flat):
    """Clearing a local-composition discriminator is not functional realism."""
    assert "necessary, not sufficient" in flat.lower()
    assert "LOCAL" in flat or "local composition" in flat.lower()


def test_the_memo_does_not_self_execute_its_own_scope_recommendation(flat):
    """Striking plan steps and retiring a scoring mode are authority calls, not executor actions."""
    assert "not" in flat and "self-executed" in flat
    assert "I have **not** edited the plan" in flat or "have not edited the plan" in flat


# ------------------------------------------------------------------ the supersession (added 2026-11-21)

def test_the_memo_is_marked_superseded_in_its_first_lines(flat, memo):
    """The whole file is withdrawn, and the banner must be the FIRST thing a reader meets.

    Position is load-bearing, not cosmetic: `doc_validity` scans the first 30 lines for a status header,
    `continuity` injection is cut at 4,000 characters, and compaction retains a measured ~17% of injected
    constraints. A correction below the fold is a correction nobody receives. The body is deliberately
    left unedited -- the sequence of corrections is itself the finding -- so the banner is the only thing
    standing between a reader and a withdrawn number.
    """
    first = "\n".join(memo.splitlines()[:30])
    assert "SUPERSEDED 2026-11-21" in first
    assert "do not quote a number from this file" in first
    assert "glm_instrument_reconciled_2026-11-21.md" in first


def test_the_banner_names_the_geometry_defect_numerically(flat):
    """Prose alone would not survive a tidy-up; the numbers are what make the record checkable."""
    assert "n_generated=400" in flat and "n=120" in flat
    assert "glm_paired_sweep_kernel.py:45" in flat


def test_the_banner_withdraws_all_three_headline_claims(flat):
    """Each of the three must be named as withdrawn, with the figure that replaces it."""
    assert "+0.0329" in flat, "the corrected kmer margin must appear"
    assert "1.129" in flat, "the fit-decoupled arm's refuting ratio must appear"
    assert "0.0409" in flat, "the sampling benchmark that replaces the decoupling claim must appear"
    # Matched on a phrase that sits within ONE line of the banner. The `flat` fixture collapses
    # whitespace, which folds the markdown blockquote prefix of the next line into the text, so a
    # phrase spanning a line break becomes "... has nothing > to remove" rather than the plain string.
    assert "seed-role split has nothing" in flat


def test_the_banner_says_what_SURVIVES_not_only_what_fell(flat):
    """A supersession that withdraws everything indiscriminately destroys a valid measurement.

    The per-seed ceilings are untouched by the geometry error, because `self_control_ceiling` reads only
    the eval half, which depends on limit and seed and NOT on n_generated.
    """
    assert "What survives" in flat
    assert "ceilings" in flat
