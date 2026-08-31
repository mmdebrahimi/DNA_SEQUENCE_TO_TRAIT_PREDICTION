"""Pin the L2 doubt layer's arithmetic and its tiering.

The measurement that motivated the surprise test is real: the full-index screen flags gentamicin
`rmtE1` at 36R/0S (the known gap) AND ciprofloxacin `qnrA1` at 4R/0S in the same table. Only one of
those means anything, and the raw signature cannot tell them apart.

Offline, pure, no fixtures.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.doubt import (NONE, STRONG, WEAK, binom_lower_tail,  # noqa: E402
                                   completeness_signal, completeness_surprise, completeness_tier,
                                   position_novelty_signal)


# --- the arithmetic ---

def test_binom_lower_tail_matches_hand_computed_values():
    assert binom_lower_tail(0, 1, 0.5) == pytest.approx(0.5)
    assert binom_lower_tail(0, 4, 0.5) == pytest.approx(0.0625)
    assert binom_lower_tail(2, 2, 0.3) == pytest.approx(1.0)      # whole distribution
    assert binom_lower_tail(0, 3, 0.0) == pytest.approx(1.0)      # p=0 -> zero successes is certain


def test_zero_trials_is_certain_not_significant():
    """A family with no labelled carriers must never look surprising."""
    assert binom_lower_tail(0, 0, 0.5) == 1.0


def test_surprise_is_none_when_there_are_no_labelled_carriers():
    """Absence of labels is UNASSESSABLE, not a clean bill of health."""
    assert completeness_surprise(0, 0, 0.5) is None
    assert completeness_tier(None, 100) == NONE


# --- the tiering, on the real numbers ---

def test_the_known_gentamicin_gap_survives_the_familywise_correction():
    """rmtE1: 36 labelled carriers, all resistant, against a cohort base S-rate of 0.517."""
    p = completeness_surprise(36, 0, 0.517)
    assert p < 1e-10, p
    assert completeness_tier(p, 131) == STRONG


def test_a_tiny_pure_family_does_NOT_survive_the_correction():
    """The defect the full-index run exposed: cipro `qnrA1` is 4R/0S -- pure, and meaningless.

    Base S-rate 0.583, so P(all four resistant) = 0.030 by chance, against ~125 uncounted families
    screened for that drug. Reporting it beside rmtE1 as the same kind of finding is what would make
    the doubt layer noise rather than signal.
    """
    p = completeness_surprise(4, 0, 0.583)
    assert 0.01 < p < 0.05, p          # nominally significant...
    assert completeness_tier(p, 125) == WEAK    # ...and still not strong


def test_an_unsurprising_family_is_tiered_none():
    """Ceftriaxone `blaTEM` at 9R/0S with base S-rate 0.202 -- not even nominally significant."""
    p = completeness_surprise(9, 0, 0.202)
    assert p > 0.05, p
    assert completeness_tier(p, 216) == NONE


def test_a_mixed_family_is_never_strong():
    """One susceptible carrier ENDS the signal -- it is evidence the exclusion is deliberate."""
    assert completeness_surprise(62, 28, 0.517) is None
    assert completeness_tier(completeness_surprise(62, 28, 0.517), 131) == NONE
    assert completeness_surprise(40, 1, 0.517) is None       # even one is enough


def test_an_enrichment_null_WOULD_have_called_a_correct_exclusion_strong():
    """Why the purity null, pinned. This is the defect the full-index run exposed.

    Gentamicin `aph(6)-Id` is 62R/28S. Under a lower-tail binomial on the observed susceptible count
    it is p ~ 5e-5 -- STRONG even after the family-wise correction -- because 28 susceptible carriers
    really are fewer than a 0.517 base rate predicts. But `aph(6)-Id` is a CORRECT exclusion: a
    streptomycin determinant that travels with gentamicin resistance by linkage. Every co-occurring
    determinant is R-enriched, so an enrichment null would flood the layer with correct exclusions.
    """
    enrichment_p = binom_lower_tail(28, 90, 0.517)
    assert enrichment_p < 0.05 / 131, enrichment_p          # the wrong null clears the bar...
    assert completeness_tier(completeness_surprise(62, 28, 0.517), 131) == NONE   # ...the right one does not


def test_the_reason_distinguishes_unassessable_from_deliberate_exclusion():
    """Both tier NONE, for opposite reasons -- collapsing them would hide which is which."""
    assert "unassessable" in completeness_signal("x", "C", 0, 0, 0.5, 50).reason
    assert "deliberate" in completeness_signal("aph(6)-Id", "STREPTOMYCIN", 62, 28, 0.517, 131).reason


def test_correction_scales_with_the_number_of_families_screened():
    """The same p-value must not be STRONG for a broad screen and a narrow one alike."""
    p = completeness_surprise(8, 0, 0.517)
    assert completeness_tier(p, 1) == STRONG          # a single pre-specified family
    assert completeness_tier(p, 500) == WEAK          # one of five hundred


# --- signals ---

def test_completeness_signal_reason_says_unassessable_when_unlabelled():
    sig = completeness_signal("x", "CLASS", 0, 0, 0.5, 50)
    assert sig.tier == NONE
    assert "unassessable" in sig.reason
    assert sig.evidence["purity_surprise_p"] is None


def test_completeness_signal_carries_its_own_correction_context():
    """A tier is unreadable without the correction it survived; both must ship with it."""
    sig = completeness_signal("rmtE1", "AMINOGLYCOSIDE", 36, 0, 0.517, 131)
    assert sig.tier == STRONG
    assert sig.evidence["n_families_tested"] == 131
    assert sig.evidence["familywise_alpha"] == 0.05


def test_position_novelty_signal_delegates_to_the_shipped_flag():
    """A novel substitution AT a catalogued HIV NNRTI position fires; a catalogued one does not."""
    from dna_decode.data.hiv_amr import NNRTI_RT_MAJOR_DRMS
    known = sorted(NNRTI_RT_MAJOR_DRMS)[0]
    wt, pos, _mut = known[0], known[1:-1], known[-1]

    novel = f"{wt}{pos}W" if not known.endswith("W") else f"{wt}{pos}Y"
    fired = position_novelty_signal([novel], "hiv-nnrti-rt")
    assert fired.tier == WEAK
    assert novel in fired.evidence["novel_substitutions"]

    quiet = position_novelty_signal([known], "hiv-nnrti-rt")
    assert quiet.tier == NONE
    assert quiet.evidence["novel_substitutions"] == []


def test_position_novelty_signal_is_quiet_on_an_empty_genotype():
    sig = position_novelty_signal([], "hiv-nnrti-rt")
    assert sig.tier == NONE
    assert sig.evidence["n_catalog_positions"] > 0, "precondition: the cell has a real catalog"


# --- AMR-arm firing rate: predeclared categories, never called a false-positive rate ---

def _spec(rows):
    import sys as _s
    _s.path.insert(0, str(ROOT / "scripts"))
    from doubt_layer_per_cell import amr_arm_specificity
    return amr_arm_specificity(rows)


def test_a_drug_with_no_confirmed_gap_is_unconfirmed_never_clean():
    """Both known gaps were invisible until independent labels arrived, so absence of a confirmed gap
    is absence of EVIDENCE. Calling it clean would assert the negative this project has twice been
    wrong about."""
    s = _spec([{"drug": "tetracycline", "status": "scored", "n_families_uncounted": 89,
                "n_strong": 0, "n_raw_signature": 0}])
    assert s["predeclared_categories"]["unconfirmed"] == ["tetracycline"]
    assert "clean" not in str(s["predeclared_categories"]).lower()
    assert "NOT a false-positive rate" in s["interpretation"]
    assert "ambiguous" in s["interpretation"].lower()


def test_the_confirmed_gap_drug_is_scored_separately_from_the_rest():
    """Pooling the drug that HAS the gap with the drugs that do not would hide both numbers."""
    s = _spec([
        {"drug": "gentamicin", "status": "scored", "n_families_uncounted": 131, "n_strong": 1,
         "n_raw_signature": 1},
        {"drug": "ciprofloxacin", "status": "scored", "n_families_uncounted": 125, "n_strong": 0,
         "n_raw_signature": 2},
    ])
    assert s["confirmed_gap_drugs"] == {"n_drugs": 1, "n_families_screened": 131, "n_strong": 1}
    assert s["unconfirmed_drugs"]["n_strong_after_correction"] == 0
    assert s["unconfirmed_drugs"]["correction_removed"] == 2


def test_an_unlabelled_drug_is_unassessable_not_counted_as_a_pass():
    """Oxacillin has no labels; folding it into either bucket would fabricate evidence."""
    s = _spec([{"drug": "oxacillin", "status": "no_labels", "n_families_uncounted": 401}])
    assert s["predeclared_categories"]["unassessable_no_labels"] == ["oxacillin"]
    assert s["unconfirmed_drugs"]["n_families_screened"] == 0


def test_the_real_run_fires_once_across_every_screened_family():
    """The measured result, pinned. If the screen ever fires on an unconfirmed drug that is a real
    finding needing adjudication -- it must fail loudly, not pass unnoticed."""
    import json
    hits = sorted(ROOT.glob("wiki/doubt_layer_per_cell_*.json"))
    if not hits:
        pytest.skip("per-cell artifact not generated on this checkout")
    s = json.loads(hits[-1].read_text(encoding="utf-8")).get("amr_arm_specificity")
    if not s:
        pytest.skip("artifact predates the specificity block")
    assert s["confirmed_gap_drugs"]["n_strong"] == 1
    assert s["unconfirmed_drugs"]["n_strong_after_correction"] == 0, (
        "the screen fired on a drug with no confirmed gap -- adjudicate it, do not silence it")
    assert s["unconfirmed_drugs"]["correction_removed"] >= 1, "the correction must be doing work"
