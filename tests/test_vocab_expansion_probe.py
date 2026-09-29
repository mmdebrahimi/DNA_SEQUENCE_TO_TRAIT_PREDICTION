"""Guards for the vocabulary-expansion probe — the measurement that STOPPED an expansion.

The probe exists because one entry (`freq`) worked and the conclusion drawn was "vocabulary coverage has
leverage; bindings do not". That generalised from n=1. Measured across all 141 remaining unlabeled numbers
the lever is exhausted, and these tests pin the two ways the measurement could have lied in the flattering
direction: a statistical modifier reading as a clean conversion, and a net-negative candidate reading as
harmless.
"""

from __future__ import annotations

from scripts.vocab_expansion_probe import MODIFIERS, probe, verdict_for


def test_a_statistical_modifier_is_REFUSED_however_good_its_numbers_look():
    """THE dangerous class. `max`/`mean`/`null` sit near numbers constantly and often appear in the matching
    path, so on raw counts they look like the best candidates available. They are the worst: prose "max 0.21"
    against a path `max_abs_r` would CONFIRM the quantity as "max" while the number is a CORRELATION. A
    confirmation that does not identify the quantity is vacuous, which is worse than no confirmation."""
    for w in ("max", "mean", "null", "delta", "gain"):
        assert w in MODIFIERS, w
        v = verdict_for(w, {"near": 9, "would_confirm": 9, "would_mismatch": 0})
        assert v.startswith("REFUSE"), (w, v)
        assert "vacuous" in v, v


def test_a_net_negative_candidate_does_NOT_read_as_harmless():
    """MY OWN BUG, pinned. The first version's final `else` printed 'no effect (no artifact hit)' for the
    confirm-but-fewer-than-mismatch case, so `purity` at 1 conversion against 4 FALSE LEADS read as benign.
    A verdict label that understates a net-negative candidate is exactly how a bad entry gets added."""
    v = verdict_for("purity", {"near": 5, "would_confirm": 1, "would_mismatch": 4})
    assert v.startswith("REFUSE"), v
    assert "net negative" in v and "4" in v, v
    # and the genuinely inert case must still be distinguishable from it
    inert = verdict_for("somenoun", {"near": 3, "would_confirm": 0, "would_mismatch": 0})
    assert "no effect" in inert, inert


def test_false_leads_are_treated_as_a_COST_not_a_neutral_outcome():
    """A false lead is strictly worse than an unlabeled number: the lead costs adjudication time while an
    unlabeled number honestly reports 'not measured'. So zero conversions with any false leads must refuse."""
    v = verdict_for("somenoun", {"near": 4, "would_confirm": 0, "would_mismatch": 3})
    assert v.startswith("REFUSE") and "0 conversions" in v, v


def test_the_live_corpus_yields_no_QUANTITY_worth_adding():
    """THE FINDING -- with the correction that made it precise.

    First measured as "zero addable entries". That held only at the coverage OF THAT MOMENT: citing
    pneumoserotype's two artifacts moved its 5 numbers into the cited set, and `concordance` immediately
    surfaced as a clean candidate (2 near / 2 would-confirm / 0 would-mismatch) and was added. So the
    vocabulary lever's headroom is **a function of which cells are cited** and RE-OPENS whenever coverage
    grows -- it is not a fixed quantity that can be exhausted once. This test is the tripwire for that.
     Every word the probe calls a clean 'candidate' on the
    real corpus is a DOMAIN NOUN or a structural word (genes / lactam / beta / aminoglycoside / star / name /
    computed / discovery), not a quantity — adding those would inflate the measurability figure while
    measuring nothing. Pinned loosely (a membership check, not an exact list) so the test survives ordinary
    corpus drift but fails if a genuine new quantity ever appears and goes unnoticed."""
    rep = probe()
    assert "error" not in rep, rep
    assert rep["n_unlabeled"] > 100, rep["n_unlabeled"]
    clean = [w for w, st in rep["words"].items()
             if verdict_for(w, st).startswith("candidate")]
    assert clean, "expected some cleanly-converting words, else the probe measures nothing"
    # None of them is a quantity. These are the observed ones; the assertion is that the clean set stays
    # inside the non-quantity vocabulary rather than that it equals this exact list.
    non_quantities = {"genes", "lactam", "beta", "aminoglycoside", "star", "name", "height", "learned",
                      "computed", "discovery", "axis", "ours", "asserted", "serotype", "coat", "dog",
                      "morphology", "owner", "reported", "panel", "arm", "cohort", "isolates", "strains",
                      # population codes and database names -- surfaced by this very test
                      "eas", "eur", "afr", "ncbi", "getrm", "clsi", "eucast",
                      # lab METHOD names and cohort-slice words -- surfaced when pneumoserotype gained
                      # citations. `quellung` is a serotyping method, `subset` names a slice of a cohort;
                      # neither is a quantity. The same event surfaced `concordance`, which IS one, and it
                      # was added to LABEL_VOCAB -- so it no longer appears here.
                      "quellung", "subset", "method", "slice", "explicit"}
    surprises = [w for w in clean if w not in non_quantities]
    assert not surprises, (
        f"a cleanly-converting word that is NOT in the known non-quantity set appeared: {surprises}. "
        "If any of these IS a genuine quantity, that is a real vocabulary candidate — measure its "
        "would_mismatch count before adding it.")
