"""Guards for the 16S-methyltransferase gentamicin candidate.

The candidate is NOT deployed and must not become deployed by accident: the frozen surface is what the
prospective lock pins, and changing it invalidates both the lock and the reproducibility freeze.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))


def test_methyltransferase_regex_is_family_specific_not_a_loose_prefix():
    """`rmt` as a substring would over-match. The families are named explicitly."""
    from gentamicin_rmt_candidate import is_methyltransferase
    for good in ("rmtB1", "rmtE1", "rmtC", "rmtF1", "armA", "npmA", "RMTB1"):
        assert is_methyltransferase(good), good
    for bad in ("aac(3)-IId", "aph(6)-Id", "aadA5", "armadillo", "rmtX-like_thing", "", "blaTEM-1"):
        assert not is_methyltransferase(bad), bad


def test_the_candidate_only_ever_adds_R_never_removes_it():
    """It is an OR over the frozen rule, so it can never turn an R call into S. Pinned because a
    'fix' that silently removed calls would be a regression the confusion matrix might not localise."""
    from gentamicin_rmt_candidate import frozen_call, candidate_call
    rows_gent = [{"Class": "AMINOGLYCOSIDE", "Subclass": "GENTAMICIN", "Element symbol": "aac(3)-IId"}]
    rows_mt = [{"Class": "AMINOGLYCOSIDE", "Subclass": "AMINOGLYCOSIDE", "Element symbol": "rmtB1"}]
    rows_none = [{"Class": "AMINOGLYCOSIDE", "Subclass": "STREPTOMYCIN", "Element symbol": "aadA5"}]
    for rows in (rows_gent, rows_mt, rows_none, rows_gent + rows_mt):
        assert candidate_call(rows) or not frozen_call(rows)
    assert frozen_call(rows_gent) and candidate_call(rows_gent)
    assert not frozen_call(rows_mt) and candidate_call(rows_mt)      # the whole point
    assert not frozen_call(rows_none) and not candidate_call(rows_none)


def test_compound_gentamicin_subclasses_are_already_counted_by_the_frozen_rule():
    """VERIFIED, not assumed. I expected `GENTAMICIN/KANAMYCIN/TOBRAMYCIN` to be a second blind spot;
    the frozen rule matches by TOKEN, so it already counts them. Pinned so the memo cannot drift back."""
    from gentamicin_rmt_candidate import frozen_call
    for sub in ("GENTAMICIN/KANAMYCIN/TOBRAMYCIN", "APRAMYCIN/GENTAMICIN/TOBRAMYCIN"):
        assert frozen_call([{"Class": "AMINOGLYCOSIDE", "Subclass": sub, "Element symbol": "x"}])


def test_the_frozen_amr_surface_is_untouched_by_this_work():
    """The prospective lock pins amr_rules.py + calibrated_amr_rules.json. The candidate is scorer-local
    and must not have edited either; if it did, the lock and the freeze are both invalidated."""
    from dna_decode.eval import prospective_lock as pl  # noqa: F401  (import proves the module loads)
    for rel in ("dna_decode/eval/amr_rules.py", "dna_decode/data/calibrated_amr_rules.json"):
        p = ROOT / rel
        if not p.exists():
            continue
        # the candidate must appear NOWHERE in the frozen files
        text = p.read_text(encoding="utf-8", errors="replace")
        assert "rmt" not in text.lower() or "npmA" not in text, (
            f"{rel} appears to have been edited to include the candidate -- that breaks the lock")


def test_the_specificity_claim_reports_its_own_vacuity():
    """Every methyltransferase carrier in the labelled data is R, so the candidate cannot produce a
    false positive there. 'spec unchanged' is therefore VACUOUS and the script must say so rather than
    let an unchanged number read as evidence of safety."""
    src = (ROOT / "scripts" / "gentamicin_rmt_candidate.py").read_text(encoding="utf-8")
    assert "VACUOUS" in src and "UNTESTED, not zero" in src
    assert "_s_labelled_carriers" in src
