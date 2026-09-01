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


def test_the_candidate_is_now_deployed_under_its_own_lock():
    """SUPERSEDED 2026-08-31, and the supersession is the point.

    This test used to assert the candidate appeared NOWHERE in the frozen surface -- correct while the
    candidate was scorer-local. The user then authorized the v2 lock, so the candidate IS the deployed
    rule and the old assertion is simply false. The honest replacement is not deletion: it is to guard
    the NEW invariant, which is that a surface change of this kind never happens WITHOUT its own lock.

    So: if amr_rules.py carries the rescue, a v2 manifest must exist, pin the live file, and carry a
    later lock_date than v1. That makes an undocumented edit to the frozen surface fail here.
    """
    import hashlib
    import json

    from dna_decode.eval.amr_rules import rule_for

    text = (ROOT / "dna_decode/eval/amr_rules.py").read_text(encoding="utf-8", errors="replace")
    deployed = bool(rule_for("gentamicin").get("symbol_rescue"))
    if not deployed:
        assert "npmA" not in text, "the rescue is in the surface but not in the rule -- half-applied edit"
        return

    manifests = sorted((ROOT / "wiki").glob("prospective_lock_manifest_*.json"))
    assert len(manifests) >= 2, "the rescue is deployed but no NEW lock manifest was minted"
    v2 = json.loads(manifests[-1].read_text(encoding="utf-8"))
    live = hashlib.sha256((ROOT / "dna_decode/eval/amr_rules.py").read_bytes()).hexdigest()
    assert v2["surface_sha256"]["dna_decode/eval/amr_rules.py"] == live, (
        "the deployed surface does not match the newest lock manifest -- the lock is stale")
    v1 = json.loads(manifests[0].read_text(encoding="utf-8"))
    assert v2["lock_date"] > v1["lock_date"], "a superseding lock must carry a later cutoff"
    assert v2.get("supersedes"), "a superseding lock must name what it retires"


def test_the_retired_v1_lock_is_preserved_not_deleted():
    """The v1 manifest is the only record of what was locked 2026-06-13..2026-08-31, and the v1
    prospective scores are interpretable ONLY against it. Deleting it would erase the audit trail."""
    assert (ROOT / "wiki" / "prospective_lock_manifest_2026-06-22.json").exists()


def test_the_specificity_claim_reports_its_own_vacuity():
    """Every methyltransferase carrier in the labelled data is R, so the candidate cannot produce a
    false positive there. 'spec unchanged' is therefore VACUOUS and the script must say so rather than
    let an unchanged number read as evidence of safety."""
    src = (ROOT / "scripts" / "gentamicin_rmt_candidate.py").read_text(encoding="utf-8")
    assert "VACUOUS" in src and "UNTESTED, not zero" in src
    assert "_s_labelled_carriers" in src


def test_the_gap_is_rmt_only_armA_is_already_counted():
    """MEASURED, not assumed: across 158 cached methyltransferase rows the AMRFinder Subclass is
    perfectly consistent per GENE -- armA -> GENTAMICIN 24/24 (already counted by the frozen rule),
    rmt* -> AMINOGLYCOSIDE 134/134 (invisible). So the rule gap is rmt-only and an armA clause is a
    no-op. The accrual memo phrasing "rmtE1/rmtE/armA-family" lumps in a gene that was never missing.
    """
    from gentamicin_rmt_candidate import in_rmt_gap, is_methyltransferase
    assert is_methyltransferase("armA") and not in_rmt_gap("armA")
    for g in ("rmtB1", "rmtE1", "rmtF1", "npmA"):
        assert in_rmt_gap(g) and is_methyltransferase(g)


def test_narrowing_to_the_gap_changes_no_score():
    """If dropping armA from the RULE changed a call, armA would not have been a no-op and the per-gene
    subclass measurement would be wrong."""
    from gentamicin_rmt_candidate import frozen_call, candidate_call
    arma = [{"Class": "AMINOGLYCOSIDE", "Subclass": "GENTAMICIN", "Element symbol": "armA"}]
    assert frozen_call(arma) and candidate_call(arma)


def test_label_hunt_draws_no_conclusion_from_a_failed_sweep():
    """A partial sweep zero means the sweep failed, not that none exist. The first version printed a
    confident "still ZERO S-labelled carriers" immediately after printing INCOMPLETE."""
    src = (ROOT / "scripts" / "gentamicin_rmt_label_hunt.py").read_text(encoding="utf-8")
    assert "NO CONCLUSION" in src and "COMPLETE sweep" in src
    assert "parse_ast_phenotypes" in src, "must reuse the fixed PD parser, not re-split the field"


def test_the_disjoint_pool_is_manifest_gated_not_selected_tsv_gated():
    """The cheap check under-covered by two thirds: 956 -> 311 once the accession manifest ran.

    Pinned because the manifest is fail-closed and the hand-rolled version is not, and because an
    unasserted edit once left the gate un-applied while the numbers looked plausible.
    """
    src = (ROOT / "scripts" / "unscored_genome_label_census.py").read_text(encoding="utf-8")
    assert "prior_accessions" in src and "INCOMPLETE_MANIFEST" in src
    assert "DISJOINT candidate pool" in src


def test_disjoint_validation_numbers_are_recorded_with_their_vacuity_caveat():
    """The rescue is measured (+0.369 sens on 131 disjoint isolates) but the specificity result is still
    vacuous on rmt -- zero S-labelled carriers in any of three datasets. Both must be stated together;
    reporting the gain without the caveat would overclaim safety."""
    memo = ROOT / "wiki" / "gentamicin_rmt_disjoint_validation_2026-08-28.md"
    if not memo.exists():
        import pytest
        pytest.skip("memo absent")
    text = memo.read_text(encoding="utf-8")
    assert "0.892" in text and "0.523" in text
    assert "arithmetic, not evidence" in text
    assert "an absence is not a bound" in text
