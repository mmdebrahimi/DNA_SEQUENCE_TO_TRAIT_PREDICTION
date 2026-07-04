"""Offline pins for the fail-closed LD-imputation pre-processor (dna_decode/imputation.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.abo_blood import call_abo_o_status  # noqa: E402
from dna_decode.imputation import LdImputer, call_with_imputation  # noqa: E402

# synthetic frozen map: one high-purity tag genotype (GG->DD@0.97) + one low-purity (AG->DI@0.60)
_IMP = LdImputer(target="rs8176719", tag="rs657152",
                 table={"GG": {"majority": "DD", "purity": 0.97, "n": 100},
                        "AG": {"majority": "DI", "purity": 0.60, "n": 50},
                        "AA": {"majority": "II", "purity": 0.95, "n": 40}},
                 min_purity=0.90)


def test_direct_passthrough_when_target_called():
    imp = _IMP.impute("DD", "GG")
    assert imp.genotype == "DD" and imp.provenance == "direct" and imp.confidence == 1.0
    imp2 = _IMP.impute("II", "")                     # target present, tag missing -> still direct
    assert imp2.genotype == "II" and imp2.provenance == "direct"


def test_imputes_when_target_uncallable_and_tag_high_purity():
    imp = _IMP.impute("--", "GG")
    assert imp.genotype == "DD" and imp.provenance.startswith("imputed:rs657152=GG@") and imp.confidence == 0.97


def test_fail_closed_low_purity():
    imp = _IMP.impute("--", "AG")                    # AG purity 0.60 < 0.90 -> ABSTAIN, not a guess
    assert imp.genotype is None and imp.provenance.startswith("abstain:low-purity")


def test_fail_closed_no_tag_and_unseen_tag():
    assert _IMP.impute("--", "--").provenance == "abstain:no-tag"
    assert _IMP.impute("--", "CC").provenance.startswith("abstain:tag-genotype-unseen")


def test_call_with_imputation_end_to_end():
    # uncallable O-deletion + high-purity tag -> imputed DD -> deterministic rule calls O
    r = call_with_imputation(call_abo_o_status, "--", "GG", _IMP)
    assert r["call"] == "O" and r["provenance"].startswith("imputed:") and r["confidence"] == 0.97
    # directly-typed non-O passes through
    d = call_with_imputation(call_abo_o_status, "DI", "AG", _IMP)
    assert d["call"] == "non-O" and d["provenance"] == "direct"
    # no confident tag -> the decoder still ABSTAINS honestly
    a = call_with_imputation(call_abo_o_status, "--", "AG", _IMP)
    assert a["call"] == "INDETERMINATE" and a["provenance"].startswith("abstain:")


def test_committed_map_loads_if_present():
    imp = LdImputer.for_target("rs8176719")
    if imp is None:
        import pytest
        pytest.skip("committed ABO map not present yet")
    assert imp.target == "rs8176719" and imp.tag == "rs657152" and len(imp.table) >= 1
