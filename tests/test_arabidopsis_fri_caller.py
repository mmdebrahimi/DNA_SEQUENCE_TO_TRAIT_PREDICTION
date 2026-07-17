"""Tests for the FRI variant caller (the flowering cell's v0.1 genome-side input).

Pins the call logic + the abstention discipline + the paper's own arithmetic. Real-data tests skip when the
CC-BY Table S1 or openpyxl is absent.
"""
from __future__ import annotations

import pytest

from dna_decode.organism_rules.arabidopsis_fri_caller import (
    LOF_CONSEQUENCES,
    NON_LOF_CONSEQUENCES,
    VERIFIED_LOF_SUBSTITUTIONS,
    FriCallerError,
    Variant,
    call_fri_from_variants,
    reference_integrity_ok,
)

_STOP = Variant(100, "S316X", "stopgain")
_FRAMESHIFT = Variant(50, "M1fs", "frameshift deletion")
_L294F = Variant(880, "L294F", "nonsynonymous SNV")     # experimentally proved LoF
_L276R = Variant(827, "L276R", "nonsynonymous SNV")     # experimentally proved LoF
_E302G = Variant(905, "E302G", "nonsynonymous SNV")     # proved FUNCTIONAL in isolation
_SYN = Variant(30, "L10L", "synonymous SNV")


# ---- the core call ---------------------------------------------------------------------------------------

def test_alt_at_a_putative_lof_lesion_calls_lof():
    c = call_fri_from_variants({_STOP: "alt", _SYN: "ref"}, "putative")
    assert c.status == "lof" and "S316X" in c.evidence[0]


def test_all_reference_calls_functional():
    assert call_fri_from_variants({_STOP: "ref", _FRAMESHIFT: "ref"}, "putative").status == "functional"


def test_a_synonymous_alt_does_not_call_lof():
    assert call_fri_from_variants({_SYN: "alt", _STOP: "ref"}, "putative").status == "functional"


# ---- abstention: a no-call is not evidence of absence ----------------------------------------------------

def test_nocall_at_a_lof_position_abstains_rather_than_defaulting_functional():
    c = call_fri_from_variants({_STOP: ".", _FRAMESHIFT: "ref"}, "putative")
    assert c.status == "unknown"
    assert c.n_nocall_lof_positions == 1
    assert "not evidence of absence" in c.note


def test_a_positive_lof_hit_wins_over_a_nocall_elsewhere():
    # We can be sure it's LoF even if another position is uncalled -- absence is the uncertain direction.
    c = call_fri_from_variants({_STOP: "alt", _FRAMESHIFT: "."}, "putative")
    assert c.status == "lof" and c.n_nocall_lof_positions == 1


def test_a_nocall_at_a_non_lof_position_is_irrelevant():
    assert call_fri_from_variants({_SYN: ".", _STOP: "ref"}, "putative").status == "functional"


def test_heterozygous_lof_is_called_and_flagged():
    c = call_fri_from_variants({_STOP: "het"}, "putative")
    assert c.status == "lof" and "heterozygous" in c.evidence[0]


# ---- the putative-vs-curated split (the whole point) -----------------------------------------------------

def test_putative_rule_cannot_see_the_verified_substitutions():
    """This IS Table S3's column behavior -- and why its a016/a089 carriers are called functional."""
    assert call_fri_from_variants({_L294F: "alt", _L276R: "alt"}, "putative").status == "functional"


def test_curated_rule_calls_the_verified_substitutions_lof():
    assert call_fri_from_variants({_L294F: "alt"}, "curated").status == "lof"
    assert call_fri_from_variants({_L276R: "alt"}, "curated").status == "lof"


def test_curated_rule_does_not_call_e302g():
    """The paper tested E302G in isolation and it STILL delayed flowering -> functional. Including it would
    be the 'every central-domain change is bad' error the paper explicitly disproves."""
    assert "E302G" not in VERIFIED_LOF_SUBSTITUTIONS
    assert call_fri_from_variants({_E302G: "alt"}, "curated").status == "functional"


def test_curated_and_putative_agree_when_no_verified_substitution_is_present():
    g = {_STOP: "alt", _SYN: "ref"}
    assert call_fri_from_variants(g, "putative").status == call_fri_from_variants(g, "curated").status


# ---- refusals --------------------------------------------------------------------------------------------

def test_unknown_consequence_raises_rather_than_guessing():
    v = Variant(1, "X1Y", "some_new_annovar_term")
    with pytest.raises(FriCallerError, match="refusing to guess"):
        call_fri_from_variants({v: "alt"}, "putative")


def test_unknown_rule_raises():
    with pytest.raises(FriCallerError, match="unknown rule"):
        call_fri_from_variants({_STOP: "alt"}, "made_up")


def test_consequence_vocabularies_are_disjoint():
    assert not (LOF_CONSEQUENCES & NON_LOF_CONSEQUENCES)


def test_integrity_guard_pins_the_verified_set():
    assert reference_integrity_ok() is True


# ---- real data -------------------------------------------------------------------------------------------

def _s1():
    from pathlib import Path
    p = Path(__file__).resolve().parent.parent / "data" / "arabidopsis" / "zhang2020" / "tpj14716-sup-0010-TableS1.xlsx"
    return p if p.exists() else None


@pytest.mark.skipif(_s1() is None, reason="Table S1 not present (browser-only fetch)")
def test_real_table_s1_reconciles_to_the_papers_stated_arithmetic():
    """171 variants; 26 loss-of-function = 4 large indels + 14 frameshift + 7 stopgain + 1 stoploss."""
    pytest.importorskip("openpyxl")
    from dna_decode.organism_rules.arabidopsis_fri_caller import load_variants

    variants, accessions, matrix = load_variants(_s1())
    assert len(variants) == 171
    assert len(accessions) == 1016
    assert sum(1 for v in variants if v.is_putative_lof) == 26
    assert reference_integrity_ok(_s1()) is True


@pytest.mark.skipif(_s1() is None, reason="Table S1 not present (browser-only fetch)")
def test_real_putative_caller_reproduces_table_s3s_own_column():
    """The caller's REAL pass: called from variants, it must agree with the source's own answer."""
    pytest.importorskip("openpyxl")
    import csv
    from pathlib import Path

    from dna_decode.organism_rules.arabidopsis_fri_caller import call_all_accessions

    s3p = Path(__file__).resolve().parent.parent / "data" / "arabidopsis" / "zhang2020" / "tpj14716-sup-0012-TableS3.tsv"
    if not s3p.exists():
        pytest.skip("Table S3 not present")
    s3 = {r["accession_id"]: r for r in csv.DictReader(s3p.open(encoding="utf-8-sig"), delimiter="\t")}
    calls = call_all_accessions(_s1(), "putative")
    dis = [a for a, c in calls.items()
           if a in s3 and c.status != "unknown"
           and c.status != ("lof" if s3[a]["deleterious_allele"] == "TRUE" else "functional")]
    assert dis == [], f"putative caller disagrees with S3's deleterious_allele on {dis}"


@pytest.mark.skipif(_s1() is None, reason="Table S1 not present (browser-only fetch)")
def test_real_curated_rule_flips_exactly_the_a016_and_a089_carriers():
    """The curated rule must be SURGICAL: it may only touch the carriers of the two verified substitutions."""
    pytest.importorskip("openpyxl")
    import csv
    from pathlib import Path

    from dna_decode.organism_rules.arabidopsis_fri_caller import call_all_accessions

    s3p = Path(__file__).resolve().parent.parent / "data" / "arabidopsis" / "zhang2020" / "tpj14716-sup-0012-TableS3.tsv"
    if not s3p.exists():
        pytest.skip("Table S3 not present")
    s3 = {r["accession_id"]: r for r in csv.DictReader(s3p.open(encoding="utf-8-sig"), delimiter="\t")}
    cur = call_all_accessions(_s1(), "curated")
    flipped = {s3[a]["allele_group"] for a, c in cur.items()
               if a in s3 and c.status == "lof" and s3[a]["deleterious_allele"] == "FALSE"}
    assert flipped == {"a016", "a089"}
