"""The ported SeqSero2 O decision procedure, tested as PURE logic -- no blastn, no DB, no Docker.

`call_o_antigen` consumes the shared engine's `per_allele` dict, so every branch can be exercised by
handing it a synthetic dict. That matters because the branches encode facts that are easy to state
backwards: the two `not_in` markers are used in OPPOSITE directions, and the marker that reads like it
should assert `1,3,19` can never produce a positive call at all.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.runner import (  # noqa: E402
    NOT_IN_310_KEY,
    NOT_IN_1319_KEY,
    WBAV_KEYS,
    WZX_310_KEY,
    WZY_946_FULL_KEY,
    call_o_antigen,
    parse_axis_antigen,
    parse_ss2_key,
)

_IDX = iter(range(1000, 9999))


def hit(ss2_key: str, cov: float, *, axis: str = "O", antigen: str | None = None,
        called: bool = True, identity: float = 99.0) -> tuple[str, dict]:
    """Build one `per_allele` entry in the real header schema `<axis>__<antigen>__<idx>__<ss2 key>`."""
    if antigen is None:
        antigen = ss2_key.split("_")[0][2:] if axis == "O" else ss2_key
    allele_id = f"{axis}__{antigen}__{next(_IDX)}__{ss2_key}"
    return allele_id, {"called": called, "percent_identity": identity, "percent_coverage": cov}


def per_allele(*entries) -> dict:
    return dict(entries)


# --- the header round-trip the whole procedure rests on --------------------------------------------

def test_the_role_field_round_trips_even_though_it_embeds_the_delimiter():
    """`O-9,46_wbaV__1002` contains `__`, so a naive `parts[3]` would truncate it to `O-9,46_wbaV` and
    every literal-key comparison in the procedure would silently miss."""
    aid = "O__9,46__239__O-9,46_wbaV__1002"
    assert parse_ss2_key(aid) == "O-9,46_wbaV__1002"
    assert parse_axis_antigen(aid) == ("O", "9,46")


def test_a_specific_gene_can_never_be_selected_as_an_antigen():
    """`SP__` entries exist only to feed the O9-vs-O2 refinement. If `parse_axis_antigen` admitted
    them, a tyr gene could win `_best_per_axis` and be emitted as an O antigen."""
    assert parse_axis_antigen("SP__tyr-O-9__325__tyr-O-9__609") is None


# --- Branch A: the O9 family ------------------------------------------------------------------------

def test_wbaV_with_wzy_gives_9_46():
    got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 95.0), hit(WZY_946_FULL_KEY, 90.0)))
    assert got == {"antigen": "9,46", "o_antigen_rule": "A_wbaV_with_wzy"}


def test_wbaV_alone_gives_plain_O9_which_the_old_caller_could_not_emit():
    """The 15 measured `9,46`-where-the-reference-says-`9` disagreements. Plain O9 has no allele in
    EITHER database -- it is reached as the DEFAULT of the wbaV branch."""
    got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 95.0)))
    assert got["antigen"] == "9"
    assert got["o_antigen_rule"] == "A_wbaV_default_O9_tyr_not_discriminable_by_blast"


def test_the_wbaV_gate_is_a_threshold_not_mere_presence():
    """Upstream gates branch A on score > 70. A wbaV hit at 55 must NOT enter it."""
    got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 55.0)))
    assert got["o_antigen_rule"].startswith("C_"), got


def test_wzy_below_its_own_gate_does_not_make_it_9_46():
    """Branch A's inner test is `wzy > 40`, separate from the wbaV gate."""
    got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 95.0), hit(WZY_946_FULL_KEY, 30.0)))
    assert got["antigen"] == "9", got


def test_the_tyr_O9_vs_O2_refinement_is_REFUSED_not_approximated():
    """THE REAL-DATA REGRESSION, pinned. Comparing the two `tyr` paralogs on BLAST metrics fired
    `O-2` on 16 isolates for 13 misses, 3 no-calls and ZERO hits -- every one a true O9 that SeqSero2
    called `9`. The references are ~99% identical to each other, so both align full length and the
    comparison is decided by noise (tyr-O-2's coverage even reads 100.2 on a gapped alignment). The
    branch default O-9 is taken instead, and the rule name says the discrimination did not run."""
    o2_looks_better = per_allele(hit(WBAV_KEYS[0], 95.0),
                                 hit("tyr-O-9__609", 100.0, axis="SP", antigen="tyr-O-9",
                                     identity=99.3),
                                 hit("tyr-O-2__608", 100.2, axis="SP", antigen="tyr-O-2",
                                     identity=99.2))
    got = call_o_antigen(o2_looks_better)
    assert got["antigen"] == "9", "a 0.2-point coverage artifact must not decide the antigen"
    assert got["o_antigen_rule"] == "A_wbaV_default_O9_tyr_not_discriminable_by_blast"


def test_the_refused_tyr_scores_are_still_reported_as_evidence():
    """Refusing to act on the numbers is not the same as hiding them -- a future k-mer-based
    refinement needs them, and a reader needs to see WHY the branch declined."""
    got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 95.0),
                                    hit("tyr-O-9__609", 100.0, axis="SP", antigen="tyr-O-9"),
                                    hit("tyr-O-2__608", 100.2, axis="SP", antigen="tyr-O-2")))
    assert set(got["tyr_scores"]) == {"tyr-O-9__609", "tyr-O-2__608"}


def test_no_branch_can_emit_O2_at_all():
    """The honest consequence, stated as a test rather than buried in a comment: a true O-2 genome
    (Paratyphi A / Nitra / Kiel / Koessen) WILL be called O-9. If a faithful k-mer refinement is ever
    added, this test is the one that should fail."""
    for cov in (40.0, 70.0, 99.9, 100.2):
        got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 95.0),
                                        hit("tyr-O-2__608", cov, axis="SP", antigen="tyr-O-2",
                                            identity=100.0)))
        assert got["antigen"] != "2", got


# --- Branch B: 3,10 vs 1,3,19 -----------------------------------------------------------------------

def _branch_b(*extra):
    return per_allele(hit(WZX_310_KEY, 95.0), hit(WZY_946_FULL_KEY, 95.0), *extra)


def test_the_reciprocal_marker_selects_3_10_positively():
    got = call_o_antigen(_branch_b(hit(NOT_IN_1319_KEY, 95.0)))
    assert got == {"antigen": "3,10", "o_antigen_rule": "B_not_in_1319_present"}


def test_1_3_19_is_reached_by_ABSENCE_not_by_a_marker_hit():
    """The asymmetry. `O-3,10_not_in_1,3,19` selects 3,10 when PRESENT; `1,3,19` is what you get when
    it is missing. A symmetric rule -- 'each not_in marker asserts its own antigen' -- would be wrong,
    and inferring one from the header name is exactly how this cell got its previous wrong causes."""
    got = call_o_antigen(_branch_b())
    assert got == {"antigen": "1,3,19", "o_antigen_rule": "B_not_in_1319_absent"}


def test_branch_b_needs_BOTH_of_its_preconditions():
    """wzx alone must fall through to the generic argmax, not call 1,3,19."""
    got = call_o_antigen(per_allele(hit(WZX_310_KEY, 95.0)))
    assert got["o_antigen_rule"].startswith("C_"), got


# --- Branch C: the argmax and its two guards --------------------------------------------------------

def test_the_1_3_19_marker_can_never_win_the_argmax_on_its_own():
    """THE MEASURED DEFECT. This 130 bp differential marker hit 60% of a diverse cohort at ~100%
    coverage against a 2% true prevalence, and the old caller emitted `1,3,19` every time. Upstream
    re-runs the argmax with this exact key excluded."""
    got = call_o_antigen(per_allele(hit(NOT_IN_310_KEY, 100.0), hit("O-7_wzy__1080", 80.0)))
    assert got["antigen"] == "7", got
    assert got["o_antigen_rule"] == "C_argmax_after_1319_exclusion"


def test_the_marker_alone_yields_NO_call_rather_than_a_confident_wrong_one():
    """With nothing else in the dict, upstream leaves the O call empty. An abstention is the honest
    output here; the old behaviour was a confident wrong antigen."""
    got = call_o_antigen(per_allele(hit(NOT_IN_310_KEY, 100.0)))
    assert got == {"antigen": None, "o_antigen_rule": "C_guard_1319_marker_only"}


def test_guard_1_rewrites_a_bare_wbaV_argmax_winner_to_O9():
    """A wbaV hit below branch A's gate can still win the argmax; with no wzy present it is O9, not
    9,46 -- the same discrimination, applied outside branch A."""
    got = call_o_antigen(per_allele(hit(WBAV_KEYS[0], 60.0)))
    assert got == {"antigen": "9", "o_antigen_rule": "C_guard_wbaV_without_wzy"}


def test_the_argmax_tiebreak_is_deterministic_and_independent_of_blast_report_order():
    """A genuine tie on BOTH metrics is broken by source index, never by whatever order blastn
    happened to report hits in -- otherwise the same genome could call differently between runs."""
    a = ("O__7__10__O-7_wzy__1080", {"called": True, "percent_identity": 99.0, "percent_coverage": 88.0})
    b = ("O__8__20__O-8_wzy__1200", {"called": True, "percent_identity": 99.0, "percent_coverage": 88.0})
    assert call_o_antigen(dict([a, b]))["antigen"] == "8"
    assert call_o_antigen(dict([b, a]))["antigen"] == "8"   # input order must not matter


def test_branch_C_ranks_identity_first_not_coverage():
    """THE OTHER REAL-DATA REGRESSION, pinned. A first cut ranked branch C on coverage alone. On six
    isolates a longer-but-worse-matching entry then out-covered the right one -- five where the
    multi-gene `O:23-gene2` beat `O-13_wzx`, one where `9,46` beat `3,10` -- each previously called
    correctly by this caller's measured identity-primary ranking. Porting upstream's BRANCHES does
    not mean porting its single-score ranking: BLAST gives two numbers and which one leads was
    already settled here by measurement."""
    got = call_o_antigen(per_allele(
        hit("O-13_wzx__1236", 92.0, antigen="13", identity=99.8),
        hit("O:23-gene2__1146", 100.0, antigen="23-gene2", identity=91.0)))
    assert got["antigen"] == "13", got
    # ...and coverage still breaks a tie on identity, which is the pre-existing ordering.
    tie = call_o_antigen(per_allele(
        hit("O-13_wzx__1236", 60.0, antigen="13", identity=99.0),
        hit("O-7_wzy__1080", 95.0, antigen="7", identity=99.0)))
    assert tie["antigen"] == "7", tie


def test_ordinary_alleles_are_unaffected():
    got = call_o_antigen(per_allele(hit("O-4_wzx__1293", 95.0), hit("O-7_wzy__1080", 60.0)))
    assert got == {"antigen": "4", "o_antigen_rule": "C_argmax"}


# --- degradation -------------------------------------------------------------------------------------

def test_uncalled_alleles_are_ignored():
    got = call_o_antigen(per_allele(hit("O-4_wzx__1293", 95.0), hit(WBAV_KEYS[0], 99.0, called=False)))
    assert got["antigen"] == "4"


def test_no_called_o_allele_abstains():
    assert call_o_antigen({})["o_antigen_rule"] == "no_o_allele_called"


def test_a_legacy_role_less_db_is_reported_not_silently_mis_branched():
    """The DB is gitignored, so a machine with an older build has 3-field headers. Every literal-key
    test would then miss and the procedure would report a confident branch-C result off an empty
    dict. It must say `legacy_db_no_role_field` instead, so the caller can fall back visibly."""
    legacy = {"O__1,3,19__126": {"called": True, "percent_identity": 99.0, "percent_coverage": 100.0}}
    assert call_o_antigen(legacy) == {"antigen": None, "o_antigen_rule": "legacy_db_no_role_field"}
