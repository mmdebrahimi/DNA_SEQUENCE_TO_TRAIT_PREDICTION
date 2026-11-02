"""Guards for the window-ordering probe.

Two properties matter and both are asserted on constructed cases rather than on real data, because
real data cannot distinguish "the orderings agree" from "the comparison is not measuring anything".

1. ORDER IS AN INPUT. The upstream extractor threads a running `last_end` through gene features in
   file order and gates each window on `start - last_end >= min_intergenic_length`, so the order
   decides which windows are emitted. A case built to make the two orderings differ must produce
   different window SETS; if it does not, the comparison is inert.
2. THE TWO UPSTREAM QUIRKS ARE REPLICATED, NOT CORRECTED. Each is pinned by a case where correcting
   it would change the result, so a future tidy-up fails loudly instead of silently measuring a rule
   nobody trained on.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.frame_window_ordering_probe import corpus_windows, coverage

GENOME = 10_000


def test_ordering_is_an_INPUT_not_a_presentation_detail():
    """The decisive non-vacuity check for the whole probe."""
    ordered = [(300, 400), (402, 500), (600, 700)]
    shuffled = [(402, 500), (300, 400), (600, 700)]

    w_sorted = set(corpus_windows(sorted(ordered), [], GENOME))
    w_file = set(corpus_windows(shuffled, [], GENOME))

    assert w_sorted != w_file, "the comparison is inert: both orderings gave the same window set"
    assert len(w_sorted ^ w_file) > 0


def test_a_gene_skipped_by_the_edge_guard_does_NOT_update_last_end():
    """Replicated quirk #1, pinned by a case where correcting it changes the answer.

    Upstream, `if start < regulatory_region_length: continue` precedes `last_end = end`, so a gene
    too close to the sequence start leaves `last_end` untouched. Here the skipped gene ends at 900
    and the next starts at 902: under the quirk `last_end` is still 0 so the gap test passes and a
    window IS emitted; under a "corrected" version `last_end` would be 900, the gap would be 2, and
    the window would be dropped.
    """
    genes = [(150, 900), (902, 1000)]
    got = corpus_windows(genes, [], GENOME)
    assert len(got) == 1, got
    assert got[0][:2] == (702, 902), got


def test_the_reverse_branch_compares_against_the_previous_genes_START():
    """Replicated quirk #2. The reverse branch tests `last_end_rev - end`, where `last_end_rev` holds
    the previous gene's START rather than its end. Asserted by shape: iterating the reverse list from
    the genome end, the emitted windows sit DOWNSTREAM of each gene end and are reverse-oriented."""
    rev = [(3000, 3100), (5000, 5100)]
    got = corpus_windows([], rev, GENOME)
    assert got, "no reverse-strand windows emitted at all"
    for a, b, strand in got:
        assert strand == "-"
        assert b - a == 200
    # the window for the gene ending at 5100 begins exactly there, not 200 before it
    assert (5100, 5300, "-") in got, got


def test_min_intergenic_gate_rejects_a_packed_gene():
    genes = [(300, 400), (403, 500)]
    got = corpus_windows(sorted(genes), [], GENOME)
    assert len(got) == 1, got
    assert got[0][:2] == (100, 300)


def test_edge_guard_rejects_a_gene_too_close_to_the_start():
    assert corpus_windows([(199, 400)], [], GENOME) == []
    assert corpus_windows([(200, 400)], [], GENOME) == [(0, 200, "+")]


def test_coverage_mask_counts_each_base_once_even_when_windows_overlap():
    """Overlapping windows must not double-count, or the covered-base total inflates."""
    w = [(100, 300, "+"), (200, 400, "+")]
    assert sum(coverage(w, GENOME)) == 300


def test_coverage_clamps_at_the_genome_boundaries():
    assert sum(coverage([(-50, 50, "+")], GENOME)) == 50
    assert sum(coverage([(GENOME - 50, GENOME + 50, "+")], GENOME)) == 50


ART = Path("wiki/frame_plan_window_ordering_2026-11-01.json")


@pytest.mark.skipif(not ART.is_file(), reason="probe has not been run on this host")
def test_the_real_run_recorded_a_divergence_and_named_its_verdict():
    """The real-data run is the deliverable. A zero difference is a legitimate outcome, so this
    asserts the artifact is internally consistent rather than asserting a particular direction."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["verdict"] in {"ORDERING_DIVERGES", "ORDERING_IMMATERIAL"}
    divergent = d["symmetric_difference_windows"] > 0 or d["tiles_changing_in_window_status"] > 0
    assert (d["verdict"] == "ORDERING_DIVERGES") == divergent
    # the sorted figures are the ones the frame plan pins, so they must be the larger-covered set
    assert d["sorted_order"]["n_windows"] >= d["file_order"]["n_windows"]
    assert d["honest_limits"], "an artifact with no stated limits is over-claiming"


# --------------------------------------------------------------------------------------------------
# gaps found by /test-epilogue
# --------------------------------------------------------------------------------------------------
def test_BOTH_gate_parameters_actually_move_the_output():
    """Every existing case runs at the module defaults (RL=200, MIN_IG=5), so neither keyword is
    proven to be wired. A parameter that is accepted and ignored turns a sweep into an invariance
    result that reads as a finding -- the probe records both values in its artifact, so a reader
    would take them as describing what was computed."""
    genes = [(300, 400), (410, 600), (700, 900)]

    base = corpus_windows(genes, [], GENOME)
    wider_region = corpus_windows(genes, [], GENOME, regulatory_region_length=100)
    wider_gap = corpus_windows(genes, [], GENOME, min_intergenic_length=50)

    assert base != wider_region, "regulatory_region_length does not move the output"
    assert base != wider_gap, "min_intergenic_length does not move the output"
    # and each moves it in its own distinct way: the region length changes window WIDTH, the
    # intergenic floor changes HOW MANY windows survive.
    assert {b - a for a, b, _ in wider_region} == {100}
    assert {b - a for a, b, _ in base} == {200}
    assert len(wider_gap) < len(base)


def test_the_reverse_edge_guard_rejects_a_gene_too_close_to_the_genome_END():
    """Mirror of `test_edge_guard_rejects_a_gene_too_close_to_the_start`, which covered only the
    forward branch. The reverse branch's guard is `genome_len - end < regulatory_region_length`, a
    different expression against a different boundary, and nothing exercised it."""
    assert corpus_windows([], [(9800, 9900)], GENOME) == []        # 100 bp of room, needs 200
    assert corpus_windows([], [(9700, 9801)], GENOME) == []        # 199 bp, still short
    assert corpus_windows([], [(9700, 9800)], GENOME) == [(9800, 10_000, "-")]  # exactly 200


def test_the_reverse_min_intergenic_gate_rejects_a_PACKED_pair():
    """`test_min_intergenic_gate_rejects_a_packed_gene` only exercises the forward gate. The reverse
    gate compares `last_end_rev - end`, where `last_end_rev` is the previous gene's START -- a
    different quantity -- so its rejection path was untested."""
    packed = corpus_windows([], [(3000, 3100), (3102, 3200)], GENOME)
    assert packed == [(3200, 3400, "-")], packed      # second gene dropped: 3102 - 3100 = 2 < 5
    spaced = corpus_windows([], [(3000, 3100), (3200, 3300)], GENOME)
    assert len(spaced) == 2, spaced                   # 3200 - 3100 = 100, both survive


def test_rc_is_an_involution_and_passes_non_ACGT_through_unchanged():
    """`rc` is how the probe finds minus-strand tiles, and a wrong `rc` halves the hit set silently
    -- the run would still print a number. The N case matters because the genome is upper-cased but
    not filtered: `str.maketrans("ACGT", "TGCA")` leaves N in place, so positions stay aligned
    rather than the string shortening and shifting every coordinate after it."""
    from scripts.frame_window_ordering_probe import rc

    assert rc("ACGT") == "ACGT"
    assert rc("AAAC") == "GTTT"
    assert rc("A") == "T" and rc("") == ""
    assert rc("ANNC") == "GNNT"          # length preserved, N not dropped
    for s in ("ACGTTGCAN", "AAGGCCTTN", ""):
        assert rc(rc(s)) == s, s


def test_the_minus_strand_position_formula_round_trips_to_FORWARD_coordinates():
    """The probe's single load-bearing coordinate conversion, and the one nothing covered.

    `main` records a minus-strand hit found at offset `p` in `rc(G)` as `n - p - TILE_LEN` in
    forward coordinates. An off-by-one there does not fail -- it shifts every minus-strand tile by
    one base, which quietly moves tiles across window boundaries and so changes
    `tiles_changing_in_window_status`, the artifact's headline number.

    Derived independently: a tile is built as the reverse complement of a KNOWN forward slice, so
    `q` is known before the formula is applied rather than being read back out of it.
    """
    import random

    from scripts.frame_window_ordering_probe import rc

    random.seed(11)
    fwd = "".join(random.choice("ACGT") for _ in range(600))
    tile_len, q = 40, 321
    tile = rc(fwd[q:q + tile_len])

    rev = rc(fwd)
    n = len(fwd)
    assert rev.count(tile) == 1, "the fixture is ambiguous; the round trip would be untestable"
    p = rev.index(tile)

    assert n - p - tile_len == q
    # the inverse check, which is what makes the coordinate meaningful: the forward genome at the
    # reported position really does hold this tile's reverse complement.
    assert fwd[q:q + tile_len] == rc(tile)
    # and the naive-but-wrong forms are genuinely different, so the assertion above is not an
    # accident of the fixture
    assert n - p != q and n - p - tile_len + 1 != q


def test_coverage_of_no_windows_is_zero_and_strand_is_IGNORED():
    """`coverage` collapses both strands onto forward coordinates -- the third tuple element is
    discarded. That is correct (an in-window test asks whether a locus is covered, not by which
    strand) and it is load-bearing: if strand were keyed, a locus covered from both strands would
    count twice and `bases_covered` would exceed the genome length."""
    assert sum(coverage([], GENOME)) == 0
    assert sum(coverage([(10, 20, "+"), (10, 20, "-")], GENOME)) == 10
    assert sum(coverage([(10, 20, "+"), (15, 25, "-")], GENOME)) == 15


@pytest.mark.skipif(not ART.is_file(), reason="probe has not been run on this host")
def test_the_real_artifact_is_internally_RECONCILED_not_just_well_formed():
    """The existing artifact test checks the verdict agrees with the two divergence fields. These
    are the arithmetic relationships that must also hold, each of which a plumbing fault would
    break while leaving the verdict consistent.

    `n_distinct_windows == n_windows` is the one that matters most: the whole measurement is a
    symmetric difference of SETS, so a duplicate window would be invisible to it and the reported
    divergence would understate.
    """
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["schema"] == "frame-plan-window-ordering-v1"
    for side in ("file_order", "sorted_order"):
        s = d[side]
        assert s["n_distinct_windows"] == s["n_windows"], (side, s)
        # a window is regulatory_region_length wide, so covered bases can never exceed n * width
        assert s["bases_covered"] <= s["n_windows"] * d["regulatory_region_length"], (side, s)
        assert s["bases_covered"] > 0 and s["bases_covered"] <= d["genome_length"]
        assert s["tiles_in_window"] <= d["uniquely_located_tiles"], (side, s)
    # the symmetric difference bounds how far the two window COUNTS can be apart
    delta = abs(d["sorted_order"]["n_windows"] - d["file_order"]["n_windows"])
    assert delta <= d["symmetric_difference_windows"], (delta, d["symmetric_difference_windows"])
    # the ordering claim the whole probe rests on: if file order WERE ascending on both strands the
    # two programs would be identical and a non-zero difference would be a bug, not a finding
    if d["verdict"] == "ORDERING_DIVERGES":
        assert not all(d["file_order_is_ascending_by_start"].values()), d[
            "file_order_is_ascending_by_start"]


PLAN = Path("plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan"
            "/technical-plan.md")


@pytest.mark.skipif(not (ART.is_file() and PLAN.is_file()),
                    reason="probe artifact or the plan that cites it is absent on this host")
def test_the_PLAN_PROSE_still_quotes_the_numbers_the_ARTIFACT_actually_holds():
    """The gap nothing else covers: the artifact's figures are RESTATED in plan prose.

    `technical-plan.md` quotes all seven of these as literals ("3,496 windows covering 625,255
    bases ... symmetric difference of 26 windows ... 146 of 43,061 uniquely-located tiles"). The
    probe writes its artifact unconditionally, so re-running it on a different annotation build --
    or after any change to `corpus_windows` -- moves the numbers while leaving the prose asserting
    the old ones, with nothing failing. That is the drift this repo's own contract-number audit
    exists to catch, and the audit only looks at cell contracts, not at plan documents.

    Asserting PRESENCE in the prose rather than freezing the artifact is deliberate: the artifact is
    a measurement and must be free to change. What must not happen is it changing SILENTLY.

    MATCHING IS WHOLE-NUMBER, and that is not a stylistic choice. The first draft of this test used
    a plain `str(v) in text`, which this repo has already measured to be the dominant defect in its
    own contract-number audit: a short token matches INSIDE a longer number. Measured here -- `27`
    (a decoy for the symmetric difference of 26) was "found" in the prose under substring matching
    and is correctly ABSENT under boundary-guarded matching. A test that cannot distinguish 26 from
    27 cannot detect the drift it exists for.
    """
    import re

    d = json.loads(ART.read_text(encoding="utf-8"))
    text = PLAN.read_text(encoding="utf-8")

    def quoted_whole(v: int) -> bool:
        """`v` appears in the prose as a whole number, in either bare or comma-grouped form."""
        return any(
            re.search(r"(?<![0-9,.])" + re.escape(form) + r"(?![0-9,.])", text)
            for form in (format(v, ","), str(v))
        )

    quoted = {
        "sorted n_windows": d["sorted_order"]["n_windows"],
        "sorted bases_covered": d["sorted_order"]["bases_covered"],
        "file n_windows": d["file_order"]["n_windows"],
        "file bases_covered": d["file_order"]["bases_covered"],
        "symmetric_difference_windows": d["symmetric_difference_windows"],
        "tiles_changing_in_window_status": d["tiles_changing_in_window_status"],
        "uniquely_located_tiles": d["uniquely_located_tiles"],
    }
    missing = {k: v for k, v in quoted.items() if not quoted_whole(v)}
    assert not missing, (
        f"{PLAN} no longer quotes these artifact figures: {missing}. Either the probe was re-run "
        f"and the numbers moved -- in which case the plan's prose is now stale and must be "
        f"corrected, not this test -- or the plan was rewritten and the citation was dropped."
    )

    # DISCRIMINATION, measured rather than assumed: each figure offset by one must NOT be found, or
    # the check above would pass against prose quoting different numbers.
    decoys = [v + 1 for v in quoted.values()] + [v - 1 for v in quoted.values()]
    hits = [v for v in decoys if quoted_whole(v)]
    assert not hits, (
        f"these near-miss values are also present in the prose: {sorted(hits)}; the citation check "
        f"above is not discriminating and would pass on drifted numbers"
    )
