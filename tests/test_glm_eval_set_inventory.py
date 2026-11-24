"""Tests for the eval-set inventory probe.

Two loads are borne here. The `set_algebra` tests are offline and synthetic: they pin that the
orthology filter cannot silently fall back to the unfiltered set, which is the failure that would
make `|C|` look healthy while measuring set B. The drift-guard test proves the `|A| == 3349` refusal
FIRES rather than merely existing -- a guard nobody has seen reject anything is not a guard.

The real-data tests are skip-if-absent so the suite stays green on a host without the refseq cache,
and the one that matters re-derives the committed 1.4291 from the current code. That number is the
whole basis of the committed refusal, so if it stops reproducing, the decomposition that calls it a
sampling artifact has lost its anchor.
"""
from __future__ import annotations

import importlib.util
import json
import math
import statistics
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
ARTIFACT = ROOT / "wiki" / "glm_eval_set_inventory_2026-11-23.json"

needs_data = pytest.mark.skipif(
    not ((MG1655 / "annotations.gff3").exists() and CORPUS.exists()),
    reason="MG1655 refseq cache or training corpus absent on this host")


def _load():
    p = ROOT / "scripts" / "glm_eval_set_inventory.py"
    spec = importlib.util.spec_from_file_location("glm_eval_set_inventory", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["glm_eval_set_inventory"] = mod
    spec.loader.exec_module(mod)
    return mod


INV = _load()


class _W:
    """A stand-in for a corpus window: just the two fields the filters read."""

    def __init__(self, seq, gene_id):
        self.seq = seq
        self.gene_id = gene_id


# ---- the set algebra, offline and synthetic --------------------------------------------------

def test_sizes_and_overlaps_on_a_hand_built_case():
    got = INV.set_algebra(["AAA", "CCC"], ["AAA", "GGG", "TTT"], ["GGG"])
    assert got["A_sequence_withheld"] == 2
    assert got["B_held_out_of_corpus"] == 3
    assert got["C_orthology_free"] == 1
    assert got["A_inter_B"] == 1
    assert got["B_inter_C"] == 1
    assert got["A_inter_C"] == 0


def test_C_must_be_a_subset_of_B_and_a_violation_RAISES():
    """C is B filtered by gene symbol. A member of C outside B means the filter ran on the wrong
    list, which would otherwise be reported as an interesting overlap rather than a bug."""
    with pytest.raises(INV.InventoryRefusal):
        INV.set_algebra(["AAA"], ["AAA"], ["GGG"])


def test_an_empty_only_set_yields_C_of_zero_and_never_falls_back_to_B():
    """The failure this guards: an empty symbol difference silently leaving C == B, so a fully
    unfiltered set would be reported as orthology-free."""
    held = [_W("AAA", "g1"), _W("CCC", "g2")]
    of = [w.seq.upper() for w in held if (({}.get(w.gene_id) or "")) in set()]
    assert of == []
    got = INV.set_algebra(["AAA"], [w.seq for w in held], of)
    assert got["C_orthology_free"] == 0
    assert got["B_held_out_of_corpus"] == 2


def test_the_symbol_difference_filter_keeps_only_mg1655_only_genes():
    held = [_W("AAA", "g1"), _W("CCC", "g2"), _W("GGG", "g3")]
    m = {"g1": "thrL", "g2": "shared", "g3": "prophageX"}
    only = {"thrL", "prophageX"}
    of = [w.seq.upper() for w in held if (m.get(w.gene_id) or "") in only]
    assert of == ["AAA", "GGG"]


def test_case_is_normalised_so_a_case_mismatch_cannot_read_as_disjointness():
    """`fit_eval_halves` returns raw `w.seq` while the scorer uppercases. Unnormalised, these two
    sets would intersect at ZERO and that would look like a finding."""
    got = INV.set_algebra(["aaa"], ["AAA"], ["AAA"])
    assert got["A_inter_B"] == 1


def test_the_filters_measured_cost_is_C_over_A():
    got = INV.set_algebra(["A" * 3, "C" * 3, "G" * 3, "T" * 3], ["A" * 3], ["A" * 3])
    assert got["C_over_A_the_filters_measured_cost"] == pytest.approx(0.25)


# ---- the drift guard must FIRE -----------------------------------------------------------------

def test_the_drift_guard_fires_on_3348():
    """Proven by making it reject, not by asserting it exists."""
    with pytest.raises(INV.InventoryRefusal) as e:
        INV.assert_expected_a(3348)
    assert "3348" in str(e.value) and "3349" in str(e.value)


def test_the_drift_guard_passes_on_the_manifest_size():
    INV.assert_expected_a(INV.EXPECTED_A_SIZE)
    assert INV.EXPECTED_A_SIZE == 3349


def test_the_drift_guard_is_not_vacuous_on_a_near_miss_either():
    for n in (0, 1, 3348, 3350, 10_000):
        with pytest.raises(INV.InventoryRefusal):
            INV.assert_expected_a(n)


# ---- the committed geometry this probe reconstructs --------------------------------------------

def test_the_reconstructed_committed_geometry_matches_the_committed_artifact():
    """The artifact records n_windows and positions_per_window; the probe must agree with it."""
    art = json.loads((ROOT / "wiki" / "glm_beyond_kmer_2026-11-22.json").read_text(encoding="utf-8"))
    assert INV.COMMITTED_N_WINDOWS == art["eval"]["n_windows"]
    assert INV.COMMITTED_N_POSITIONS == art["eval"]["positions_per_window"]


def test_uniform_constant_agrees_with_the_gate_module():
    from dna_decode.glm.eval_admissibility import UNIFORM_CE_NATS_DNA
    assert INV.UNIFORM == UNIFORM_CE_NATS_DNA == math.log(4)


# ---- real data -----------------------------------------------------------------------------------

@needs_data
def test_set_A_recomputes_to_the_manifest_size():
    """The manifest records only `withheld_union_size`; set A itself is not persisted anywhere."""
    from dna_decode.glm.training_corpus import falsifier_eval_windows
    a = falsifier_eval_windows(INV.MG1655, length=INV.LENGTH, limit=INV.SCORING_LIMIT,
                               seeds=range(INV.N_SEEDS), kind=INV.KIND)
    manifest = json.loads(
        (ROOT / "wiki" / "glm_training_corpus_manifest_2026-11-21.json").read_text(encoding="utf-8"))
    assert len(a) == manifest["leak_gate"]["withheld_union_size"] == 3349


@needs_data
def test_the_committed_1_4291_reproduces_from_the_current_code():
    """The anchor of the whole decomposition. If this drifts, calling 1.4291 a sampling artifact
    loses its basis, because the claim rests on reproducing it and then varying one factor."""
    from dna_decode.glm.baselines import MarkovGenerator
    from scripts.glm_beyond_kmer import (ORDER, ecoli_fit_windows, load_corpus_windows,
                                         mg1655_candidate_sets)
    cset = set(load_corpus_windows())
    chain = MarkovGenerator(order=ORDER, seed=0).fit(ecoli_fit_windows(cset, length=INV.LENGTH))
    of = mg1655_candidate_sets(cset, length=INV.LENGTH)["orthology_free"]
    d = INV.decompose_committed_cell(chain, of, length=INV.LENGTH)
    assert d["cells"]["committed_6_windows_8_sampled_positions"] == 1.4291


@needs_data
def test_the_decomposition_crosses_uniform_when_only_the_sample_size_changes():
    """The committed cell is above uniform; hold the windows fixed, score every position, and it is
    below. That contrast is the finding, so it is pinned rather than described."""
    from dna_decode.glm.baselines import MarkovGenerator
    from scripts.glm_beyond_kmer import (ORDER, ecoli_fit_windows, load_corpus_windows,
                                         mg1655_candidate_sets)
    cset = set(load_corpus_windows())
    chain = MarkovGenerator(order=ORDER, seed=0).fit(ecoli_fit_windows(cset, length=INV.LENGTH))
    of = mg1655_candidate_sets(cset, length=INV.LENGTH)["orthology_free"]
    c = INV.decompose_committed_cell(chain, of, length=INV.LENGTH)["cells"]
    assert c["committed_6_windows_8_sampled_positions"] > INV.UNIFORM
    assert c["same_6_windows_ALL_positions"] < INV.UNIFORM
    assert c["ALL_windows_ALL_positions"] < INV.UNIFORM


@needs_data
def test_set_C_is_far_larger_than_the_committed_run_scored():
    """The premise this probe exists to test: the filter was believed to leave too little to score."""
    from scripts.glm_beyond_kmer import load_corpus_windows, mg1655_candidate_sets
    of = mg1655_candidate_sets(set(load_corpus_windows()), length=INV.LENGTH)["orthology_free"]
    assert len(of) == 648
    # The scorer refuses when len(of) < --n-windows, so its own default would have run.
    assert len(of) >= 120


# ---- the artifact ---------------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_the_artifact_parses_as_json():
    """This repo has shipped an unparseable hand-edited artifact before."""
    json.loads(ARTIFACT.read_text(encoding="utf-8"))


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_the_artifact_emits_no_generator_verdict():
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert art["emits_no_verdict_about_the_generator"] is True
    assert "verdict" not in art
    blob = json.dumps(art).upper()
    assert "BEYOND_KMER_CONSERVATIVE_POSITIVE" not in blob
    assert "WORSE_THAN_ORDER3" not in blob


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_every_number_in_the_headline_appears_in_a_structured_field():
    """The contract-number-audit discipline: prose figures must be re-derivable from fields."""
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    d = art["committed_cell_decomposition"]
    s = art["sizes_and_overlaps"]
    # The headline quotes BOTH, because the scorer slices the list while statistics want distinct
    # windows. Conflating them is what this assertion caught the first time it ran.
    assert s["C_list_entries"] == 648
    assert s["C_orthology_free"] == 634
    assert s["C_duplicate_entries"] == 14
    assert s["B_list_entries"] == 3923
    assert s["B_held_out_of_corpus"] == 3906
    assert s["B_duplicate_entries"] == 17
    assert d["cells"]["committed_6_windows_8_sampled_positions"] == 1.4291
    assert d["cells"]["ALL_windows_ALL_positions"] == 1.3571
    assert d["uniform_ce_nats"] == 1.3863
    assert d["n_predictions_in_the_committed_cell"] == 48
    assert art["committed_run_geometry_reconstructed"]["windows_that_were_AVAILABLE"] == 648


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_the_slices_the_scorer_actually_uses_are_duplicate_free():
    """The duplicate finding is a LATENT hazard only because these two slices are clean. If either
    stops being clean, a scored window would be counted twice."""
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    ss = art["slice_safety"]
    assert ss["first_6_is_duplicate_free"] is True
    assert ss["first_120_is_duplicate_free"] is True
    assert ss["distinct_in_first_120"] == 120
    # And the full list is NOT duplicate-free -- otherwise the guard above is vacuous.
    assert ss["first_648_is_duplicate_free"] is False


# ================================================================================================
# `null_ce_on` and `null_side_leak_contrast`: the two scoring helpers, previously untested.
#
# Both are pure given a chain, and `markov_nll` reads ONLY `chain._counts`, so a two-line stub
# exercises them exactly. That matters because `null_side_leak_contrast` produces the 0.0007-nat
# null-side leak figure that `wiki/glm_beyond_kmer_two_sets_2026-11-23.json` quotes and reasons
# from, and its SIGN decides which party the near-ortholog leak favours. A flipped sign there would
# invert a published claim while every number still looked plausible.
# ================================================================================================

from types import SimpleNamespace  # noqa: E402

# p(A | AAA) = 97/100, so one such prediction costs -ln(0.97) nats. Every other context is absent
# from the table, so it backs off to the add-one uniform and costs exactly ln(4).
PREDICTABLE = SimpleNamespace(_counts={"AAA": [97, 1, 1, 1]})
EMPTY_TABLE = SimpleNamespace(_counts={})

EASY = "A" * 10          # every scoreable position is A|AAA -> -ln(0.97)
HARD = "ACGTACGTAC"      # no position's context is in the table -> ln(4) throughout

# A window whose positions do NOT all cost the same, which is what makes "scores ALL positions"
# a checkable claim: over positions 3..9 the costs are -ln(0.97), -ln(0.97), -ln(0.01) and then
# ln(4) four times, so the all-position mean is a value no subset of positions reproduces.
HETEROGENEOUS = "AAAAACGTAC"
HET_ALL_POSITION_MEAN = 1.4588


def test_null_ce_on_computes_a_hand_checkable_cross_entropy():
    got = INV.null_ce_on(PREDICTABLE, [EASY], length=10)
    assert got["null_ce_nats"] == round(-math.log(0.97), 4) == 0.0305
    assert got["n_windows"] == 1


def test_null_ce_on_really_averages_EVERY_position_not_a_subset():
    """THE claim that distinguishes this helper from the scorer, on a window whose positions cost
    different amounts. The homogeneous fixtures above cannot test it: if the implementation scored
    a sample of positions they would return the identical value, so they pin the arithmetic while
    this pins the GEOMETRY. 1.4588 is the mean over all seven positions; the first three alone give
    1.5554 and the last four give 1.3863, so neither a head slice nor a tail slice reaches it."""
    got = INV.null_ce_on(PREDICTABLE, [HETEROGENEOUS], length=10)
    assert got["null_ce_nats"] == HET_ALL_POSITION_MEAN
    hand = [-math.log(0.97), -math.log(0.97), -math.log(0.01)] + [math.log(4)] * 4
    assert got["null_ce_nats"] == round(statistics.fmean(hand), 4)
    assert got["n_positions_per_window"] == 7 == len(hand)


def test_null_ce_on_means_over_WINDOWS_after_meaning_over_positions():
    """Two windows of differing difficulty: the reported figure is the mean of the per-window means.
    With the heterogeneous window beside the easy one, a pooled mean over all positions would give
    the same answer only because both windows have equal position counts -- so this pins the pair
    rather than claiming to separate the two estimators."""
    got = INV.null_ce_on(PREDICTABLE, [HETEROGENEOUS, EASY], length=10)
    expected = statistics.fmean([HET_ALL_POSITION_MEAN, -math.log(0.97)])
    assert got["null_ce_nats"] == pytest.approx(expected, abs=5e-5)
    assert got["n_windows"] == 2


def test_the_leak_contrast_also_averages_every_position_of_each_half():
    """Same claim on the leak contrast, whose halves were otherwise only ever homogeneous."""
    got = INV.null_side_leak_contrast(PREDICTABLE, [HETEROGENEOUS, HARD], [HARD], length=10)
    assert got["A_minus_C_ortholog_bearing"]["null_ce_nats"] == HET_ALL_POSITION_MEAN
    assert got["A_inter_C_orthology_free"]["null_ce_nats"] == 1.3863


def test_null_ce_on_scores_every_position_from_ORDER_to_length():
    """The geometry the two-set artifact's `why_147` field leans on: 150 - order 3 = 147."""
    from scripts.glm_beyond_kmer import ORDER
    for length in (10, 20, 150):
        got = INV.null_ce_on(PREDICTABLE, ["A" * length], length=length)
        assert got["n_positions_per_window"] == length - ORDER
    assert INV.null_ce_on(PREDICTABLE, ["A" * 150], length=150)["n_positions_per_window"] == 147


def test_null_ce_on_raises_when_length_OVERRUNS_the_sequences():
    """`length` is an ARGUMENT here while the scorer derives it from the checkpoint config, so the
    two can disagree. BOTH directions must refuse, and the under-run direction is the one that
    matters: it used to score fewer positions silently and report the short count as
    `n_positions_per_window`, i.e. a quietly different measurement under the full-geometry label.
    The overrun direction previously surfaced only as an `IndexError` from indexing; it now refuses
    by name, before any scoring, with the expected length in the message."""
    # OVERRUN: windows shorter than `length`
    with pytest.raises(INV.InventoryRefusal) as over:
        INV.null_ce_on(PREDICTABLE, ["A" * 10], length=150)
    assert "150" in str(over.value)
    # UNDER-RUN: windows longer than `length` -- the direction that used to pass silently
    with pytest.raises(INV.InventoryRefusal) as under:
        INV.null_ce_on(PREDICTABLE, ["A" * 150], length=10)
    assert "10" in str(under.value)
    # and a matching length is still scored normally, so the guard is not blanket
    assert INV.null_ce_on(PREDICTABLE, ["A" * 10], length=10)["n_windows"] == 1


def test_null_ce_on_returns_None_rather_than_crashing_on_an_empty_set():
    """A `statistics.fmean([])` would raise; a 0.0 would read as a perfect null. Neither happens."""
    got = INV.null_ce_on(PREDICTABLE, [], length=10)
    assert got["null_ce_nats"] is None
    assert got["n_windows"] == 0


def test_an_unseen_context_backs_off_to_uniform_and_the_gate_then_REFUSES():
    """End-to-end on the two modules: a chain that has learned nothing scores exactly uniform, and
    the gate's strict comparison must refuse it rather than admitting a set on a chance-level null.
    This is the composition the probe performs per set, and it was untested."""
    from dna_decode.glm.eval_admissibility import REFUSED_NULL_BELOW_UNIFORM, admissibility
    got = INV.null_ce_on(EMPTY_TABLE, [HARD, EASY], length=10)
    assert got["null_ce_nats"] == 1.3863 == round(INV.UNIFORM, 4)
    gate = admissibility(null_ce_nats=got["null_ce_nats"])
    assert gate.state == REFUSED_NULL_BELOW_UNIFORM
    assert not gate.admissible


def test_null_ce_on_consumes_no_randomness():
    """The documented difference from the scorer: no sampling, so no seed. If a draw crept in, the
    probe figures would stop being comparable across runs while still looking clean."""
    import random
    state = random.getstate()        # restored, so this cannot perturb any later test
    try:
        random.seed(1234)
        first = INV.null_ce_on(PREDICTABLE, [EASY, HARD], length=10)
        random.seed(999)
        assert INV.null_ce_on(PREDICTABLE, [EASY, HARD], length=10) == first
    finally:
        random.setstate(state)
    assert "no seed involved" in first["scoring"]


# ---- the leak contrast, where the SIGN is the finding -------------------------------------------

def test_the_leak_is_POSITIVE_when_the_null_scores_the_ortholog_bearing_half_BETTER():
    """The documented convention: positive = the null gained from the leak. Here A\\C is the EASY
    window (low cross-entropy = better) and A&C is the HARD one, so the sign must be positive."""
    got = INV.null_side_leak_contrast(PREDICTABLE, [EASY, HARD], [HARD], length=10)
    assert got["A_minus_C_ortholog_bearing"]["null_ce_nats"] == 0.0305
    assert got["A_inter_C_orthology_free"]["null_ce_nats"] == 1.3863
    assert got["null_side_leak_nats"] > 0
    assert got["null_side_leak_nats"] == round(math.log(4) + math.log(0.97), 4)


def test_the_leak_is_NEGATIVE_when_the_bearing_half_is_the_harder_one():
    """The mirror, so the test above cannot pass on a sign-blind implementation."""
    got = INV.null_side_leak_contrast(PREDICTABLE, [EASY, HARD], [EASY], length=10)
    assert got["A_minus_C_ortholog_bearing"]["null_ce_nats"] == 1.3863
    assert got["A_inter_C_orthology_free"]["null_ce_nats"] == 0.0305
    assert got["null_side_leak_nats"] < 0


def test_the_sign_convention_string_matches_the_direction_the_code_computes():
    """Prose and arithmetic must agree; this artifact reader is told how to read the number."""
    got = INV.null_side_leak_contrast(PREDICTABLE, [EASY, HARD], [HARD], length=10)
    assert "positive" in got["sign_convention"] and "BEARING half better" in got["sign_convention"]
    assert got["null_side_leak_nats"] > 0


def test_an_empty_half_yields_None_and_NOT_a_zero_leak():
    """A 0.0 would read as 'measured, and the leak is nil' -- a finding -- when nothing was
    measured at all. Each half's own CE is None too, so the emptiness is visible, not inferred."""
    both_free = INV.null_side_leak_contrast(PREDICTABLE, [EASY], [EASY], length=10)
    assert both_free["A_minus_C_ortholog_bearing"]["n_windows"] == 0
    assert both_free["A_minus_C_ortholog_bearing"]["null_ce_nats"] is None
    assert both_free["null_side_leak_nats"] is None

    none_free = INV.null_side_leak_contrast(PREDICTABLE, [EASY, HARD], [], length=10)
    assert none_free["A_inter_C_orthology_free"]["null_ce_nats"] is None
    assert none_free["null_side_leak_nats"] is None


def test_the_two_halves_partition_set_A_and_nothing_else_enters_either():
    """A window in C but NOT in A belongs to neither half -- the contrast is defined strictly
    inside the withheld set, so an outside window must not be silently scored."""
    outsider = "GGGGGGGGGG"
    got = INV.null_side_leak_contrast(PREDICTABLE, [EASY, HARD], [HARD, outsider], length=10)
    assert got["A_minus_C_ortholog_bearing"]["n_windows"] == 1
    assert got["A_inter_C_orthology_free"]["n_windows"] == 1


def test_the_leak_contrast_uppercases_so_a_case_mismatch_cannot_empty_a_half():
    """The same trap the module documents for `set_algebra`: `fit_eval_halves` returns raw `w.seq`.
    Unnormalised, a lowercase A against an uppercase C makes A&C empty, which would report
    `null_side_leak_nats: None` and read as 'this could not be measured'."""
    got = INV.null_side_leak_contrast(PREDICTABLE, [EASY.lower(), HARD.lower()], [HARD], length=10)
    assert got["A_inter_C_orthology_free"]["n_windows"] == 1
    assert got["null_side_leak_nats"] is not None


def test_the_leak_contrast_is_deterministic_and_says_it_uses_no_sampling():
    a = INV.null_side_leak_contrast(PREDICTABLE, [EASY, HARD], [HARD], length=10)
    b = INV.null_side_leak_contrast(PREDICTABLE, [HARD, EASY], [HARD], length=10)
    assert a["null_side_leak_nats"] == b["null_side_leak_nats"]
    assert "no sampling, no seed" in a["scoring"]


# ---- the position draw inside the decomposition ------------------------------------------------

def test_the_decomposition_advances_ONE_rng_across_windows():
    """The identity that makes the committed 1.4291 reproduce. `decompose_committed_cell` must draw
    each window positions from a SINGLE `random.Random(seed)` advanced across windows -- if it were
    re-seeded per window, every window would get the SAME positions and the committed figure would
    not come out of the first cell."""
    import random
    draws = []
    real = random.Random.sample

    def spy(self, population, k):
        got = real(self, population, k)
        draws.append(list(got))
        return got

    seqs = [("ACGT" * 5)[:20] for _ in range(6)]
    try:
        random.Random.sample = spy
        INV.decompose_committed_cell(EMPTY_TABLE, seqs, length=20, n_windows=6, n_positions=8)
    finally:
        random.Random.sample = real

    # Two sampled cells are computed: the first six windows, then all six -> 12 draws of 8.
    assert len(draws) == 12 and all(len(d) == 8 for d in draws)
    # Compared against the DRAWS of a single advanced RNG, not against object identity: a
    # per-window `random.Random(seed)` is a short-lived temporary and CPython reuses its address,
    # so an id()-based check can collapse to one value and pass on the very bug it names.
    expect = random.Random(0)
    assert draws[:6] == [expect.sample(list(range(3, 20)), 8) for _ in range(6)]
    # NON-VACUITY: a per-window reseed would make all six draws IDENTICAL, so the assertion above
    # is only meaningful because they are not.
    assert len({tuple(d) for d in draws[:6]}) > 1, "the draws are identical, so nothing was pinned"


# ---- set_algebra fields that nothing read -------------------------------------------------------

def test_A_is_subset_of_B_is_COMPUTED_and_not_asserted():
    """The module documents that A and B are built with different limits so nesting is NOT
    guaranteed. `C_is_subset_of_B` is hardcoded True (a violation raises instead), so this field is
    the only one reporting a genuine observation -- it must be able to come out False."""
    assert INV.set_algebra(["AAA"], ["AAA", "CCC"], ["AAA"])["A_is_subset_of_B"] is True
    assert INV.set_algebra(["AAA", "GGG"], ["AAA"], [])["A_is_subset_of_B"] is False


def test_the_filters_cost_is_None_rather_than_a_crash_when_A_is_empty():
    got = INV.set_algebra([], [], [])
    assert got["C_over_A_the_filters_measured_cost"] is None
    assert got["A_sequence_withheld"] == 0


# ---- the drift guard is WIRED INTO main, not just defined --------------------------------------

def test_main_REFUSES_with_exit_3_when_set_A_drifts_and_writes_no_artifact(tmp_path, monkeypatch):
    """`assert_expected_a` is unit-tested, but nothing proved `main` consults it before scoring.
    A guard that exists and is never called is the vacuous shape this repo has paid for. No real
    data is touched: the inputs are stubbed and the refusal must come before any set is built."""
    fake_mg = tmp_path / "GCF_fake"
    fake_mg.mkdir()
    (fake_mg / "annotations.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    monkeypatch.setattr(INV, "MG1655", fake_mg)
    monkeypatch.setattr(INV, "load_corpus_windows", lambda *a, **k: ["ACGT" * 10])
    monkeypatch.setattr(INV, "falsifier_eval_windows", lambda *a, **k: {"ACGT" * 10})

    def exploding(*a, **k):  # pragma: no cover - must never be reached
        raise AssertionError("main built the candidate sets despite a drifted set A")

    monkeypatch.setattr(INV, "mg1655_candidate_sets", exploding)
    out = tmp_path / "must_not_exist.json"
    assert INV.main(["--out", str(out)]) == 3
    assert not out.exists()


def test_main_REFUSES_with_exit_2_when_the_mg1655_annotation_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(INV, "MG1655", tmp_path / "nothing_here")
    out = tmp_path / "must_not_exist.json"
    assert INV.main(["--out", str(out)]) == 2
    assert not out.exists()


def test_the_default_out_is_a_committed_artifact_so_tests_must_redirect():
    """Pins the hazard rather than the fix: `main` default --out overwrites the committed inventory
    artifact, whose numbers the tests above pin. Every test here passes --out."""
    assert INV.OUT == ARTIFACT


# ---- the corpus reader and the chain scorer, both shared with the scorer ------------------------

def test_load_corpus_windows_skips_comments_and_blanks_and_uppercases(tmp_path):
    """Documented as 'the ONE way every consumer reads this file': a probe that strips differently
    sees a different `held` set, and `held` is the basis of every candidate set. Untested until now
    except indirectly through the real corpus."""
    from scripts.glm_beyond_kmer import load_corpus_windows
    p = tmp_path / "corpus.txt"
    p.write_text("# a header comment\nacgt\n\n  GGGG  \n#another\ntTtT\n", encoding="utf-8")
    assert load_corpus_windows(p) == ["ACGT", "GGGG", "TTTT"]


def test_markov_nll_reads_the_fitted_counts_on_a_real_chain():
    """Exact arithmetic on a real `MarkovGenerator`: 'AAAA' fitted once gives p(A|AAA) = 2/5 under
    add-one smoothing, so the cost is -ln(0.4)."""
    from dna_decode.glm.baselines import MarkovGenerator
    from scripts.glm_beyond_kmer import markov_nll
    chain = MarkovGenerator(order=3, seed=0).fit(["AAAA"])
    assert chain._counts["AAA"] == [2, 1, 1, 1]
    assert markov_nll(chain, "AAAA", 3) == pytest.approx(-math.log(0.4))
    # the four bases at one context form a distribution, so the implied probabilities sum to 1
    total = sum(math.exp(-markov_nll(chain, "AAA" + b, 3)) for b in "ACGT")
    assert total == pytest.approx(1.0)


def test_markov_nll_does_not_consume_the_chain_own_rng():
    """It must read `_counts`, never sample: `_next` would draw a base and measure sampling noise
    instead of the assigned probability. A consumed RNG is the observable signature of that."""
    from dna_decode.glm.baselines import MarkovGenerator
    from scripts.glm_beyond_kmer import markov_nll
    chain = MarkovGenerator(order=3, seed=0).fit(["ACGTACGTACGT"])
    before = chain.rng.getstate()
    for p in range(3, 12):
        markov_nll(chain, "ACGTACGTACGT", p)
    assert chain.rng.getstate() == before


def test_mg1655_candidate_sets_applies_the_symbol_filter_to_the_HELD_list(monkeypatch):
    """Two existing tests re-implement this filter inline, so neither would notice the real function
    changing. This calls it, with the three data-reading dependencies stubbed.

    The failure it guards: `of` must be built from `held` (corpus-absent windows), not from all
    windows -- otherwise a window the model trained on could enter the eval set."""
    from scripts import glm_beyond_kmer as gbk

    wins = [_W("AAAACCCC", "g_only"), _W("GGGGTTTT", "g_shared"),
            _W("CCCCAAAA", "g_only"), _W("TTTTGGGG", "g_unmapped")]
    monkeypatch.setattr(gbk, "build_natural_set", lambda *a, **k: (wins, None))
    monkeypatch.setattr(gbk, "symbols",
                        lambda p: {"thrL", "prophageX"} if "GCF_000005845.2" in str(p)
                        else {"thrL"})
    monkeypatch.setattr(gbk, "id_to_symbol",
                        lambda p: {"g_only": "prophageX", "g_shared": "thrL"})

    # "CCCCAAAA" is IN the corpus, so it must be excluded from `held` and therefore from `of`.
    got = gbk.mg1655_candidate_sets({"CCCCAAAA"}, length=8)
    assert got["only"] == {"prophageX"}
    assert [w.seq for w in got["held"]] == ["AAAACCCC", "GGGGTTTT", "TTTTGGGG"]
    # only the MG1655-only symbol survives; the shared and the unmapped gene are both dropped
    assert got["orthology_free"] == ["AAAACCCC"]


def test_ecoli_fit_windows_keeps_only_windows_that_are_IN_the_corpus(monkeypatch):
    """The strong null must never be fit on a window it will score. The filter direction is the
    opposite of `held`, so an inverted comparison would silently fit the null on the eval set."""
    from scripts import glm_beyond_kmer as gbk
    wins = [_W("aaaacccc", "g1"), _W("GGGGTTTT", "g2")]
    monkeypatch.setattr(gbk, "build_natural_set", lambda *a, **k: (wins, None))
    # uppercased before the membership test, so a lowercase draft window still matches
    assert gbk.ecoli_fit_windows({"AAAACCCC"}, length=8) == ["AAAACCCC"]
    assert gbk.ecoli_fit_windows(set(), length=8) == []
