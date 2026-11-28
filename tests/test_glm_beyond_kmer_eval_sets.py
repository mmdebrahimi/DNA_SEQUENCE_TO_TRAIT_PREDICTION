"""Tests for the parameterised scorer and the two-set grading.

The load-bearing test is `test_the_default_path_still_reproduces_the_committed_numbers`. The scorer
produced a committed artifact, so it is the one file here where a behaviour change IS a measurement
change -- the gate was refactored out of it and the window selection was parameterised, and both
had to leave the committed path untouched.

Second in weight is `test_positions_are_drawn_ONCE_per_window`. The scorer's whole pairing rests on
the same sampled positions serving every arm and both nulls within a window; if someone moved the
draw inside the arm loop, every delta would silently become a comparison of different positions and
nothing would error.
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.glm_beyond_kmer import (  # noqa: E402
    EVAL_SETS,
    SET_LABELS,
    SET_STRAIN_NOTES,
    EvalSetUnavailable,
    select_eval_windows,
)

CKPT = ROOT / "data" / "glm_diffusion_ckpt" / "latest.pt"
MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
COMMITTED = ROOT / "wiki" / "glm_beyond_kmer_2026-11-22.json"
TWO_SETS = ROOT / "wiki" / "glm_beyond_kmer_two_sets_2026-11-23.json"
RULE = ROOT / "wiki" / "glm_two_set_reading_rule_2026-11-24.json"

needs_model = pytest.mark.skipif(
    not (CKPT.exists() and (MG1655 / "annotations.gff3").exists() and CORPUS.exists()),
    reason="checkpoint, MG1655 cache or corpus absent on this host")


# ---- the eval-set selector, pure and offline ---------------------------------------------------

OF = [f"SEQ{i:03d}" for i in range(40)] + ["SEQ000"]      # 41 entries, 40 distinct, 1 duplicate
WITHHELD = {f"SEQ{i:03d}" for i in range(20, 60)}          # overlaps OF on SEQ020..SEQ039


def test_the_default_set_keeps_the_committed_coordinate_order_slice():
    """`of[:n]` is the committed selection and the reproduction pin depends on it not moving."""
    seqs, prov = select_eval_windows("orthology_free", of_seqs=OF, withheld_seqs=WITHHELD,
                                     n_windows=6, selection_seed=0)
    assert seqs == OF[:6]
    assert prov["deduplicated"] is False
    assert prov["selection_seed"] is None
    assert prov["n_available"] == len(OF)


def test_every_other_set_is_a_seeded_sample_of_distinct_sequences():
    seqs, prov = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                                     n_windows=10, selection_seed=0)
    assert len(seqs) == 10
    assert len(set(seqs)) == 10, "a sampled set must not repeat a window"
    assert prov["deduplicated"] is True
    assert prov["selection_seed"] == 0
    assert prov["n_available"] == len(WITHHELD)


def test_the_sample_is_deterministic_for_a_given_seed():
    a = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=10, selection_seed=0)[0]
    b = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=10, selection_seed=0)[0]
    assert a == b


def test_a_different_seed_gives_a_different_sample():
    """Otherwise the seed is decoration and the artifact's selection_seed means nothing."""
    a = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=10, selection_seed=0)[0]
    b = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=10, selection_seed=7)[0]
    assert a != b


def test_the_sample_is_NOT_the_lexicographic_prefix():
    """A sorted prefix over-selects A-rich windows, so it is a biased sample of DNA. This pins that
    the selection is a random draw rather than a head of the sorted order."""
    seqs = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                               n_windows=10, selection_seed=0)[0]
    assert seqs != sorted(WITHHELD)[:10]


def test_ortholog_bearing_is_withheld_minus_orthology_free():
    seqs = select_eval_windows("ortholog_bearing", of_seqs=OF, withheld_seqs=WITHHELD,
                               n_windows=5, selection_seed=0)[0]
    expected_pool = {s.upper() for s in WITHHELD} - {s.upper() for s in OF}
    assert set(seqs) <= expected_pool
    # and it must be DISJOINT from the orthology-free arm, or the contrast is diluted
    assert not (set(seqs) & {s.upper() for s in OF})


def test_withheld_orthology_free_is_the_intersection():
    seqs, prov = select_eval_windows("withheld_orthology_free", of_seqs=OF, withheld_seqs=WITHHELD,
                                     n_windows=5, selection_seed=0)
    expected_pool = {s.upper() for s in WITHHELD} & {s.upper() for s in OF}
    assert set(seqs) <= expected_pool
    assert prov["n_available"] == len(expected_pool)


def test_the_two_contrast_arms_are_disjoint_by_construction():
    """The whole reason the contrast is A\\C vs A&C rather than A vs C."""
    bearing = set(select_eval_windows("ortholog_bearing", of_seqs=OF, withheld_seqs=WITHHELD,
                                      n_windows=5, selection_seed=0)[0])
    free = set(select_eval_windows("withheld_orthology_free", of_seqs=OF, withheld_seqs=WITHHELD,
                                   n_windows=5, selection_seed=0)[0])
    assert not (bearing & free)


def test_an_unavailable_set_REFUSES_rather_than_returning_fewer():
    with pytest.raises(EvalSetUnavailable) as e:
        select_eval_windows("withheld_orthology_free", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=10_000, selection_seed=0)
    assert "10000" in str(e.value)


def test_an_unknown_set_name_refuses():
    with pytest.raises(EvalSetUnavailable):
        select_eval_windows("not_a_set", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=1, selection_seed=0)


# ---- the per-set notes must describe the set they ship with ------------------------------------

def test_every_eval_set_has_a_label_and_a_strain_note():
    """A missing key would KeyError at artifact-build time, after the expensive scoring run."""
    for s in EVAL_SETS:
        assert s in SET_LABELS and SET_LABELS[s]
        assert s in SET_STRAIN_NOTES and SET_STRAIN_NOTES[s]


def test_the_ortholog_bearing_note_does_not_claim_the_opposite_of_what_it_is():
    """The committed artifact hardcoded 'eval genes have NO draft counterpart'. That is FALSE for
    the ortholog-bearing arm, whose every window has one by construction."""
    note = SET_STRAIN_NOTES["ortholog_bearing"]
    assert "NO draft counterpart" not in note
    assert "DOES" in note and "near-ortholog" in note


def test_the_sequence_withheld_note_states_that_it_is_not_orthology_free():
    note = SET_STRAIN_NOTES["sequence_withheld"]
    assert "NOT orthology-free" in note


# ---- real runs ----------------------------------------------------------------------------------

@needs_model
def test_the_default_path_still_reproduces_the_committed_numbers(tmp_path):
    """THE behaviour-preservation pin. The gate was extracted and selection parameterised; the
    committed path's numbers and verdict must not have moved."""
    from scripts import glm_beyond_kmer as gbk
    out = tmp_path / "repro.json"
    argv = ["glm_beyond_kmer.py", "--n-windows", "6", "--positions", "8", "--out", str(out)]
    old = sys.argv
    try:
        sys.argv = argv
        assert gbk.main() == 0
    finally:
        sys.argv = old

    got = json.loads(out.read_text(encoding="utf-8"))
    pri = got["cells"]["left_only|markov_ecoli"]
    committed = json.loads(COMMITTED.read_text(encoding="utf-8"))
    cpri = committed["cells"]["left_only|markov_ecoli"]

    assert pri["diffusion_ce_nats"] == cpri["diffusion_ce_nats"] == 1.4098
    assert pri["null_ce_nats"] == cpri["null_ce_nats"] == 1.4291
    assert pri["paired_delta_nats"] == cpri["paired_delta_nats"] == 0.0192
    assert got["verdict"] == committed["verdict"] == "EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM"
    assert got["verdict_detail"] == committed["verdict_detail"]
    # and the extracted gate agrees with the inline original's decision
    assert got["admissibility"]["state"] == "REFUSED_BOTH_BELOW_UNIFORM"


@needs_model
def test_positions_are_drawn_ONCE_per_window(tmp_path, monkeypatch):
    """The pairing identity. The same sampled positions must serve every arm and both nulls inside
    a window, so the draw happens exactly once per window -- not once per arm."""
    from scripts import glm_beyond_kmer as gbk
    calls = []
    real = random.Random.sample

    def counting_sample(self, population, k):
        calls.append(k)
        return real(self, population, k)

    monkeypatch.setattr(random.Random, "sample", counting_sample)
    out = tmp_path / "pairing.json"
    old = sys.argv
    try:
        sys.argv = ["glm_beyond_kmer.py", "--n-windows", "4", "--positions", "5",
                    "--out", str(out)]
        assert gbk.main() == 0
    finally:
        sys.argv = old

    # The default set is NOT sampled (it is the coordinate-order slice), so every `sample` call is
    # a position draw: exactly one per window, each asking for `--positions` positions.
    assert calls == [5, 5, 5, 5], f"expected one 5-position draw per window, got {calls}"


@needs_model
def test_the_withheld_path_uses_falsifier_eval_windows_and_not_the_symbol_filter(tmp_path,
                                                                                monkeypatch):
    """`--eval-set sequence_withheld` must draw from the DESIGNED withhold, not from the
    symbol-filtered orthology-free list."""
    import dna_decode.glm.training_corpus as tc
    seen = {"called": False}
    real = tc.falsifier_eval_windows

    def spy(*args, **kwargs):
        seen["called"] = True
        return real(*args, **kwargs)

    monkeypatch.setattr(tc, "falsifier_eval_windows", spy)
    from scripts import glm_beyond_kmer as gbk
    out = tmp_path / "withheld.json"
    old = sys.argv
    try:
        sys.argv = ["glm_beyond_kmer.py", "--eval-set", "sequence_withheld",
                    "--n-windows", "4", "--positions", "5", "--out", str(out)]
        assert gbk.main() == 0
    finally:
        sys.argv = old
    assert seen["called"], "the withheld path did not consult falsifier_eval_windows"
    got = json.loads(out.read_text(encoding="utf-8"))
    assert got["eval"]["set_key"] == "sequence_withheld"
    assert got["eval"]["provenance"]["n_available"] == 3349


@needs_model
def test_a_run_does_not_touch_the_checkpoint(tmp_path):
    """The single irreplaceable artifact: ~3 h of compute, gitignored, one copy."""
    before = (CKPT.stat().st_size, CKPT.stat().st_mtime_ns)
    from scripts import glm_beyond_kmer as gbk
    old = sys.argv
    try:
        sys.argv = ["glm_beyond_kmer.py", "--n-windows", "3", "--positions", "4",
                    "--out", str(tmp_path / "x.json")]
        assert gbk.main() == 0
    finally:
        sys.argv = old
    assert (CKPT.stat().st_size, CKPT.stat().st_mtime_ns) == before
    assert before[0] == 9821617


@needs_model
def test_the_default_out_is_the_committed_artifact_so_tests_must_redirect():
    """Pins the hazard rather than the fix: the default --out OVERWRITES the committed artifact, so
    every test above passes --out. If the default ever changes, this says so."""
    from scripts.glm_beyond_kmer import OUT
    assert OUT == COMMITTED


# ---- the two-set grading ------------------------------------------------------------------------

def test_grade_contrast_branches_are_all_reachable():
    from scripts.glm_two_set_contrast import grade_contrast
    m = 0.01
    assert grade_contrast([0.02, 0.03], m) == "LEAK_LARGE_FILTER_LOAD_BEARING"
    assert grade_contrast([-0.03, -0.02], m) == "LEAK_INVERTED_FILTER_SELECTS_EASIER_WINDOWS"
    assert grade_contrast([-0.005, 0.005], m) == "LEAK_SMALL_FILTER_OVERPRICED"
    assert grade_contrast([-0.03, -0.002], m) == "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN"


def test_a_ci_touching_the_margin_is_unresolved_not_small():
    """Boundary: `LEAK_SMALL` requires the interval INSIDE the margin on both sides."""
    from scripts.glm_two_set_contrast import grade_contrast
    assert grade_contrast([-0.011, 0.005], 0.01) == "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN"
    assert grade_contrast([-0.01, 0.01], 0.01) == "LEAK_SMALL_FILTER_OVERPRICED"


def test_the_two_sample_bootstrap_is_deterministic_and_refuses_an_empty_arm():
    from scripts.glm_two_set_contrast import ContrastRefusal, two_sample_bootstrap
    x, y = [0.1, 0.2, 0.3], [0.0, 0.05, 0.1]
    a = two_sample_bootstrap(x, y, b=200, seed=0)
    b = two_sample_bootstrap(x, y, b=200, seed=0)
    assert a == b
    assert a["observed"] == pytest.approx(0.15, abs=1e-9)
    with pytest.raises(ContrastRefusal):
        two_sample_bootstrap([], y, b=10, seed=0)


def test_paired_deltas_refuses_an_artifact_without_per_window_rows():
    from scripts.glm_two_set_contrast import ContrastRefusal, paired_deltas
    with pytest.raises(ContrastRefusal):
        paired_deltas({"cells": {}})


@pytest.mark.skipif(not TWO_SETS.exists(), reason="two-set artifact not yet generated")
def test_the_two_set_artifact_parses_and_its_prose_numbers_are_in_structured_fields():
    art = json.loads(TWO_SETS.read_text(encoding="utf-8"))
    q1 = art["QUESTION_1_DOES_THE_COMMITTED_REFUSAL_STAND"]
    q2 = art["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]
    assert q1["branch"] == "ADMISSIBLE"
    assert q1["refusal_stands"] is False
    assert q1["all_four_sets_admissible"] is True
    # the sign flip, which the prose leans on
    assert q1["the_committed_cell_for_comparison"]["paired_delta_nats"] == 0.0192
    assert q1["per_set"]["orthology_free"]["paired_delta_nats"] == -0.0191
    assert q1["per_set"]["orthology_free"]["n_predictions"] == 3600
    assert q1["the_committed_cell_for_comparison"]["n_predictions"] == 48
    # the contrast
    assert q2["equivalence_margin_nats"] == 0.01
    assert q2["branch"] == "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN"
    # and the decomposition identity must hold
    d = q2["decomposition"]
    assert (d["delta_null_bearing_minus_free"] - d["delta_model_bearing_minus_free"]) == \
        pytest.approx(q2["two_sample_bootstrap"]["observed"], abs=1e-4)


@pytest.mark.skipif(not TWO_SETS.exists(), reason="two-set artifact not yet generated")
def test_the_artifact_does_not_claim_a_beyond_kmer_verdict():
    art = json.loads(TWO_SETS.read_text(encoding="utf-8"))
    blob = json.dumps(art).upper()
    assert "BEYOND_KMER_CONSERVATIVE_POSITIVE" not in blob
    assert any("NOT a beyond-k-mer verdict" in h for h in art["honest_limits"])


@pytest.mark.skipif(not TWO_SETS.exists(), reason="two-set artifact not yet generated")
def test_the_artifact_records_that_it_does_not_retire_the_expensive_arm():
    art = json.loads(TWO_SETS.read_text(encoding="utf-8"))
    assert "authority" in art["what_this_does_not_retire"].lower()


def test_the_reading_rule_was_frozen_before_any_model_number():
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    assert rule["frozen_before_any_model_number_at_this_geometry_existed"] is True
    assert rule["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]["equivalence_margin_nats"] == 0.01
    # four branches each, including an explicit unresolved one
    assert len(rule["QUESTION_1_DOES_THE_COMMITTED_REFUSAL_STAND"]["branches"]) == 4
    q2b = rule["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]["branches"]
    assert len(q2b) == 4
    assert "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN" in q2b


# ---- the declared arm must actually bind --------------------------------------------------------

def test_the_generator_criterion_does_NOT_squat_the_reserved_ceiling_prefix():
    """REGRESSION GUARD, written because this exact mistake was made here.

    `scripts/declare_successor_arms.py::ceiling_criterion` requires EXACTLY ONE criterion whose key
    starts with `ceiling:` and raises otherwise. Declaring the generator's criterion under that
    prefix made two, which took down 125 tests across 8 arms-ledger files at once -- none of which
    this plan touched, and none of which a targeted test run would have exercised.

    The `ceiling:` namespace belongs to the promoter-activity portfolio. A criterion about a
    different question gets a different prefix.
    """
    from scripts.declare_successor_arms import CRITERION_PREFIX, ceiling_criterion
    led = json.loads((ROOT / "soraya_runs" / "arms-ledger.json").read_text(encoding="utf-8"))
    reserved = [k for k in led["criteria"] if k.startswith(CRITERION_PREFIX)]
    assert len(reserved) == 1, f"the reserved {CRITERION_PREFIX!r} namespace holds {reserved}"
    # and the consumer that enforces it must resolve without raising
    assert ceiling_criterion(led) == reserved[0]
    # the generator criterion exists, and is NOT in that namespace
    gen = [k for k in led["criteria"] if "order-3 Markov chain" in k]
    assert len(gen) == 1 and not gen[0].startswith(CRITERION_PREFIX)
    # NON-VACUITY: the consumer must actually reject a second reserved criterion, otherwise this
    # whole test is asserting a condition nothing enforces.
    from scripts.declare_successor_arms import DeclareError
    squatted = {"criteria": dict(led["criteria"], **{f"{CRITERION_PREFIX} a second one": {}})}
    with pytest.raises(DeclareError):
        ceiling_criterion(squatted)


RE_PRETRAIN_ARM = "representative-re-pretrain-excluding-draft-ecoli"


def test_the_re_pretrain_arm_is_accounted_for_and_never_evaporates():
    """R6 as an IMPLICATION, converted from an exact-state pin on 2026-11-28.

    **The property R6 rests on is that an expensive declared arm cannot EVAPORATE** when a cheap
    result weakens its motivation -- it must still be in the ledger, either blocking the criterion
    because it has not run, or carrying a measurement because it has. The old form asserted
    `len(planned) == 1` and that the one planned arm was this one, which is that property PLUS a
    snapshot of the day the arm was declared and unrun.

    The snapshot is what went red: the arm RAN. It is now `measured` at 0.6535 with
    `wiki/glm_representative_repretrain_graded_2026-11-25.json` as its graded evidence, so
    `planned` is legitimately empty and the old assertion could only be kept green by writing
    `planned` back into the ledger -- i.e. by lying about an arm that has results on disk. That is
    the strongest possible form of R6's intent being satisfied: the arm did not evaporate, it was
    executed.

    This is the same conversion, for the same reason, that
    `tests/test_arms_ledger_stack_arm.py::test_the_verdict_is_blocked_IF_AND_ONLY_IF_an_arm_is_untried`
    already applied to the promoter-activity criterion on 2026-11-02: the implication never moved,
    only the arm set did.

    Asserted in BOTH directions, so neither a `planned` arm that already carries a metric nor a
    `measured` arm with no number can pass.
    """
    led = json.loads((ROOT / "soraya_runs" / "arms-ledger.json").read_text(encoding="utf-8"))
    by_crit = {c: m for c, m in led["criteria"].items() if RE_PRETRAIN_ARM in m}

    # ANTI-EVAPORATION: the whole point of R6. Present, exactly once, on the generator criterion.
    assert len(by_crit) == 1, (
        f"the expensive re-pretrain arm must appear exactly once; found it on {sorted(by_crit)}. "
        "An arm that has left the ledger has evaporated, which is what R6 forbids.")
    crit, arms = next(iter(by_crit.items()))
    assert "order-3 Markov chain" in crit, f"arm sits on an unexpected criterion: {crit!r}"

    arm = arms[RE_PRETRAIN_ARM]
    assert arm["cost_tier"] == "large", "the arm's cost tier is what made it R6-relevant"
    assert arm.get("status") in ("planned", "measured"), (
        f"status {arm.get('status')!r}: an expensive declared arm is either still blocking or "
        "measured, never anything else")

    if arm["status"] == "planned":
        # It has not run, so it must still block and must not be carrying a result.
        assert arm.get("metric") is None, (
            "an arm that still reads `planned` must not already carry a metric -- that combination "
            "means a measurement was recorded without advancing the status")
    else:
        # It ran, so the number and its evidence must both exist.
        assert isinstance(arm.get("metric"), (int, float)), (
            "a `measured` arm must carry a numeric metric; `measured` with no number is the "
            "mutate-then-save failure this ledger has already been bitten by")
        graded = ROOT / "wiki" / "glm_representative_repretrain_graded_2026-11-25.json"
        assert graded.exists(), (
            f"the arm reads `measured` but its graded evidence {graded.name} is absent; a status "
            "without an artifact behind it is an unverifiable claim")


def test_the_evaporation_guard_is_NON_VACUOUS():
    """Removing the arm must fail the guard above, or that guard asserts nothing.

    The check that matters is presence, and a presence check is exactly the kind that passes
    trivially if the lookup is wrong. So the removal is simulated against the real ledger shape.
    """
    led = json.loads((ROOT / "soraya_runs" / "arms-ledger.json").read_text(encoding="utf-8"))
    stripped = {c: {n: v for n, v in m.items() if n != RE_PRETRAIN_ARM}
                for c, m in led["criteria"].items()}
    assert not [c for c, m in stripped.items() if RE_PRETRAIN_ARM in m], (
        "the simulated removal did not actually remove the arm, so the guard it validates proves "
        "nothing")


# ================================================================================================
# `select_eval_windows` EDGE CASES.
#
# The tests above cover the happy paths of all four sets. These cover the boundaries where the
# selector's two branches differ from each other, which is where a between-set comparison could
# quietly stop being apples-to-apples.
# ================================================================================================

def test_requesting_exactly_the_pool_size_succeeds_on_both_branches():
    """Boundary on `len(pool) < n_windows`. Off by one here would refuse a fully-available set, or
    worse, admit a request for one window more than exists."""
    seqs, prov = select_eval_windows("orthology_free", of_seqs=OF, withheld_seqs=WITHHELD,
                                     n_windows=len(OF), selection_seed=0)
    assert seqs == OF and prov["n_scored"] == len(OF)

    seqs, prov = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                                     n_windows=len(WITHHELD), selection_seed=0)
    assert set(seqs) == {s.upper() for s in WITHHELD}
    assert prov["n_scored"] == len(WITHHELD)

    with pytest.raises(EvalSetUnavailable):
        select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                            n_windows=len(WITHHELD) + 1, selection_seed=0)


def test_the_default_path_needs_no_withheld_set_at_all():
    """The cost claim in `main`: set A costs ten natural-set builds, so it is only computed for a
    set that needs it and the default path is handed an EMPTY set. If the default branch ever read
    `withheld_seqs`, the committed path would silently acquire that cost -- or crash."""
    seqs, prov = select_eval_windows("orthology_free", of_seqs=OF, withheld_seqs=set(),
                                     n_windows=6, selection_seed=0)
    assert seqs == OF[:6]
    assert prov["n_available"] == len(OF)


def test_the_default_path_ignores_the_selection_seed_entirely():
    """It is a fixed coordinate-order slice, so the seed must not perturb it -- otherwise the
    reproduction pin would depend on a value the committed run never recorded."""
    base = select_eval_windows("orthology_free", of_seqs=OF, withheld_seqs=set(),
                               n_windows=6, selection_seed=0)[0]
    for seed in (1, 7, 99):
        assert select_eval_windows("orthology_free", of_seqs=OF, withheld_seqs=set(),
                                   n_windows=6, selection_seed=seed)[0] == base


def test_the_default_path_PRESERVES_a_duplicate_while_every_other_set_removes_it():
    """The deliberate asymmetry, and the hazard the inventory artifact records: the orthology-free
    LIST carries 14 repeated sequences, and the committed slice is kept byte-exact rather than
    silently deduplicated. `of[:6]` and `of[:120]` happen to be clean, which is the only reason
    this is latent. A slice that reaches the duplicate scores one window twice."""
    seqs, prov = select_eval_windows("orthology_free", of_seqs=OF, withheld_seqs=set(),
                                     n_windows=len(OF), selection_seed=0)
    assert len(seqs) == 41 and len(set(seqs)) == 40
    assert seqs.count("SEQ000") == 2
    assert prov["deduplicated"] is False
    # ...while a sampled set cannot repeat a window, whatever its pool looked like.
    dup_pool = list(WITHHELD) + list(WITHHELD)
    s2, p2 = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=dup_pool,
                                 n_windows=20, selection_seed=0)
    assert len(s2) == len(set(s2)) == 20 and p2["deduplicated"] is True


def test_the_default_path_passes_sequences_through_unchanged():
    """No uppercasing on the committed branch: `mg1655_candidate_sets` already uppercased them, and
    re-normalising here could only mask a problem upstream."""
    lower = ["acgt", "ggtt", "ttaa"]
    seqs = select_eval_windows("orthology_free", of_seqs=lower, withheld_seqs=set(),
                               n_windows=3, selection_seed=0)[0]
    assert seqs == lower


def test_the_sampled_branch_uppercases_so_a_case_mismatch_cannot_mis_split_the_arms():
    """The case trap the inventory module documents, in the selector that builds the contrast arms.
    Unnormalised, a lowercase withheld set against an uppercase orthology-free list makes A&C empty
    and A\\C the whole of A -- so the 'ortholog-bearing' arm would silently become all of set A and
    the 'orthology-free' arm would vanish, which reads as a set-definition finding, not a bug."""
    lower_withheld = [s.lower() for s in sorted(WITHHELD)[:10]]      # all present in OF, uppercased
    free, prov = select_eval_windows("withheld_orthology_free", of_seqs=OF,
                                     withheld_seqs=lower_withheld, n_windows=10,
                                     selection_seed=0)
    assert prov["n_available"] == 10 and set(free) == {s.upper() for s in lower_withheld}
    # and the complementary arm must therefore be EMPTY, not the whole of A
    with pytest.raises(EvalSetUnavailable):
        select_eval_windows("ortholog_bearing", of_seqs=OF, withheld_seqs=lower_withheld,
                            n_windows=1, selection_seed=0)


def test_the_sampled_branch_draws_from_the_SORTED_pool():
    """Pins reproducibility across processes. A set's iteration order varies with the hash seed, so
    sampling an unsorted set would give a different draw in a different interpreter while every
    same-process determinism test still passed."""
    import random as _r
    pool = sorted({s.upper() for s in WITHHELD})
    expected = _r.Random(0).sample(pool, 10)
    got = select_eval_windows("sequence_withheld", of_seqs=OF, withheld_seqs=WITHHELD,
                              n_windows=10, selection_seed=0)[0]
    assert got == expected


def test_the_two_contrast_arms_availabilities_sum_to_the_withheld_set():
    """A&C and A\\C partition A, so their `n_available` must add up. A gap would mean windows are
    being dropped from the contrast without anything saying so."""
    def avail(name):
        return select_eval_windows(name, of_seqs=OF, withheld_seqs=WITHHELD,
                                   n_windows=1, selection_seed=0)[1]["n_available"]
    assert avail("ortholog_bearing") + avail("withheld_orthology_free") == avail(
        "sequence_withheld")


# ================================================================================================
# `paired_deltas`: the SIGN convention, which decides the direction of every published claim.
# ================================================================================================

def _rows(*pairs):
    """Rows with DECOY values in the arms and nulls that must not be read."""
    return {"per_window": [{"left_only": a, "markov_ecoli": n,
                            "bidirectional": 9.0, "markov_pooled": -9.0, "uniform": 1.3863}
                           for a, n in pairs]}


def test_paired_deltas_is_POSITIVE_when_the_model_has_LOWER_cross_entropy():
    """Matches the scorer's `d = [mk - ce]`: positive = the model is better. An inverted sign here
    would flip the two-set artifact's whole reading while every magnitude stayed identical."""
    from scripts.glm_two_set_contrast import paired_deltas
    got = paired_deltas(_rows((1.0, 1.2)))
    assert got == [pytest.approx(0.2)]
    # and NEGATIVE when the null is better
    assert paired_deltas(_rows((1.2, 1.0))) == [pytest.approx(-0.2)]


def test_paired_deltas_reads_the_primary_cell_by_default_and_ignores_the_decoys():
    """Defaults must be the primary cell `left_only|markov_ecoli`; the decoy values above would
    produce wildly different deltas if the wrong keys were read."""
    from scripts.glm_two_set_contrast import PRIMARY_ARM, PRIMARY_NULL, paired_deltas
    assert (PRIMARY_ARM, PRIMARY_NULL) == ("left_only", "markov_ecoli")
    assert paired_deltas(_rows((1.0, 1.2))) == [pytest.approx(0.2)]


def test_paired_deltas_honours_an_explicit_arm_and_null():
    from scripts.glm_two_set_contrast import paired_deltas
    art = _rows((1.0, 1.2))
    assert paired_deltas(art, "bidirectional", "markov_pooled") == [pytest.approx(-18.0)]
    assert paired_deltas(art, "left_only", "uniform") == [pytest.approx(0.3863)]


def test_paired_deltas_returns_one_delta_per_row_in_row_order():
    """`n` and the bootstrap's resampling both depend on this being one-to-one with the rows."""
    from scripts.glm_two_set_contrast import paired_deltas
    got = paired_deltas(_rows((1.0, 1.2), (1.0, 0.9), (2.0, 2.0)))
    assert got == [pytest.approx(0.2), pytest.approx(-0.1), pytest.approx(0.0)]


def test_paired_deltas_refuses_an_EMPTY_per_window_list_too():
    """Distinct from the missing-key case: a present-but-empty list would otherwise flow on and
    make the bootstrap refuse later with a message about arms rather than about the artifact."""
    from scripts.glm_two_set_contrast import ContrastRefusal, paired_deltas
    with pytest.raises(ContrastRefusal):
        paired_deltas({"per_window": []})


def test_paired_deltas_raises_loudly_for_an_arm_the_artifact_never_scored():
    """It must not silently yield nothing or zero for an unscored arm."""
    from scripts.glm_two_set_contrast import paired_deltas
    with pytest.raises(KeyError):
        paired_deltas(_rows((1.0, 1.2)), "an_arm_that_was_never_scored")


# ================================================================================================
# `two_sample_bootstrap`: the contrast's own sampling distribution.
# ================================================================================================

def test_the_bootstrap_refuses_an_empty_SECOND_arm_as_well():
    """Only the first arm was covered, so a one-sided guard would have passed."""
    from scripts.glm_two_set_contrast import ContrastRefusal, two_sample_bootstrap
    with pytest.raises(ContrastRefusal):
        two_sample_bootstrap([0.1, 0.2], [], b=10, seed=0)


def test_swapping_the_arms_negates_the_observed_contrast():
    """Pins which arm is x. The rule reads positive as 'the model did better on the
    ortholog-BEARING arm', so bearing must be the first argument -- and `main` passes it first."""
    from scripts.glm_two_set_contrast import two_sample_bootstrap
    x, y = [0.10, 0.12, 0.08], [0.0, 0.01, -0.01]
    fwd = two_sample_bootstrap(x, y, b=100, seed=0)["observed"]
    rev = two_sample_bootstrap(y, x, b=100, seed=0)["observed"]
    assert fwd == pytest.approx(-rev)
    assert fwd > 0


def test_constant_arms_give_a_degenerate_interval_at_the_observed_difference():
    """With no within-arm variance there is nothing to resample, so the CI must collapse onto the
    point estimate rather than producing spurious width."""
    from scripts.glm_two_set_contrast import two_sample_bootstrap
    got = two_sample_bootstrap([0.05] * 8, [0.0] * 8, b=200, seed=0)
    assert got["observed"] == pytest.approx(0.05)
    assert got["ci95"] == [pytest.approx(0.05), pytest.approx(0.05)]


def test_the_interval_is_ordered_and_brackets_the_observed_difference():
    from scripts.glm_two_set_contrast import two_sample_bootstrap
    got = two_sample_bootstrap([0.3, 0.1, 0.2, 0.25], [0.0, 0.05, -0.05, 0.02], b=500, seed=0)
    lo, hi = got["ci95"]
    assert lo <= got["observed"] <= hi


def test_the_bootstrap_echoes_its_own_geometry_so_the_artifact_is_self_describing():
    from scripts.glm_two_set_contrast import two_sample_bootstrap
    got = two_sample_bootstrap([0.1, 0.2, 0.3], [0.0, 0.1], b=321, seed=5)
    assert (got["b"], got["seed"], got["n_x"], got["n_y"]) == (321, 5, 3, 2)


def test_a_smaller_arm_widens_the_interval():
    """Each arm is resampled at its OWN n, so the CI must carry per-arm power. If both arms were
    resampled at a shared or fixed size, a tiny arm would report a falsely tight interval."""
    from scripts.glm_two_set_contrast import two_sample_bootstrap
    spread = [0.0, 0.1, 0.2, 0.3, -0.1, -0.2, 0.15, -0.15]
    wide = two_sample_bootstrap(spread[:3], [0.0] * 3, b=800, seed=0)["ci95"]
    tight = two_sample_bootstrap(spread * 8, [0.0] * 64, b=800, seed=0)["ci95"]
    assert (wide[1] - wide[0]) > (tight[1] - tight[0])


# ---- grade_contrast boundaries the existing tests leave open -----------------------------------

def test_a_ci_only_PARTLY_beyond_the_POSITIVE_margin_is_unresolved():
    """The existing unresolved case tests the negative side only. The rule says LEAK_LARGE requires
    the interval ENTIRELY above the margin."""
    from scripts.glm_two_set_contrast import grade_contrast
    assert grade_contrast([0.005, 0.03], 0.01) == "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN"


def test_a_ci_whose_endpoint_sits_exactly_ON_the_margin_is_not_resolved():
    """`lo == margin` is not 'entirely above' and `hi == -margin` is not 'entirely below', so both
    must fall through to unresolved. Asserted against the rule's own wording in the frozen JSON."""
    from scripts.glm_two_set_contrast import grade_contrast
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    branches = rule["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]["branches"]
    assert "entirely ABOVE" in branches["LEAK_LARGE_FILTER_LOAD_BEARING"]
    assert grade_contrast([0.01, 0.03], 0.01) == "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN"
    assert grade_contrast([-0.03, -0.01], 0.01) == "UNRESOLVED_CONTRAST_SPANS_THE_MARGIN"


# ================================================================================================
# `glm_two_set_contrast.main`: the grader's refusals and its derived fields.
#
# No test executed `main` at all -- the artifact-reading tests above stay green however the grading
# code changes. It needs no checkpoint, so synthetic run artifacts exercise it completely.
# ================================================================================================

CONTRAST_SETS = ("orthology_free", "sequence_withheld",
                 "ortholog_bearing", "withheld_orthology_free")


def _fake_run(*, state="ADMISSIBLE", admissible=True, delta=0.0, n_windows=120, positions=30,
              seed=0, selection_seed=0, rows=10):
    """A minimal run artifact carrying every field `main` reads."""
    return {
        "checkpoint_step": 6000,
        "verdict": "UNRESOLVED",
        "run_args": {"n_windows": n_windows, "positions": positions, "seed": seed,
                     "selection_seed": selection_seed},
        "admissibility": {"state": state, "admissible": admissible,
                          "uniform_ce_nats": 1.3862943611198906},
        "eval": {"n_predictions": n_windows * positions, "n_windows": n_windows,
                 "positions_per_window": positions, "provenance": {"n_available": 500}},
        "cells": {f"{arm}|{null}": {"diffusion_ce_nats": 1.30, "null_ce_nats": 1.30 + delta,
                                    "paired_delta_nats": delta, "ci95": [-0.01, 0.01],
                                    "resolved": False, "n_windows": n_windows}
                  for arm in ("left_only", "bidirectional")
                  for null in ("markov_ecoli", "markov_pooled", "uniform")},
        "per_window": [{"left_only": 1.0, "bidirectional": 1.0, "markov_ecoli": 1.0 + delta,
                        "markov_pooled": -9.0, "uniform": 1.3863} for _ in range(rows)],
    }


def _grade(tmp_path, arts, *, omit=None):
    """Lay the artifacts out and run `main`, always writing to tmp_path."""
    from scripts.glm_two_set_contrast import main
    runs = tmp_path / "runs"
    runs.mkdir(exist_ok=True)
    for s, art in arts.items():
        if s == omit:
            continue
        (runs / f"run_{s}_120x30.json").write_text(json.dumps(art), encoding="utf-8")
    repro = tmp_path / "repro.json"
    repro.write_text(json.dumps(_fake_run(state="REFUSED_BOTH_BELOW_UNIFORM", admissible=False,
                                          delta=0.0192, n_windows=6, positions=8)),
                     encoding="utf-8")
    out = tmp_path / "graded.json"
    rc = main(["--runs-dir", str(runs), "--repro", str(repro), "--out", str(out)])
    return rc, (json.loads(out.read_text(encoding="utf-8")) if out.exists() else None), out


def test_the_grader_refuses_a_missing_run_artifact_and_writes_nothing(tmp_path):
    rc, got, out = _grade(tmp_path, {s: _fake_run() for s in CONTRAST_SETS},
                          omit="ortholog_bearing")
    assert rc == 2 and got is None and not out.exists()


def test_the_grader_refuses_when_the_SEED_differs_across_sets(tmp_path):
    """The guard that keeps a between-set difference from being a geometry difference. The seed is
    part of the compared tuple, so a re-run of one arm under a different seed must not be graded."""
    arts = {s: _fake_run() for s in CONTRAST_SETS}
    arts["sequence_withheld"] = _fake_run(seed=9)
    rc, got, out = _grade(tmp_path, arts)
    assert rc == 3 and not out.exists()


def test_the_grader_refuses_when_the_SELECTION_seed_differs_across_sets(tmp_path):
    """Different window draws per arm would make the contrast partly a sampling difference."""
    arts = {s: _fake_run() for s in CONTRAST_SETS}
    arts["ortholog_bearing"] = _fake_run(selection_seed=3)
    rc, _, out = _grade(tmp_path, arts)
    assert rc == 3 and not out.exists()


def test_the_grader_refuses_a_geometry_that_is_not_the_frozen_120x30(tmp_path):
    rc, _, out = _grade(tmp_path, {s: _fake_run(n_windows=60) for s in CONTRAST_SETS})
    assert rc == 3 and not out.exists()


def test_the_grader_grades_a_LARGE_leak_when_the_bearing_arm_does_better(tmp_path):
    """End-to-end through paired_deltas -> two_sample_bootstrap -> grade_contrast, which nothing
    composed before. Bearing better by 0.05 is far beyond the 0.01 margin."""
    arts = {s: _fake_run() for s in CONTRAST_SETS}
    arts["ortholog_bearing"] = _fake_run(delta=0.05)
    arts["withheld_orthology_free"] = _fake_run(delta=0.0)
    rc, got, _ = _grade(tmp_path, arts)
    assert rc == 0
    q2 = got["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]
    assert q2["branch"] == "LEAK_LARGE_FILTER_LOAD_BEARING"
    assert q2["two_sample_bootstrap"]["observed"] == pytest.approx(0.05)
    assert q2["ortholog_bearing_mean_delta"] == pytest.approx(0.05)
    assert q2["orthology_free_mean_delta"] == pytest.approx(0.0)


def test_the_grader_grades_an_INVERTED_leak_when_the_free_arm_does_better(tmp_path):
    """The mirror, so the test above cannot pass on a sign-insensitive grader."""
    arts = {s: _fake_run() for s in CONTRAST_SETS}
    arts["ortholog_bearing"] = _fake_run(delta=0.0)
    arts["withheld_orthology_free"] = _fake_run(delta=0.05)
    rc, got, _ = _grade(tmp_path, arts)
    assert rc == 0
    q2 = got["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]
    assert q2["branch"] == "LEAK_INVERTED_FILTER_SELECTS_EASIER_WINDOWS"
    assert q2["two_sample_bootstrap"]["observed"] == pytest.approx(-0.05)


def test_the_Q1_fields_are_DERIVED_from_the_gate_states_and_not_hardcoded(tmp_path):
    """`refusal_stands`, `branch` and `all_four_sets_admissible` have only ever been observed at
    one value each -- all-admissible -- so nothing proved they could come out the other way. Here
    the set the committed run actually scored REFUSES and the other three admit."""
    arts = {s: _fake_run() for s in CONTRAST_SETS}
    arts["orthology_free"] = _fake_run(state="REFUSED_MODEL_BELOW_UNIFORM", admissible=False)
    rc, got, _ = _grade(tmp_path, arts)
    assert rc == 0
    q1 = got["QUESTION_1_DOES_THE_COMMITTED_REFUSAL_STAND"]
    assert q1["refusal_stands"] is True
    assert q1["all_four_sets_admissible"] is False
    assert q1["branch"] == "REFUSED_MODEL_BELOW_UNIFORM"
    # and that branch must be one the frozen rule actually anticipates
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    assert q1["branch"] in rule["QUESTION_1_DOES_THE_COMMITTED_REFUSAL_STAND"]["branches"]


def test_a_refusal_on_a_NON_primary_set_leaves_the_refusal_not_standing(tmp_path):
    """Q1 is keyed on `orthology_free` specifically, because that is the set the committed run
    scored. A refusal elsewhere lowers `all_four_sets_admissible` without reviving the refusal."""
    arts = {s: _fake_run() for s in CONTRAST_SETS}
    arts["sequence_withheld"] = _fake_run(state="REFUSED_MODEL_BELOW_UNIFORM", admissible=False)
    rc, got, _ = _grade(tmp_path, arts)
    q1 = got["QUESTION_1_DOES_THE_COMMITTED_REFUSAL_STAND"]
    assert rc == 0
    assert q1["refusal_stands"] is False
    assert q1["all_four_sets_admissible"] is False
    assert q1["branch"] == "ADMISSIBLE"


def test_the_graded_artifact_carries_the_margin_from_the_frozen_rule(tmp_path):
    """The margin must be READ from the frozen rule, never restated in the grader."""
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    margin = rule["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]["equivalence_margin_nats"]
    _, got, _ = _grade(tmp_path, {s: _fake_run() for s in CONTRAST_SETS})
    assert got["QUESTION_2_HOW_BIG_IS_THE_MODEL_SIDE_LEAK"]["equivalence_margin_nats"] == margin


def test_the_grader_default_out_is_a_committed_artifact_so_tests_must_redirect():
    """Same hazard as the scorer's: the default --out overwrites the committed two-set artifact
    whose numbers the tests above pin. Every test here passes --out."""
    from scripts.glm_two_set_contrast import OUT
    assert OUT == TWO_SETS


# ================================================================================================
# The persisted per-window rows, which exist ONLY so the between-set contrast is computable later.
# ================================================================================================

@needs_model
def test_the_scorer_persists_per_window_rows_that_reconcile_with_every_cell(tmp_path):
    """The rows are the whole reason the scorer was modified: the committed artifact kept only cell
    summaries, which made a two-sample bootstrap on a difference of paired deltas impossible after
    the fact. Nothing checked that what is persisted actually reproduces the cells it ships with --
    so a row written from the wrong arm would be invisible until a grader quietly mis-graded."""
    import statistics

    from scripts import glm_beyond_kmer as gbk
    from scripts.glm_two_set_contrast import paired_deltas

    out = tmp_path / "rows.json"
    old = sys.argv
    try:
        sys.argv = ["glm_beyond_kmer.py", "--n-windows", "3", "--positions", "4",
                    "--boot", "50", "--out", str(out)]
        assert gbk.main() == 0
    finally:
        sys.argv = old

    art = json.loads(out.read_text(encoding="utf-8"))
    rows = art["per_window"]
    assert len(rows) == 3
    for row in rows:
        assert set(row) == {"left_only", "bidirectional", "markov_pooled", "markov_ecoli",
                            "uniform"}
    # every cell's paired delta must be recoverable from the rows alone
    for key, cell in art["cells"].items():
        arm, null = key.split("|")
        assert statistics.fmean(paired_deltas(art, arm, null)) == pytest.approx(
            cell["paired_delta_nats"], abs=1e-4), key
    # and the run is self-describing, which the committed artifact was not
    # Exact-shape equality is kept deliberately: it is what catches a field silently vanishing.
    # `ckpt` is held out of the comparison only because it is an absolute, machine-specific path,
    # and is then asserted on its own rather than dropped.
    run_args = dict(art["run_args"])
    ckpt = run_args.pop("ckpt")
    assert run_args == {"n_windows": 3, "positions": 4, "seed": 0, "boot": 50,
                        "eval_set": "orthology_free", "selection_seed": 0,
                        "ckpt_is_committed": True,
                        "eval_sets_and_nulls_pinned_to": "glm_training_corpus_2026-11-21.txt"}
    assert Path(ckpt).resolve() == gbk.CKPT.resolve()
    assert art["eval"]["provenance"]["set"] == "orthology_free"
