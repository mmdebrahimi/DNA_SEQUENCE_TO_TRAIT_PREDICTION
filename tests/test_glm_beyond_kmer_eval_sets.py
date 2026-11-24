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

def test_the_re_pretrain_arm_is_declared_and_blocks_its_criterion():
    """R6: a declared-but-unrun expensive arm BLOCKS its criterion instead of evaporating when a
    cheap result weakens its motivation."""
    led = json.loads((ROOT / "soraya_runs" / "arms-ledger.json").read_text(encoding="utf-8"))
    arms = {(c, n): v for c, m in led["criteria"].items() for n, v in m.items()}
    planned = {k: v for k, v in arms.items() if v.get("status") == "planned"}
    assert len(planned) == 1, f"expected exactly one planned arm, got {sorted(planned)}"
    (crit, name), arm = next(iter(planned.items()))
    assert name == "representative-re-pretrain-excluding-draft-ecoli"
    assert arm["cost_tier"] == "large"
    assert "order-3 Markov chain" in crit
