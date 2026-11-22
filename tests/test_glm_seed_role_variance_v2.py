"""The v2 decomposition must have seven genuinely different arms, and each must move what its name says.

The whole experiment is a four-column seed table, so that is where it can silently stop measuring
anything. Two arms sharing a plan would produce two identical SDs and a decomposition that reads as
perfectly sensible; an arm called `split_only` whose split seed was accidentally pinned would report a
spread caused by something else under a name asserting otherwise. Both are pinned here directly, without
running a single discriminator.

The v1 defect is also pinned: v1 could not have caught its own geometry error because it never compared
its all-coupled arm against the grid. `test_the_gate_against_the_grid_is_not_optional` asserts the gate
exists and refuses, and `test_the_v1_geometry_is_recorded_as_the_defect` pins the n=120-against-n=400
difference in the artifact so a later tidy-up cannot quietly drop the reason v2 exists.
"""

from __future__ import annotations

import importlib.util
import json
import statistics
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_seed_role_variance_v2.py"
OUT = ROOT / "wiki" / "glm_seed_role_variance_v2_2026-11-21.json"

ROLES = ("fit_corpus", "eval_corpus", "generation", "split")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_seedvar2", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the geometry that caused v2

def test_the_committed_geometry_is_n120_not_n400(mod):
    """The single value v1 got wrong. v1 used 400; the kernel froze 120."""
    assert mod.N_GENERATED == 120
    assert mod.LIMIT == 1200
    assert mod.LENGTH == 150
    assert mod.POSITIONAL_BINS == 6
    assert mod.MARKOV_ORDER == 3


# ------------------------------------------------------------------ the arm wiring

def test_arm_A_moves_all_four_roles_together(mod):
    assert mod.seed_plan("A_all_coupled", 3) == (3, 3, 3, 3)
    assert mod.seed_plan("A_all_coupled", 9) == (9, 9, 9, 9)


def test_arm_A2_keeps_one_logical_seed_but_separates_the_rng_streams(mod):
    """The rival v1 never tested: three `random.Random` objects seeded with the SAME integer.

    A2 must still be driven by one seed per replicate (so the logical coupling is intact) while giving
    generation and split distinct seeds (so the streams are independent). If A2 ever collapsed to A, the
    stream-sharing rival would silently stop being tested while still being reported.
    """
    a = mod.seed_plan("A2_streams_decoupled", 4)
    assert a[0] == a[1] == 4, "A2 must keep the corpus at the replicate seed"
    assert a[2] != 4 and a[3] != 4, "A2 must NOT reuse the replicate seed for generation or split"
    assert len({a[1], a[2], a[3]}) == 3, "the three streams must be distinct"
    assert a != mod.seed_plan("A_all_coupled", 4), "A2 collapsed onto A"


def test_arm_F_is_the_only_arm_whose_fit_and_eval_corpora_differ(mod):
    """Arm F is the direct test of the mechanism v1 asserted and could not have measured.

    In every v1 arm the fit and eval halves came from ONE corpus call, so the fit-set matching was
    constant across arms and could not explain a difference between them.
    """
    differ = [a for a in mod.ARMS
              if any(p[0] != p[1] for p in (mod.seed_plan(a, s) for s in range(6)))]
    assert differ == ["F_fit_decoupled"], differ
    f = mod.seed_plan("F_fit_decoupled", 2)
    assert f[0] != f[1] and f[1] == 2


def test_each_single_role_arm_pins_everything_else(mod):
    seeds = list(range(8))
    assert mod.arm_varies_only("B_generation_only", seeds) == {"generation"}
    assert mod.arm_varies_only("D_corpus_only", seeds) == {"fit_corpus", "eval_corpus"}
    assert mod.arm_varies_only("E_split_only", seeds) == {"split"}


def test_arm_C_fuses_corpus_and_split_which_is_why_D_and_E_exist(mod):
    """v1's arm C moved corpus AND split, so it bounded their combined contribution only."""
    seeds = list(range(8))
    assert mod.arm_varies_only("C_corpus_and_split", seeds) == {"fit_corpus", "eval_corpus", "split"}


def test_an_unknown_arm_is_refused_not_defaulted(mod):
    with pytest.raises(ValueError):
        mod.seed_plan("Z_invented", 3)


def test_the_seven_arms_are_not_accidentally_identical(mod):
    assert mod.arms_are_distinct(list(range(10))) is True


def test_distinctness_collapses_detectably_when_there_is_nothing_to_vary(mod):
    """Non-vacuity for the guard: on a single-seed list equal to the anchor several arms coincide."""
    assert mod.arms_are_distinct([7], anchor=7) is False


def test_every_arm_name_matches_the_roles_it_actually_moves(mod):
    """The design table and the arm names must agree, or an arm reports under a false label."""
    seeds = list(range(10))
    expect = {
        "A_all_coupled": set(ROLES),
        "A2_streams_decoupled": set(ROLES),
        "B_generation_only": {"generation"},
        "C_corpus_and_split": {"fit_corpus", "eval_corpus", "split"},
        "D_corpus_only": {"fit_corpus", "eval_corpus"},
        "E_split_only": {"split"},
        "F_fit_decoupled": set(ROLES),
    }
    for arm, want in expect.items():
        assert mod.arm_varies_only(arm, seeds) == want, arm


# ------------------------------------------------------------------ the bootstrap

def test_identical_spreads_give_a_ratio_interval_containing_one(mod):
    """Control: two samples with the SAME spread must NOT be reported as differing.

    The second sample is the first shifted by a constant, because standard deviation is invariant under
    a location shift -- so the true ratio is exactly 1.0 and any `excludes_one` here would be a false
    positive. An earlier version of this test hand-wrote two 'similar looking' lists whose spreads
    actually differed by 32%; the helper correctly flagged them and the TEST was wrong, which is its own
    reminder that eyeballing two lists is not a way to equalise a variance.
    """
    xs = [0.50, 0.52, 0.48, 0.51, 0.49, 0.53, 0.47, 0.50, 0.52, 0.48] * 3
    ys = [v + 0.04 for v in xs]
    assert statistics.stdev(xs) == pytest.approx(statistics.stdev(ys), abs=1e-12)
    out = mod.bootstrap_sd_ratio(xs, ys)
    assert out["point"] == pytest.approx(1.0, abs=1e-6)
    assert out["excludes_one"] is False, out


def test_a_clearly_wider_sample_is_flagged_as_wider(mod):
    tight = [0.500, 0.501, 0.499, 0.500, 0.501, 0.499, 0.500, 0.501, 0.499, 0.500] * 3
    wide = [0.40, 0.60, 0.42, 0.58, 0.44, 0.56, 0.46, 0.54, 0.38, 0.62] * 3
    out = mod.bootstrap_sd_ratio(wide, tight)
    assert out["point"] > 5 and out["excludes_one"] is True, out


def test_the_ratio_point_estimate_is_the_real_sd_ratio(mod):
    """The bootstrap must report the observed ratio, not a resampled mean of ratios."""
    xs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    ys = [0.30, 0.31, 0.29, 0.32, 0.28, 0.33]
    out = mod.bootstrap_sd_ratio(xs, ys)
    assert out["point"] == pytest.approx(statistics.stdev(xs) / statistics.stdev(ys), abs=1e-3)


# ------------------------------------------------------------------ the gate

def test_the_gate_against_the_grid_is_not_optional(mod):
    """v1 shipped an arm claiming to reproduce the runner, with no code that checked.

    The gate must be present in the source AND must refuse rather than warn, because a warn-then-continue
    gate is what let the off-geometry numbers be read as the status quo in the first place.
    """
    src = SCRIPT.read_text(encoding="utf-8")
    assert "GATE_FAILED_NO_NUMBERS_REPORTED" in src
    assert "grid_nulls" in src
    assert mod.GATE_TOL <= 1e-3


def test_more_replicates_than_v1(mod):
    """A ratio of two SDs at n=10 cannot separate these hypotheses; v1 said so and used 10 anyway."""
    assert mod.N_REPLICATES >= 30


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="v2 has not been run on this host")
def test_the_artifact_records_the_gate_as_passed_with_per_seed_rows():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["gate"]["passed"] is True
    assert len(rec["gate"]["rows"]) >= 20
    assert all(r["reconciles"] for r in rec["gate"]["rows"])


@pytest.mark.skipif(not OUT.exists(), reason="v2 has not been run on this host")
def test_the_v1_geometry_is_recorded_as_the_defect():
    """The reason v2 exists must survive a tidy-up, numerically and not just in prose."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    why = rec["why_it_supersedes"]
    assert "400" in why and "120" in why
    assert "0.0140" in why and "0.0410" in why
    assert rec["supersedes"].endswith("glm_seed_role_variance_2026-11-21.json")


@pytest.mark.skipif(not OUT.exists(), reason="v2 has not been run on this host")
def test_every_arm_produced_more_than_one_distinct_value():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for arm, modes in rec["summary"].items():
        for mode, s in modes.items():
            assert s["n_distinct"] >= 2, (arm, mode, s["n_distinct"])


@pytest.mark.skipif(not OUT.exists(), reason="v2 has not been run on this host")
def test_every_arm_to_arm_claim_carries_an_interval():
    """A point ratio is a coarse reading; the artifact must ship the CI beside it."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode, c in rec["comparisons"].items():
        for arm, r in c["vs_A"].items():
            assert r["ci95"] is not None and len(r["ci95"]) == 2, (mode, arm)
            assert r["excludes_one"] in (True, False), (mode, arm)


@pytest.mark.skipif(not OUT.exists(), reason="v2 has not been run on this host")
def test_the_artifact_scopes_itself_to_the_markov_probe():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert "NOT the HF" in rec["probe_candidate"]
    assert any("HF generator" in lim for lim in rec["honest_limits"])


@pytest.mark.skipif(not OUT.exists(), reason="v2 has not been run on this host")
def test_the_v1_anchor_sweep_is_not_silently_inherited():
    """v1's anchor sweep ran at the wrong geometry, so citing it here would launder a stale control."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert any("NOT inherited" in lim for lim in rec["honest_limits"])
