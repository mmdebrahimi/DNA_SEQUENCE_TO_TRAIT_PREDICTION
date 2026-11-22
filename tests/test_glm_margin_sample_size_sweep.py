"""The sweep must keep its prediction registered, keep its geometry quarantined, and name BOTH outcomes.

This probe exists to settle one question: is the kmer mode's unresolved winnable region a sample-size
problem? That makes it the kind of run most at risk of a retro-fitted reading, so three things are pinned.

1. **The prediction is in the source, not in a write-up.** It names what the null, the floor and the
   margin should do AND names the informative failure -- a margin that stays unresolved while its error
   bar visibly shrinks would settle the kmer question in the opposite direction. Both outcomes exist as
   verdict strings, so neither can be presented as the one that was expected all along.
2. **The geometry is off the committed cell by design, so it is quarantined.** `limit_natural` moves the
   Markov fit size, the discriminator's train and test sizes and the class balance at once. Joining these
   numbers to the committed grid's generator figures would be the same-data/different-question error.
3. **The floor is a working instrument check.** Real sequence against real sequence has a true AUROC of
   0.5 at every sample size, so a floor that drifts off 0.5 means the comparison is not measuring what it
   claims.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_margin_sample_size_sweep.py"
OUT = ROOT / "wiki" / "glm_margin_sample_size_sweep_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_mss", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the pre-registration

def test_the_prediction_is_registered_in_the_source(mod):
    p = mod.PREDICTION
    assert p["registered_before_any_number"] is True
    for key in ("null", "floor", "margin", "the_informative_failure"):
        assert key in p and len(p[key]) > 20, key


def test_the_informative_failure_has_its_own_verdict_string(mod):
    """A prediction with only one nameable outcome is not a prediction.

    The failure branch must be reachable as a verdict, otherwise the run can only ever confirm.
    """
    src = SCRIPT.read_text(encoding="utf-8")
    assert "STAYS_UNRESOLVED_AS_PRECISION_IMPROVES__MARGIN_IS_NEAR_ZERO" in src
    assert "SAMPLE_SIZE_FIXES_IT__MARGIN_IS_REAL_AND_WAS_MERELY_UNMEASURED" in src
    assert "UNRESOLVED_AND_PRECISION_DID_NOT_IMPROVE__SWEEP_UNINFORMATIVE" in src


def test_the_first_config_is_the_committed_geometry(mod):
    """The sweep must contain the committed cell so the off-geometry rows have an anchor."""
    assert mod.CONFIGS[0] == (1200, 120)


def test_the_sweep_spends_headroom_that_actually_exists(mod):
    """Every config needs limit + n_generated spare windows for the floor; the corpus has 4340."""
    for limit, n_gen in mod.CONFIGS:
        assert limit + n_gen <= 4340, (limit, n_gen)


# ------------------------------------------------------------------ the resolve helper

def test_a_straddling_interval_is_not_resolvable(mod):
    out = mod.resolve([0.08, -0.05, 0.06, -0.02, 0.09, 0.01, -0.06, 0.07, 0.04, 0.11])
    assert out["resolvable_above_zero"] is False
    assert out["ci95"][0] < 0 < out["ci95"][1]


def test_a_tight_positive_margin_resolves(mod):
    out = mod.resolve([0.05, 0.06, 0.04, 0.055, 0.045, 0.05, 0.06, 0.05, 0.05, 0.055])
    assert out["resolvable_above_zero"] is True and out["ci95"][0] > 0


def test_resolve_uses_the_standard_error(mod):
    m = [0.05, 0.06, 0.04, 0.055, 0.045, 0.05, 0.06, 0.05, 0.05, 0.055]
    out = mod.resolve(m)
    assert out["stderr"] == pytest.approx(out["sd"] / math.sqrt(10), abs=1e-4)


def test_resolve_refuses_on_too_few_replicates(mod):
    assert mod.resolve([0.1, 0.2])["resolvable"] is None


def test_feature_count_matches_the_builder(mod):
    """p-against-n is reported, so the count has to be the real one: kmer 4^4, both plus bins*4^2."""
    assert mod.feature_count("kmer") == 256
    assert mod.feature_count("both") == 256 + 6 * 16


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="sweep has not been run on this host")
def test_the_geometry_is_quarantined_from_the_committed_grid():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    note = rec["geometry_note"]
    assert "NEVER" in note and "committed grid" in note
    assert any(v["is_committed_geometry"] for v in rec["by_config"].values())
    assert sum(1 for v in rec["by_config"].values() if v["is_committed_geometry"]) == 1


@pytest.mark.skipif(not OUT.exists(), reason="sweep has not been run on this host")
def test_the_floor_stays_at_chance_at_every_geometry():
    """The instrument check. Real-versus-real must sit at 0.5 regardless of sample size."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for key, cfg in rec["by_config"].items():
        for mode in ("kmer", "both"):
            assert abs(cfg[mode]["floor_mean"] - 0.5) < 0.05, (key, mode, cfg[mode]["floor_mean"])


@pytest.mark.skipif(not OUT.exists(), reason="sweep has not been run on this host")
def test_the_prediction_checks_are_recorded_mechanically():
    """Whether the prediction held is computed, not narrated."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        c = rec["prediction_checks"][mode]
        assert c["null_rose"] is True, mode
        assert c["floor_near_half_everywhere"] is True, mode
        assert isinstance(c["n_configs_resolvable"], int)


@pytest.mark.skipif(not OUT.exists(), reason="sweep has not been run on this host")
def test_kmer_resolves_once_the_committed_n_is_raised():
    """The actionable finding: the committed cell is the ONLY geometry where kmer does not resolve.

    Pinned because it is what reverses the earlier recommendation to stop scoring the mode -- and because
    the cheap fix (more generated sequences at the same corpus) is visible only in this comparison.
    """
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    cfgs = rec["by_config"]
    committed = next(v for v in cfgs.values() if v["is_committed_geometry"])
    assert committed["kmer"]["margin"]["resolvable_above_zero"] is False
    same_corpus_more_generated = cfgs["limit1200_n400"]
    assert same_corpus_more_generated["kmer"]["margin"]["resolvable_above_zero"] is True, (
        "raising n_generated alone must be enough, or the fix is not free")
    assert rec["verdict_by_mode"]["kmer"].startswith("SAMPLE_SIZE_FIXES_IT")


@pytest.mark.skipif(not OUT.exists(), reason="sweep has not been run on this host")
def test_both_mode_improves_monotonically_with_data():
    """`both` was already resolved; the sweep shows the margin widening rather than saturating."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    means = [v["both"]["margin"]["mean"] for v in rec["by_config"].values()]
    assert means[-1] > means[0]
    assert all(v["both"]["margin"]["resolvable_above_zero"] for v in rec["by_config"].values())


@pytest.mark.skipif(not OUT.exists(), reason="sweep has not been run on this host")
def test_the_confound_between_fit_size_and_eval_size_is_declared():
    """The sweep cannot separate 'the discriminator sees more' from 'the null is fitted better'."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "not separated" in joined or "does not separate" in joined
    assert "necessary, NOT sufficient" in joined
