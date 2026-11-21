"""The per-seed margin must be computed honestly, and the probe must refuse when it cannot anchor.

Three things are pinned here, and the second is the one that matters most:

  1. `margin = null - ceiling`, NEVER clamped. A `max(0, ...)` would turn "the null sits below its own
     noise floor, so this mode has no winnable region" into "the margin is zero", which reads as a
     near-miss. That edit is the single most likely way this finding gets silently inverted.
  2. The seed-7 anchor gate is NON-VACUOUS -- perturbing the ceiling makes the probe refuse. A gate that
     cannot fail is not a gate, and this repo has shipped two of those (the HIV operator-prefixed fold
     filter and the AMRFinder POINT-row filter, both of which removed exactly nothing).
  3. The seed list is READ from the committed grid, not restated, so a grid built on different seeds
     cannot be scored against the wrong ten.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_seed_role_ceilings.py"
OUT = ROOT / "wiki" / "glm_seed_role_ceilings_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_seedceil", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the unclamped margin

def test_the_margin_is_a_plain_difference(mod):
    assert mod.margin(0.6679, 0.5865) == pytest.approx(0.0814)
    assert mod.margin(0.5009, 0.5262) == pytest.approx(-0.0253)


def test_a_ceiling_above_the_null_yields_a_NEGATIVE_margin_with_no_clamping(mod):
    """The load-bearing one. A null below its own noise floor is the FINDING, not an error."""
    m = mod.margin(0.4623, 0.5401)
    assert m < 0, "the margin was clamped; a negative margin is the finding"
    assert m == pytest.approx(-0.0778)


# ------------------------------------------------------------------ the anchor gate is real

def test_the_seed_list_is_read_from_the_grid_not_restated(mod):
    grid = {"grid": [{"seed": 3, "null": {"kmer": 0.1, "both": 0.2}},
                     {"seed": 11, "null": {"kmer": 0.3, "both": 0.4}}]}
    got = mod.seeds_and_nulls_from_grid(grid)
    assert sorted(got) == [3, 11], "the probe used its own seed list instead of the grid's"
    assert got[11] == {"kmer": 0.3, "both": 0.4}


def test_a_seed_carrying_two_different_nulls_is_REFUSED(mod):
    """Collapsing to one null per seed is exact ONLY because all models in a seed share it.

    If that ever stops holding, averaging over the difference would be silent, so it refuses instead.
    """
    grid = {"grid": [{"seed": 1, "null": {"kmer": 0.50, "both": 0.60}},
                     {"seed": 1, "null": {"kmer": 0.51, "both": 0.60}}]}
    with pytest.raises(SystemExit) as e:
        mod.seeds_and_nulls_from_grid(grid)
    assert "distinct null pairs" in str(e.value)


def test_the_committed_anchor_values_are_the_ones_the_repo_recorded(mod):
    """Pinned so a future edit cannot quietly relax the gate by moving the target."""
    assert mod.COMMITTED == {"kmer": 0.5447, "both": 0.5865}
    assert mod.ANCHOR_SEED == 7
    assert mod.GATE_TOL <= 1e-3, "a loose tolerance makes the anchor gate decorative"


def test_the_anchor_gate_is_non_vacuous(mod, monkeypatch, tmp_path):
    """Perturb the measured ceiling by 0.01 and the probe must REFUSE rather than write.

    Checked by substituting the measurement, not by editing the committed constant -- a test that moved
    the constant would prove only that the constant is read.
    """
    real = mod.ceilings_for_seed

    def perturbed(seed: int):
        out = real(seed) if seed != mod.ANCHOR_SEED else None
        if out is None:
            out = {m: {"ceiling": mod.COMMITTED[m] + 0.01, "values": [0.5] * 5, "n_eval": 600}
                   for m in mod.MODES}
        return out

    monkeypatch.setattr(mod, "ceilings_for_seed", perturbed)
    monkeypatch.setattr(mod, "OUT", tmp_path / "must-not-appear.json")
    rc = mod.main()
    assert rc == 3, f"the probe returned {rc} on a broken anchor instead of refusing"
    assert not (tmp_path / "must-not-appear.json").exists(), "it wrote an artifact despite refusing"


# ------------------------------------------------------------------ the verdict wording

def test_a_negative_mean_margin_reads_as_no_winnable_region(mod):
    s = {"margin_mean": -0.0121, "seeds_with_negative_margin": 6, "seeds_with_positive_margin": 4}
    assert mod.verdict_for_mode(s) == "NO_WINNABLE_REGION_ON_AVERAGE"


def test_a_unanimously_positive_margin_reads_as_winnable_at_every_seed(mod):
    s = {"margin_mean": 0.0964, "seeds_with_negative_margin": 0, "seeds_with_positive_margin": 10}
    assert mod.verdict_for_mode(s) == "WINNABLE_REGION_AT_EVERY_SEED"


def test_a_positive_mean_with_sign_flips_is_called_SEED_DEPENDENT_not_winnable(mod):
    """A positive average that flips sign across seeds must not read as a clean winnable region."""
    s = {"margin_mean": 0.02, "seeds_with_negative_margin": 3, "seeds_with_positive_margin": 7}
    assert mod.verdict_for_mode(s) == "WINNABLE_REGION_SEED_DEPENDENT"


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_committed_artifact_is_internally_consistent():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["anchor_gate"]["passed"] is True
    assert len(rec["rows"]) == 20, "10 seeds x 2 modes"
    for r in rec["rows"]:
        assert r["margin"] == pytest.approx(r["null"] - r["ceiling"], abs=1e-4), r
    # The two modes must NOT be reported with the same verdict by accident.
    assert rec["verdict_by_mode"]["kmer"] != rec["verdict_by_mode"]["both"], (
        "both modes returned the same verdict; the per-mode split is the whole point of this probe")


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_artifact_records_WHY_the_number_was_missing():
    """The schema gap is the reason the per-seed margin never existed, and it belongs in the record."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    why = rec["why_this_was_missing"]
    assert "no ceiling field" in why
    assert "0.5447" in why, "the single-seed ceiling that was over-generalised must be named"
