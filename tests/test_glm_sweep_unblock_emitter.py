"""`wiki/glm_generator_sampling_sweep.json` is a TRIGGER FILE, so its emitter needs guards.

WHY (2026-10-31). That exact filename is the `file-exists` expiry test for constraint
`gd-finetune-gated-on-sampling-sweep-not-compute`, and both CLAUDE.md and the 2026-10-30 handoff
say in as many words: do not create it BY HAND. Writing it on an unverified grid would fire a LoRA
fine-tune's unblock trigger on a sweep this project's own gate excludes. On 2026-10-29 the sweep
ran and the file was DELIBERATELY withheld for exactly that reason -- 2 of 3 grid points returned
`excluded:requested-but-unverified` because the kernel never passed `--verify-knob`.

So the discharge condition lives in code, and these tests pin the two refusals. The second one is
the non-obvious guard and the reason this file exists: a grid of ALL-DEFAULT points is
"all eligible" -- `eligible:default-sampling` is a genuine eligibility -- and yet it attests
NOTHING about the knob. A naive `all(eligible)` check would emit the trigger off a grid where the
temperature knob was never requested, let alone exercised. That is the same
nothing-to-verify-vs-unverified collapse `sweep_gates` was written to prevent, re-entering through
the emitter.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "emit_glm_sweep_unblock_artifact.py"
CANONICAL = ROOT / "wiki" / "glm_generator_sampling_sweep.json"

_VERIFIED = {"requested_non_default": True, "knob_wired": True,
             "knob_verified_end_to_end": True,
             "verification": {"verdict": "TEMPERATURE_IS_LIVE", "live": True}}
_UNVERIFIED = {"requested_non_default": True, "knob_wired": True,
               "knob_verified_end_to_end": False}
_DEFAULT = {"requested_non_default": False, "knob_wired": False,
            "knob_verified_end_to_end": False}


def _emitter():
    spec = importlib.util.spec_from_file_location("_emit", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _cell(tmp: Path, name: str, temperature, knob: dict) -> str:
    p = tmp / name
    p.write_text(json.dumps({
        "knob_contract": knob,
        "sampling": {"temperature": temperature, "top_k": 0,
                     "is_recorded_default": temperature == 1.0},
        "results": [],
    }), encoding="utf-8")
    return str(p)


def _summary(tmp: Path) -> str:
    p = tmp / "summary.json"
    p.write_text(json.dumps({"record": "glm-sampling-sweep-v1", "grid": []}), encoding="utf-8")
    return str(p)


def _run(m, tmp: Path, cells: list):
    return m.main(["--summary", _summary(tmp), "--cells", *cells,
                   "--source-kernel", "test", "--out", str(tmp / "trigger.json")])


def test_refuses_when_any_non_default_point_is_unverified(tmp_path):
    """The 2026-10-29 state: the knob was requested, wired, and never attested."""
    m = _emitter()
    cells = [_cell(tmp_path, "t1.json", 1.0, _DEFAULT),
             _cell(tmp_path, "t08.json", 0.8, _UNVERIFIED),
             _cell(tmp_path, "t06.json", 0.6, _UNVERIFIED)]
    with pytest.raises(SystemExit) as e:
        _run(m, tmp_path, cells)
    assert "not plan-eligible" in str(e.value)
    assert not (tmp_path / "trigger.json").exists(), "refused but still wrote the trigger"


def test_refuses_an_all_default_grid_even_though_every_point_is_eligible(tmp_path):
    """THE load-bearing guard. Every point is eligible; the grid still attests nothing.

    `eligible:default-sampling` is a real eligibility, so `all(eligible)` passes here. But no
    point ever REQUESTED a non-default knob, so nothing exercised the lever the trigger exists to
    prove is live. Emitting on this grid would discharge the constraint off a sweep that never
    swept.
    """
    m = _emitter()
    cells = [_cell(tmp_path, "a.json", 1.0, _DEFAULT),
             _cell(tmp_path, "b.json", 1.0, _DEFAULT)]
    from dna_decode.glm.sweep_gates import ELIGIBLE, verdict_eligibility
    assert all(verdict_eligibility(_DEFAULT) in ELIGIBLE for _ in cells), (
        "premise of this test is gone: a default point is no longer eligible")
    with pytest.raises(SystemExit) as e:
        _run(m, tmp_path, cells)
    assert "no grid point carries a VERIFIED non-default knob" in str(e.value)
    assert not (tmp_path / "trigger.json").exists()


def test_emits_on_a_verified_grid(tmp_path):
    m = _emitter()
    cells = [_cell(tmp_path, "t1.json", 1.0, _DEFAULT),
             _cell(tmp_path, "t08.json", 0.8, _VERIFIED),
             _cell(tmp_path, "t06.json", 0.6, _VERIFIED)]
    assert _run(m, tmp_path, cells) == 0
    out = json.loads((tmp_path / "trigger.json").read_text(encoding="utf-8"))
    att = out["knob_attestation"]
    assert att["n_points"] == 3
    assert att["n_eligible"] == 3
    assert att["n_verified_non_default"] == 2


def test_refuses_to_overwrite_without_force(tmp_path):
    m = _emitter()
    cells = [_cell(tmp_path, "t1.json", 1.0, _DEFAULT),
             _cell(tmp_path, "t08.json", 0.8, _VERIFIED)]
    assert _run(m, tmp_path, cells) == 0
    with pytest.raises(SystemExit) as e:
        _run(m, tmp_path, cells)
    assert "already exists" in str(e.value)


def test_the_committed_trigger_records_a_verified_grid():
    """The live artifact, not a fixture: it must not have been emitted off an unverified sweep."""
    if not CANONICAL.is_file():
        pytest.skip("trigger not emitted in this checkout")
    d = json.loads(CANONICAL.read_text(encoding="utf-8"))
    att = d["knob_attestation"]
    assert att["n_eligible"] == att["n_points"]
    assert att["n_verified_non_default"] >= 1
    for p in att["points"]:
        assert p["eligibility"] in ("eligible:verified", "eligible:default-sampling"), p
        if p["requested_non_default"]:
            assert p["knob_verified_end_to_end"] is True, p
    # It must not read as an unblock of G-D -- the edge is held by a user call.
    assert "AUTHORITY" in d["discharges"]["still_blocked_because"]
