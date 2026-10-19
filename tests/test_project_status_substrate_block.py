"""`scripts/project_status.py::_print_substrate_block` -- the renderer had NO automated test.

It was verified by running the script and reading the output, which is how the defect below survived
code review: on 2026-10-13 the block was added with `sym['detail'][:150]` and `w[:160]`, sixty lines
BELOW a comment in the same function that says, verbatim, *"WRAPPED, never truncated. A hard [:150]
cut this very block's warning mid-sentence ... silently amputating the load-bearing half of a caveat
is the same failure as omitting it."*

MEASURED consequence before the fix:
  - on the real 4-rung asymmetric census, `[:150]` dropped **size_lever** from the warning entirely
    -- the single field whose asymmetry caused the incident the block exists to prevent -- and cut a
    rung name mid-token at `'sigma|70_designed'`;
  - on the live census, `[:160]` cut all three ceiling warnings mid-word at `square-ro|oted)`,
    losing the Spearman-Brown derivation that is the whole point of the warning.

So the mechanism built to stop a rung-scoped fact being read as a project-level wall was itself
amputating the name of the rung-scoped field. The guards here are the source-level one
(`test_the_substrate_block_carries_no_hard_character_cut`, which fires even when the content happens
to be short) plus the two content guards that reproduce the real shapes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import scripts.project_status as ps  # noqa: E402
from dna_decode.eval.substrate_symmetry import SYMMETRY_REQUIRED, audit  # noqa: E402

CENSUS = ROOT / "wiki" / "mpra_substrate_census_2026-10-12.json"
SCRIPT = ROOT / "scripts" / "project_status.py"


def _render(monkeypatch, capsys, census_path: Path) -> str:
    """Render the block alone against a chosen census path."""
    monkeypatch.setattr(ps, "_CENSUS", census_path)
    ps._print_substrate_block()
    return capsys.readouterr().out


def _write(tmp_path: Path, payload) -> Path:
    p = tmp_path / "census.json"
    p.write_text(payload if isinstance(payload, str) else json.dumps(payload), encoding="utf-8")
    return p


#: The pre-amendment census shape, which is the shape the validator exists to refuse: `size_lever`
#: and `best_measured_arm` on peak_tiles ONLY, `largest_n_measured` nowhere. Reconstructed here (not
#: read from the repaired file) because repairing the artifact removed the only instance of the
#: defect, so the regression guard has to carry its own copy of it.
HISTORICAL_ASYMMETRIC = {
    "rungs": {
        "peak_tiles": {"size_lever": "EXHAUSTED", "best_measured_arm": {"spearman": 0.3968},
                       "ceiling_for_a_sequence_model": 0.941},
        "fragments_LB": {"ceiling_for_a_sequence_model": 0.7928, "n": 321123},
        "fragments_M9": {"ceiling_for_a_sequence_model": 0.8014, "n": 318457},
        "sigma70_designed": {"ceiling_for_a_sequence_model": None, "designed_rows": 10898},
    }
}


# ==================================================================================================
# FAIL-SOFT: every degraded path must say WHY, and must never read as evidence of exhaustion
# ==================================================================================================
def test_absent_census_prints_its_reason_and_main_still_exits_0(monkeypatch, capsys, tmp_path):
    """Orientation is exit-0-always, so the only thing that can go wrong is a CLEAN silence. A reader
    who sees nothing under the heading must not conclude the rungs are spent."""
    monkeypatch.setattr(ps, "_CENSUS", tmp_path / "does_not_exist.json")
    rc = ps.main()
    out = capsys.readouterr().out
    assert rc == 0
    assert "UNAVAILABLE" in out
    assert "does_not_exist.json" in out, "the reason must name the file a reader has to go find"
    assert "Do NOT infer that a rung is exhausted" in out


def test_unparseable_census_names_the_exception_and_refuses_the_exhausted_inference(
        monkeypatch, capsys, tmp_path):
    """A broken artifact is indistinguishable from an empty one unless the renderer says so -- the
    documented 'a test that errors looks like a finding' shape."""
    out = _render(monkeypatch, capsys, _write(tmp_path, "{not json"))
    assert "UNAVAILABLE" in out
    assert "JSONDecodeError" in out, "the exception TYPE is the diagnosis; a bare message is not"
    assert "This is NOT evidence that a rung is exhausted." in out


def test_a_census_with_no_rungs_key_refuses_loudly_rather_than_printing_an_empty_table(
        monkeypatch, capsys, tmp_path):
    """Documents which branch a rungs-less census reaches: the blanket `except` (KeyError), NOT a
    clean zero-row table. An empty table under the heading would be the silence this block exists to
    break, so reaching the loud branch is the correct outcome and is pinned here."""
    out = _render(monkeypatch, capsys, _write(tmp_path, {"schema": "mpra-substrate-census-v1"}))
    assert "UNAVAILABLE" in out and "KeyError" in out
    assert "This is NOT evidence that a rung is exhausted." in out


def test_main_still_exits_0_when_the_census_is_unparseable(monkeypatch, capsys, tmp_path):
    """The blanket `except Exception` is deliberate (orientation must not die on an import), so pin
    that it actually holds end to end -- a crash here would block the one script every session is
    told to run before making a scope claim."""
    monkeypatch.setattr(ps, "_CENSUS", _write(tmp_path, "{not json"))
    assert ps.main() == 0


# ==================================================================================================
# THE TRUNCATION GUARD -- the defect this file was written for
# ==================================================================================================
def test_the_substrate_block_carries_no_hard_character_cut(monkeypatch):
    """THE durable guard, source-level and content-independent.

    A content assertion only fires when the string happens to be long enough; a reintroduced
    `[:200]` would pass every other test in this file while silently amputating a longer census.
    The sibling guard `test_the_source_carries_no_hardcoded_regime_figure` works the same way, and
    this file's own comment at the regime block already forbids the pattern in prose.
    """
    import re

    src = SCRIPT.read_text(encoding="utf-8")
    body = src.split("def _print_substrate_block")[1]
    # COMMENTS ARE EXEMPT: the comment explaining the defect quotes the slice that caused it, and a
    # guard that fires on the text documenting it gets suppressed, which is how a guard stops being
    # one.
    code = "\n".join(ln for ln in body.splitlines() if not ln.strip().startswith("#"))
    cuts = re.findall(r"\[:\s*\d{2,}\s*\]", code)
    assert not cuts, f"a fixed character cut is back in the substrate block: {cuts}"
    assert "textwrap" in code, "the warnings must be WRAPPED; that is what replaced the cut"


def test_the_symmetry_warning_keeps_EVERY_asymmetric_field_and_rung_name(
        monkeypatch, capsys, tmp_path):
    """Non-vacuity on the REAL defect shape. `[:150]` kept `best_measured_arm` and dropped
    `size_lever` -- and `size_lever` is the field whose asymmetry caused 321,123 unused fragments to
    be reported as 'no rows available'. The warning about the incident was amputating the
    incident."""
    out = " ".join(_render(monkeypatch, capsys, _write(tmp_path, HISTORICAL_ASYMMETRIC)).split())
    assert "WARNING - census symmetry asymmetric:" in out
    for f in ("size_lever", "best_measured_arm"):
        assert f in out, f"{f} was dropped from the rendered warning"
    for rung in ("fragments_LB", "fragments_M9", "sigma70_designed"):
        assert rung in out, f"{rung} lost or cut mid-token in the rendered warning"


def test_the_ceiling_warning_keeps_the_derivation_it_exists_to_state(
        monkeypatch, capsys, tmp_path):
    """`[:160]` cut every live ceiling warning at `square-ro|oted)`. The step-up derivation IS the
    content -- without it the line reads as an unexplained pair of numbers, and the misquote it
    prevents (0.4583 as a ceiling instead of 0.7928) is exactly a reader taking the raw agreement."""
    census = {"rungs": {
        "fragments_LB": {"single_replicate_agreement": 0.4583,
                         "ceiling_for_a_sequence_model": 0.7928,
                         "n": 321123, "largest_n_measured": 30000,
                         **{f: 1 for f in SYMMETRY_REQUIRED if f != "ceiling_for_a_sequence_model"}},
        "fragments_M9": {"single_replicate_agreement": 0.4731,
                         "ceiling_for_a_sequence_model": 0.8014,
                         "n": 318457, "largest_n_measured": 30000,
                         **{f: 1 for f in SYMMETRY_REQUIRED if f != "ceiling_for_a_sequence_model"}},
    }}
    out = " ".join(_render(monkeypatch, capsys, _write(tmp_path, census)).split())
    assert "0.4583" in out and "0.7928" in out
    assert "(Spearman-Brown step-up, square-rooted)" in out, "the derivation was cut"


# ==================================================================================================
# THE RENDERER'S PURPOSE, and its contract with audit()
# ==================================================================================================
@pytest.mark.skipif(not CENSUS.exists(), reason="census artifact absent (gitignored / fresh clone)")
def test_the_live_census_prints_under_used_NEXT_TO_exhausted(monkeypatch, capsys):
    """The stated purpose of the block: a reader sees every rung side by side and so cannot read one
    rung's `exhausted` as a project-level wall. If only the exhausted rungs rendered, the block would
    reproduce the original defect while looking like a fix."""
    out = _render(monkeypatch, capsys, CENSUS)
    for rung in json.loads(CENSUS.read_text(encoding="utf-8"))["rungs"]:
        assert rung in out, f"rung {rung} is missing from the rendered table"
    assert "under-used:" in out, "no under-used rung rendered -- the contrast is the point"
    assert "exhausted" in out
    assert "Never inherit one." in out


@pytest.mark.skipif(not CENSUS.exists(), reason="census artifact absent (gitignored / fresh clone)")
def test_the_renderer_reads_only_keys_that_audit_actually_emits(monkeypatch, capsys):
    """Couples the two changed files. The renderer indexes `audit()`'s dict by hand, and its blanket
    `except Exception` converts a renamed key into `UNAVAILABLE - KeyError` -- i.e. the block
    degrades into the exact silence it exists to break, exit 0, nothing red. Assert the keys exist,
    AND that the real render does not land in the UNAVAILABLE branch."""
    rep = audit(json.loads(CENSUS.read_text(encoding="utf-8")))
    assert {"symmetry", "utilization", "ceiling_warnings"} <= set(rep)
    assert {"verdict", "detail"} <= set(rep["symmetry"])
    for u in rep["utilization"]:
        assert {"rung", "verdict", "rows_available", "largest_n_measured"} <= set(u)
    assert "UNAVAILABLE" not in _render(monkeypatch, capsys, CENSUS)
