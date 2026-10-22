"""The wall memo must stay a DIAGNOSIS and must not be readable as a RESULT.

Two guards earn their keep here.

`test_the_sidecar_is_parseable_json` exists because this repo committed an invalid hand-edited JSON
artifact for a day and cited it as evidence -- a scope note written with Python implicit string
concatenation. Nothing read it, so nothing complained.

`test_the_memo_quotes_none_of_the_frozen_bars_headline_figures` is the stronger one. The memo describes
a pre-registered run that HAS NOT BEEN EXECUTED; a reader who found one of the bar's numbers in it
could reasonably conclude an arm had run. The numbers the memo IS allowed to quote are the ones this
work measured -- the tile count and the noise ceiling -- and those are asserted present so the test
cannot pass by the memo simply saying nothing.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MEMO = REPO / "wiki" / "tile_chain_substrate_wall_2026-10-21.md"
SIDECAR = REPO / "wiki" / "tile_chain_substrate_wall_2026-10-21.json"
PROBE = REPO / "wiki" / "promcalc_toolchain_probe_2026-10-21.json"
BAR = REPO / "wiki" / "glm_tile_hybrid_primary_acceptance_bar.json"

WALL_CLASSES = {"code-closable", "data-absent", "external-wall"}


def test_both_artifacts_exist():
    assert MEMO.is_file() and SIDECAR.is_file()


def test_the_sidecar_is_parseable_json():
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    assert rec["schema"] == "tile-chain-substrate-wall-v1"


def test_exactly_two_dependencies_with_DISTINCT_classes_from_the_closed_vocabulary():
    """One undifferentiated 'substrate missing' is the reading this memo exists to replace."""
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    deps = rec["dependencies"]
    assert len(deps) == 2
    classes = [d["wall_class"] for d in deps]
    assert set(classes) <= WALL_CLASSES
    assert len(set(classes)) == 2, "two dependencies sharing one class would defeat the point"


def test_every_dependency_says_what_is_needed_and_from_whom():
    """The field that separates 'write more code' from 'a drive, a host or an install must appear'."""
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    for d in rec["dependencies"]:
        assert d["needed_from_whom"].strip()


def test_the_memo_names_both_dependencies_and_both_classes():
    text = MEMO.read_text(encoding="utf-8")
    assert "code-closable" in text and "data-absent" in text
    assert "promcalc" in text and "peak-tile" in text


def test_the_memo_records_the_correction():
    """The substrate was asserted present on the strength of the wrong directory."""
    text = MEMO.read_text(encoding="utf-8").lower()
    assert "correction" in text
    assert "wrong" in text
    assert json.loads(SIDECAR.read_text(encoding="utf-8"))["correction"]["why_it_was_wrong"].strip()


def test_the_memo_quotes_none_of_the_frozen_bars_headline_figures():
    """A memo carrying the bar's numbers could be read as an arm having run. None has."""
    text = MEMO.read_text(encoding="utf-8")
    forbidden = ["0.4254", "0.4847", "0.4693", "0.4807", "0.4217", "0.3968", "0.3112", "0.5366"]
    present = [f for f in forbidden if f in text]
    assert not present, f"memo quotes the bar's figures: {present}"


def test_the_memo_DOES_carry_the_numbers_this_work_measured():
    """NON-VACUITY for the test above: a memo that quoted nothing would also pass it."""
    text = MEMO.read_text(encoding="utf-8")
    assert "44106" in text, "the reproduced tile count is the code-closable half's evidence"
    assert "0.9410" in text, "the reproduced noise ceiling is the other half of that evidence"


def test_the_memo_states_that_no_arm_was_run():
    text = MEMO.read_text(encoding="utf-8")
    assert "not a result" in text.lower()
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    assert "No arm" in rec["scope"]
    assert any("not a result" in s.lower() for s in rec["what_this_does_not_establish"])


def test_the_satisfiable_verdict_is_not_presented_as_a_reproduction():
    """The bar's anchor gate is a separate and stricter question; REFUSED_ANCHOR_MISMATCH stays live."""
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    joined = " ".join(rec["what_this_does_not_establish"])
    assert "REFUSED_ANCHOR_MISMATCH" in joined


def test_the_isolation_question_is_RESOLVED_and_its_resolution_travels():
    """Pin that the question is ANSWERED with evidence, not that the answer is a given value.

    This asserted `isolation_required is True` until 2026-10-22, on the strength of stage B's
    docstring ("promoter_calculator pulls sklearn 1.9.1"). That was measured FALSE: the distribution
    is a py3-none-any pure-Python wheel with unpinned deps, and a dry-run resolution would install
    exactly one package while leaving sklearn alone. The old test was pinning a relayed premise.

    The durable guard is therefore: the flag is a BOOL (so the question is settled rather than open),
    and the amendment carries the refuted claim plus its evidence. A sidecar that quietly dropped
    either would let the next reader re-adopt the premise that just cost an authority fork.
    """
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    dep = next(d for d in rec["dependencies"] if d["wall_class"] == "data-absent")
    assert isinstance(dep["isolation_required"], bool)

    am = rec["amendment_2"]
    assert "sklearn 1.9.1" in am["claim_refuted"]
    assert "dry-run" in am["evidence"]
    assert "promotercalculator" in am["trap_recorded"]
    # the wall must still be recorded as real, just reclassified -- not quietly dissolved
    assert "REFUSED_ANCHOR_MISMATCH" in am["wall_reclassified"]


def test_the_stale_and_sound_halves_of_the_bars_note_are_BOTH_recorded():
    """A reader must discard neither clause wholesale nor trust both."""
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    audit = rec["frozen_bar_execution_notes_audit"]
    assert "STALE" in audit["no_uv"]
    assert "RIGHT" in audit["no_d_data_tree"]
    assert audit["bar_edited"] is False


@pytest.mark.skipif(not BAR.is_file(), reason="frozen bar artifact absent")
def test_the_frozen_bar_was_not_modified():
    """The memo is allowed to audit a frozen pre-registration; it is not allowed to edit one."""
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert bar["frozen"] == "2026-10-17"
    assert bar["frozen_before_any_number_was_seen"] is True
    assert bar["execution_notes"]["status"].startswith("PRE-REGISTRATION ONLY")
    assert "no uv" in bar["execution_notes"]["not_runnable_on"], \
        "the bar still carries its stale clause; this memo audits it rather than editing it"


@pytest.mark.skipif(not PROBE.is_file(), reason="probe artifact absent")
def test_the_memo_cites_a_probe_artifact_that_exists_and_agrees():
    rec = json.loads(SIDECAR.read_text(encoding="utf-8"))
    dep = next(d for d in rec["dependencies"] if d["wall_class"] == "data-absent")
    probe = json.loads(PROBE.read_text(encoding="utf-8"))
    assert dep["probe_artifact"] == "wiki/promcalc_toolchain_probe_2026-10-21.json"
    assert dep["probe_verdict"] == probe["verdict"]
