"""Tests for the DERIVED sourmash-vs-Mash crosscheck bar.

The bar's whole claim to being a commitment rather than a formality is that it is
derived from the committed Mash table's own arithmetic. These tests pin the four
properties that make that true: the near-bar set is computed, the verdicts come from
the imported gate, an incomplete table is refused, and a table with no negative
cannot calibrate the TRIP clause.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from dna_decode.eval import rejection_gates
from scripts.freeze_crosscheck_bar import (
    BAR_PATH,
    IncompleteLineageTable,
    SCHEMA,
    VacuousBar,
    assert_complete,
    build_bar,
    near_bar_pairs,
    non_vacuity,
    pairs_from,
)

REPO = Path(__file__).resolve().parent.parent


def _cell(cohort: str, organism: str, drug: str, counts: dict[str, tuple[int, int]],
          *, partial: bool = False, missing: int = 0) -> dict:
    return {
        "cohort": cohort,
        "organism": organism,
        "drug": drug,
        "partial": partial,
        "n_genomes_missing": missing,
        "lineage_tier_emitted": True,
        "thresholds": {
            thr: {"effective_lineage_N_R": r, "effective_lineage_N_S": s}
            for thr, (r, s) in counts.items()
        },
    }


def _table(cells: list[dict]) -> dict:
    return {"_schema": "provdisjoint-lineage-metrics-v1", "date": "2026-06-12",
            "cells": cells}


# --------------------------------------------------------------------------- #
# the near-bar set is COMPUTED, not hardcoded
# --------------------------------------------------------------------------- #

def test_near_bar_set_is_empty_when_every_count_is_large():
    data = _table([
        _cell("a", "Org", "drugA", {"0.001": (40, 38), "0.005": (33, 31)}),
        _cell("b", "Org2", "drugB", {"0.001": (25, 27), "0.005": (22, 20)}),
    ])
    pairs = pairs_from(data)
    assert len(pairs) == 4
    assert near_bar_pairs(pairs) == [], "no count is within +/-2 of the bar"


def test_near_bar_set_includes_a_count_at_the_bar():
    data = _table([
        _cell("a", "Org", "drugA", {"0.001": (40, 38), "0.005": (3, 31)}),
    ])
    pairs = pairs_from(data)
    near = near_bar_pairs(pairs)
    assert [p.minority for p in near] == [3]
    assert near[0].threshold == "0.005"


def test_near_bar_window_tracks_the_gate_constant(monkeypatch):
    """Raising the bar must move the computed near-bar set, not leave it stale."""
    data = _table([_cell("a", "Org", "d", {"0.001": (14, 20), "0.005": (2, 9)})])
    pairs = pairs_from(data)
    assert sorted(p.minority for p in near_bar_pairs(pairs)) == [2]
    monkeypatch.setattr(rejection_gates, "MIN_EFFECTIVE_LINEAGES", 12)
    assert sorted(p.minority for p in near_bar_pairs(pairs)) == [14]


# --------------------------------------------------------------------------- #
# verdicts come THROUGH the imported gate
# --------------------------------------------------------------------------- #

def test_verdicts_move_when_the_gate_constant_moves(monkeypatch):
    """If this fails, the generator has reimplemented the gate locally."""
    data = _table([_cell("a", "Org", "d", {"0.001": (9, 20), "0.005": (5, 7)})])
    before = [p.mash_verdict for p in pairs_from(data)]
    assert before == [rejection_gates.PASS, rejection_gates.PASS]

    monkeypatch.setattr(rejection_gates, "MIN_EFFECTIVE_LINEAGES", 10)
    after = [p.mash_verdict for p in pairs_from(data)]
    assert after == [rejection_gates.TRIP, rejection_gates.TRIP], (
        "verdicts did not follow MIN_EFFECTIVE_LINEAGES -- the gate is not imported"
    )


def test_minority_class_is_what_the_gate_sees():
    """G8 consumes the SMALLER of the R/S effective-lineage counts."""
    data = _table([_cell("a", "Org", "d", {"0.001": (2, 99)})])
    pair = pairs_from(data)[0]
    assert pair.minority == 2
    assert pair.mash_verdict == rejection_gates.TRIP


# --------------------------------------------------------------------------- #
# an incomplete table is refused
# --------------------------------------------------------------------------- #

def test_partial_cell_is_refused():
    data = _table([_cell("a", "Org", "d", {"0.001": (2, 9)}, partial=True)])
    with pytest.raises(IncompleteLineageTable):
        assert_complete(data)
    with pytest.raises(IncompleteLineageTable):
        build_bar(data, Path("x.json"), b"{}")


def test_missing_genomes_is_refused():
    data = _table([_cell("a", "Org", "d", {"0.001": (2, 9)}, missing=3)])
    with pytest.raises(IncompleteLineageTable):
        build_bar(data, Path("x.json"), b"{}")


def test_absent_lineage_counts_are_refused():
    data = _table([_cell("a", "Org", "d", {})])
    data["cells"][0]["thresholds"] = {"0.001": {"effective_lineage_N_R": 5}}
    with pytest.raises(IncompleteLineageTable):
        pairs_from(data)


# --------------------------------------------------------------------------- #
# non-vacuity
# --------------------------------------------------------------------------- #

def test_all_pass_table_cannot_calibrate_the_trip_clause():
    """A table with no TRIP has no reachable negative, so the bar would be theatre."""
    data = _table([
        _cell("a", "Org", "dA", {"0.001": (40, 38), "0.005": (33, 31)}),
        _cell("b", "Org2", "dB", {"0.001": (25, 27), "0.005": (22, 20)}),
    ])
    with pytest.raises(VacuousBar):
        non_vacuity(pairs_from(data))
    with pytest.raises(VacuousBar):
        build_bar(data, Path("x.json"), b"{}")


def test_non_vacuity_records_the_base_rate():
    data = _table([
        _cell("a", "Org", "dA", {"0.001": (40, 38), "0.005": (33, 31)}),
        _cell("b", "Org2", "dB", {"0.001": (2, 27), "0.005": (22, 20)}),
    ])
    vac = non_vacuity(pairs_from(data))
    assert vac["n_pairs"] == 4
    assert vac["n_pass"] == 3
    assert vac["n_trip"] == 1
    assert vac["pass_base_rate"] == "3/4"
    assert vac["bar"] == rejection_gates.MIN_EFFECTIVE_LINEAGES
    assert "NEARLY VACUOUS" in vac["why_this_matters"]


# --------------------------------------------------------------------------- #
# the committed artifact
# --------------------------------------------------------------------------- #

def _committed_bar() -> dict:
    return json.loads((REPO / BAR_PATH).read_text(encoding="utf-8"))


def test_committed_bar_parses_and_is_frozen_before_any_number():
    bar = _committed_bar()
    assert bar["_schema"] == SCHEMA
    assert bar["frozen_before_any_agreement_number_existed"] is True
    assert bar["derived_from"]["path"] == "wiki/provdisjoint_lineage_metrics.json"
    assert len(bar["derived_from"]["sha256"]) == 64


def test_committed_bar_sha256_matches_the_live_lineage_table():
    import hashlib

    bar = _committed_bar()
    raw = (REPO / bar["derived_from"]["path"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == bar["derived_from"]["sha256"]


def test_committed_bar_has_the_twenty_pair_reference_table():
    bar = _committed_bar()
    assert len(bar["mash_reference_table"]) == 20
    assert bar["non_vacuity"]["pass_base_rate"] == "19/20"
    assert bar["non_vacuity"]["n_trip"] == 1
    assert bar["non_vacuity"]["median_minority_count"] == 14


def test_committed_bar_names_the_three_near_bar_pairs_and_the_one_trip():
    bar = _committed_bar()
    clause2 = next(c for c in bar["clauses"] if c["id"] == 2)
    clause3 = next(c for c in bar["clauses"] if c["id"] == 3)
    assert clause2["near_bar_set_is_computed"] is True
    labels = {p["label"]: p["minority"] for p in clause2["near_bar_pairs"]}
    assert labels == {
        "Klebsiella x ciprofloxacin @ 0.005": 2,
        "Escherichia_coli_Shigella x ciprofloxacin @ 0.005": 4,
        "Escherichia_coli_Shigella x gentamicin @ 0.005": 5,
    }
    assert [p["label"] for p in clause3["trip_pairs"]] == [
        "Klebsiella x ciprofloxacin @ 0.005"
    ]


def test_committed_bar_demotes_ari_and_cluster_count_with_a_reason():
    bar = _committed_bar()
    fields = bar["reported_not_gated"]["fields"]
    assert "adjusted_rand_index" in fields
    assert "n_clusters" in fields
    assert bar["reported_not_gated"]["why"]
    gated = {c["name"] for c in bar["clauses"]}
    assert not any("ari" in n or "rand" in n for n in gated)


def test_committed_bar_labels_clause_1_as_nearly_vacuous():
    """The weak clause must say it is weak, in the artifact, not only in a memo."""
    bar = _committed_bar()
    clause1 = next(c for c in bar["clauses"] if c["id"] == 1)
    assert "VACUOUS" in clause1["strength"].upper()
    assert "19 of 20" in clause1["why_weak"]


def test_committed_bar_imports_the_constant_rather_than_restating_it():
    bar = _committed_bar()
    assert bar["gate_source"]["function"] == "g8_dedup_collapses_balance"
    assert bar["gate_source"]["min_effective_lineages"] == \
        rejection_gates.MIN_EFFECTIVE_LINEAGES


def test_committed_bar_carries_honest_limits():
    bar = _committed_bar()
    limits = " ".join(bar["honest_limits"]).lower()
    assert "not mash's bottom-sketch" in limits
    assert "single trip" in limits
    assert "not cabbage" in limits


# --------------------------------------------------------------------------- #
# the generator must reproduce its own artifact
# --------------------------------------------------------------------------- #

# The ONE key the committed bar carries that `build_bar` does not emit. It was added
# after the freeze, deliberately, to record that clause 3 is implied by clause 1.
# Paths the GENERATOR now emits that the frozen artifact predates. Named, never a
# generic "extra keys are fine", so the exception set cannot quietly grow. It is
# empty of `clause_logic` on purpose: the generator was taught to emit that (so a
# `--force` re-freeze can no longer silently drop it), which is why the reconcile
# gate below checks ADDITIVE-ONLY rather than exact equality.
GENERATOR_ADDITIONS = {
    "/clause_logic/why_build_bar_emits_it",
    "/derived_from/sha256_lf",
    "/derived_from/eol_note",
}


def _additive_diff(frozen: dict, rebuilt: dict, path: str = "") -> list[tuple[str, str]]:
    """Every difference as `(path, ADDED|REMOVED|CHANGED)`, recursing into dicts."""
    out: list[tuple[str, str]] = []
    sentinel = object()
    for key in sorted(set(frozen) | set(rebuilt)):
        a, b = frozen.get(key, sentinel), rebuilt.get(key, sentinel)
        here = path + "/" + key
        if isinstance(a, dict) and isinstance(b, dict):
            out += _additive_diff(a, b, here)
        elif a is sentinel:
            out.append((here, "ADDED"))
        elif b is sentinel:
            out.append((here, "REMOVED"))
        elif a != b:
            out.append((here, "CHANGED"))
    return out


def test_build_bar_reproduces_the_committed_artifact():
    """The reconcile gate.

    The bar's claim to be DERIVED holds only while the generator still reproduces it.
    Without this, a hand-edit of a clause requirement, a moved reference number or a
    changed near-bar set is invisible -- this repo has already recorded a hand-edited
    `wiki/` JSON surviving to a commit because nothing re-parsed it.

    The gate is ADDITIVE-ONLY and the additions are NAMED, rather than a generic
    "extra keys are fine", which is what keeps the exception set from growing. It
    bites in BOTH directions: a re-freeze that CHANGED or DROPPED any frozen value
    fails, and so does one that adds an unnamed key. `clause_logic` used to sit in
    the exception set because the generator never emitted it -- which meant a
    `--force` re-freeze silently deleted it. The generator now emits it, so the
    annotation is pinned here instead of merely excused.
    """
    raw = (REPO / "wiki" / "provdisjoint_lineage_metrics.json").read_bytes()
    regenerated = build_bar(
        json.loads(raw.decode("utf-8")),
        Path("wiki/provdisjoint_lineage_metrics.json"),
        raw,
    )
    committed = _committed_bar()

    delta = _additive_diff(committed, regenerated)
    changed = [pth for pth, kind in delta if kind == "CHANGED"]
    removed = [pth for pth, kind in delta if kind == "REMOVED"]
    added = {pth for pth, kind in delta if kind == "ADDED"}

    assert not changed, f"a re-freeze would CHANGE a frozen value: {changed}"
    assert not removed, f"a re-freeze would DROP a frozen value: {removed}"
    assert added == GENERATOR_ADDITIONS, (
        "a re-freeze adds something unnamed; add it to GENERATOR_ADDITIONS only after "
        f"confirming it weakens nothing: {sorted(added ^ GENERATOR_ADDITIONS)}"
    )
    # clause_logic is EMITTED now, so it can no longer vanish in a re-freeze
    assert "clause_logic" in regenerated
    assert regenerated["clause_logic"]["clause_3_is_implied_by_clause_1"] is True


def test_the_clause_logic_reference_does_not_dangle():
    """The scorer and its tests point a reader at `clause_logic`; it must be there.

    A prose reference to a block that no longer exists is the `--help` over-claim
    pattern: an orientation string refuted by the artifact it describes.
    """
    bar = _committed_bar()
    assert "clause_logic" in bar, "the block the scorer cites is missing"
    assert bar["clause_logic"]["clause_3_is_implied_by_clause_1"] is True
    scorer = (REPO / "scripts" / "sourmash_vs_mash_crosscheck.py").read_text(
        encoding="utf-8")
    assert "clause_logic" in scorer


# --------------------------------------------------------------------------- #
# a frozen bar is not silently re-frozen
# --------------------------------------------------------------------------- #

def test_freeze_refuses_to_overwrite_an_existing_bar(tmp_path, capsys):
    """The most dangerous behaviour in the generator: re-freezing a pre-commitment.

    The whole argument for the bar is that it predates any agreement number. A silent
    rewrite would retire that claim while leaving `frozen_before_any_agreement_number_existed`
    reading true -- so the refusal is load-bearing and was untested.
    """
    from scripts.freeze_crosscheck_bar import main as freeze_main

    existing = tmp_path / "bar.json"
    existing.write_text('{"sentinel": true}', encoding="utf-8")
    assert freeze_main(["--out", str(existing)]) == 1
    assert "REFUSED" in capsys.readouterr().out
    assert json.loads(existing.read_text(encoding="utf-8")) == {"sentinel": True}, (
        "the refusal must leave the existing bar byte-intact"
    )


def test_force_is_what_rewrites_it_and_it_drops_the_post_freeze_annotation(tmp_path):
    """Non-vacuity for the refusal, and the recorded cost of overriding it.

    `--force` genuinely writes, so the refusal above is the only thing standing
    between a frozen bar and a rewrite. It used to also LOSE `clause_logic`, because
    the generator never emitted it -- a rewrite that silently dropped a pre-commitment
    annotation. The generator now emits it, and this test is the non-vacuity proof:
    delete the emission and the assertions below fail.
    """
    from scripts.freeze_crosscheck_bar import main as freeze_main

    out = tmp_path / "bar.json"
    out.write_text('{"sentinel": true}', encoding="utf-8")
    assert freeze_main(["--out", str(out), "--force"]) == 0
    rewritten = json.loads(out.read_text(encoding="utf-8"))
    assert rewritten["_schema"] == SCHEMA
    assert rewritten["non_vacuity"]["pass_base_rate"] == "19/20"
    # A re-freeze PRESERVES the annotation -- this is the fixed behaviour.
    assert "clause_logic" in rewritten, "a --force re-freeze dropped clause_logic"
    assert rewritten["clause_logic"]["clause_3_is_implied_by_clause_1"] is True
    for key in ("discovered", "note", "why_this_is_recorded_rather_than_fixed"):
        assert rewritten["clause_logic"][key], f"clause_logic/{key} lost in a re-freeze"
    # and it records BOTH digests of its own derived-from bytes
    raw = (REPO / "wiki" / "provdisjoint_lineage_metrics.json").read_bytes()
    assert rewritten["derived_from"]["sha256"] == hashlib.sha256(raw).hexdigest()
    assert rewritten["derived_from"]["sha256_lf"] == hashlib.sha256(
        raw.replace(bytes((13, 10)), bytes((10,)))
    ).hexdigest()


def test_to_lf_is_idempotent_and_never_makes_crcrlf():
    from scripts.freeze_crosscheck_bar import _to_lf

    crlf = b"a" + bytes((13, 10)) + b"b" + bytes((13, 10))
    assert _to_lf(crlf) == b"a" + bytes((10,)) + b"b" + bytes((10,))
    assert _to_lf(_to_lf(crlf)) == _to_lf(crlf)
    assert bytes((13,)) not in _to_lf(crlf)
