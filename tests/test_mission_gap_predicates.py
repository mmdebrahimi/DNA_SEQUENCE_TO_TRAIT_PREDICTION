"""Mission-gap predicate guards for `project_state/trait-to-edit/big-idea.md`.

Why this file exists: five of the mission's ten predicates named artifacts that no
producer ever writes, so those gaps could never close on any successful measurement.
`file-exists` is a cheap predicate and that is exactly its hazard -- a misnamed
argument is indistinguishable from "not done yet" by inspection.

The guard is a THREE-STATE CENSUS, not a reachability assertion. Two weaker designs
were tried and rejected on measurement:

  * "every argument resolves on disk or is named by a saved plan" is WRONG -- four
    predicates (both `glm_generator_*`, both `trait_to_edit_*`) are forward-looking
    work with no saved producer plan, so the assertion would fail on correct content.
  * matching against plan PROSE is worse than useless: a plan that DISCUSSES a
    misname (this reconciliation plan does) makes that misname read as reachable.
    A plan `Files:` metadata line does not help either -- the CABBAGE_G8 plan
    declares `scripts/cabbage_lineage_metrics.py`, the script, never the JSON.

The producer-of-record that IS verifiable is CODE: `dna_decode/eval/cabbage_lineage.py`
carries `SIDECAR_NAME = "cabbage_lineage_metrics.json"`. So a predicate is
has-producer-plan iff some module under `dna_decode/` or `scripts/` names the file.
`tests/` is deliberately outside the scan, which is why the misname fixtures below
cannot contaminate it -- no exclusion list is needed.

Non-vacuity is pinned by construction: the five real historical misnames are
regression fixtures and must each classify `no-producer-yet`.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")

REPO = Path(__file__).resolve().parent.parent
MISSION = REPO / "project_state" / "trait-to-edit" / "big-idea.md"
PRODUCER_ROOTS = ("dna_decode", "scripts")

# The five artifact names the mission USED to carry, none of which any producer writes.
# Keeping them here is what makes the census non-vacuous.
HISTORICAL_MISNAMES = (
    "wiki/cabbage_g8_effective_lineages.json",
    "wiki/cabbage_scored_cell_2026.json",
    "wiki/f3_deconfound_calibration_1001g.json",
    "wiki/f3_ibs_crosscheck_1001g.json",
    "wiki/glm_generator_baseline_panel.json",
)

MET = "met"
HAS_PRODUCER = "has-producer-plan"
NO_PRODUCER = "no-producer-yet"


def _mission_text() -> str:
    return MISSION.read_text(encoding="utf-8")


def _parsed():
    import big_idea

    return big_idea, big_idea.parse_big_idea(_mission_text())


def _producer_corpus() -> str:
    """Every committed producer module, concatenated.

    Matching uses the FULL basename including the extension -- deliberate and
    load-bearing. A bare-stem search would let `f3_deconfound_calibration_1001g.json`
    match inside `f3_deconfound_calibration_1001g_2026-10-10.json`, and a misname
    would wear the costume of a real producer. That is the same substring-strictness
    defect already recorded for `scripts/contract_number_audit.py`.
    """
    parts = []
    for root in PRODUCER_ROOTS:
        base = REPO / root
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*.py")):
            try:
                parts.append(path.read_text(encoding="utf-8", errors="replace"))
            except OSError:
                continue
    return "\n".join(parts)


def classify_file_predicate(arg: str, producer_corpus: str) -> str:
    """met | has-producer-plan | no-producer-yet for one `file-exists` argument."""
    if (REPO / arg).exists():
        return MET
    return HAS_PRODUCER if Path(arg).name in producer_corpus else NO_PRODUCER


def _criteria(gap):
    return list(gap.endpoint_criteria)


def test_mission_parses_and_validates():
    big_idea, bi = _parsed()
    big_idea.validate(bi)  # raises BigIdeaError on a malformed predicate list
    assert len(bi.gaps) == 4
    assert sum(g.strategy_budget for g in bi.gaps) == 12, "Sigma budget must stay 12"


def test_predicate_count_is_ten():
    _, bi = _parsed()
    total = sum(len(_criteria(g)) for g in bi.gaps)
    assert total == 10, f"expected 10 predicates, found {total}"


def test_gap2_is_fully_met_after_repointing():
    """Both de-confounder artifacts exist; the gap was unregisterable, not unmet."""
    _, bi = _parsed()
    gap = next(g for g in bi.gaps if g.id == "deconfounder-calibrated")
    args = [c.arg for c in _criteria(gap)]
    assert len(args) == 2
    for arg in args:
        assert (REPO / arg).exists(), f"gap-2 predicate still unresolvable: {arg}"


def test_gap1_first_predicate_is_still_open_and_has_a_producer():
    """The repoint fixes REACHABILITY, not doneness -- the measurement is still owed."""
    _, bi = _parsed()
    gap = next(g for g in bi.gaps if g.id == "amr-label-substrate-validated")
    first = _criteria(gap)[0].arg
    assert first == "wiki/cabbage_lineage_metrics.json"
    assert not (REPO / first).exists(), "gap-1 p1 should still be OPEN"
    assert classify_file_predicate(first, _producer_corpus()) == HAS_PRODUCER


def test_gap1_names_the_source_diverse_cell_and_records_the_risk():
    _, bi = _parsed()
    gap = next(g for g in bi.gaps if g.id == "amr-label-substrate-validated")
    args = [c.arg for c in _criteria(gap)]
    assert "wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json" in args
    text = _mission_text()
    # The choice is only defensible if the flattering-selection risk is stated.
    assert "Flattering-selection risk" in text
    assert "source_diverse" in text
    assert "96-98% USA" in text


def test_file_predicate_census_is_pinned():
    """3 met / 1 has-producer / 5 no-producer across the 9 file-exists predicates."""
    _, bi = _parsed()
    corpus = _producer_corpus()
    states: dict[str, str] = {}
    for gap in bi.gaps:
        for crit in _criteria(gap):
            if crit.kind == "file-exists":
                states[crit.arg] = classify_file_predicate(crit.arg, corpus)
    assert len(states) == 9, states
    counts = {
        state: sum(1 for v in states.values() if v == state)
        for state in (MET, HAS_PRODUCER, NO_PRODUCER)
    }
    assert counts == {MET: 3, HAS_PRODUCER: 1, NO_PRODUCER: 5}, (counts, states)


def test_historical_misnames_all_classify_as_no_producer():
    """Non-vacuity: the census can only catch a misname if a real one fails it."""
    corpus = _producer_corpus()
    for name in HISTORICAL_MISNAMES:
        assert classify_file_predicate(name, corpus) == NO_PRODUCER, (
            f"{name} was a real misname and must not read as reachable"
        )


def test_test_exit_0_predicate_actually_passes():
    """The one non-file predicate, run for real -- it is the 4th met predicate."""
    _, bi = _parsed()
    preds = [c for g in bi.gaps for c in _criteria(g) if c.kind == "test-exit-0"]
    assert len(preds) == 1
    assert preds[0].arg == "pytest tests/test_glm_edit.py -q"
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/test_glm_edit.py", "-q",
         "-p", "no:cacheprovider"],
        cwd=REPO, capture_output=True, text=True, timeout=900,
    )
    assert proc.returncode == 0, proc.stdout[-2000:]


def test_four_of_ten_predicates_are_met():
    """The authoritative count. Was reported as 0/9; it is 4/10."""
    _, bi = _parsed()
    met = 0
    for gap in bi.gaps:
        for crit in _criteria(gap):
            if crit.kind == "file-exists":
                met += (REPO / crit.arg).exists()
            elif crit.kind == "test-exit-0":
                met += 1  # asserted for real in the test above
    assert met == 4, f"expected 4/10 met, got {met}/10"
