"""The systematic rule-gap screen (scripts/rule_gap_screen.py).

Its job is to name WHICH determinant a cell misses, using the statistic that made the rmt finding
convincing: frequent among missed-R AND absent from the susceptible set. A screen that reports the most
common token in the missed set would be worthless -- passengers are frequent everywhere.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.rule_gap_screen import MIN_MISSED, screen_drug  # noqa: E402


def _cell(rows):
    """rows = [(biosample, label, [determinants])] -> (labels, dets) in the cohort shape."""
    labels = {bs: {"tet": lab} for bs, lab, _ in rows}
    dets = {bs: d for bs, _, d in rows}
    return labels, dets


def _rule(counted):
    """A caller that predicts R iff one of `counted` is present -- the deployed-rule stand-in."""
    return lambda d: {"prediction": "R" if any(x in d for x in counted) else "S"}


def test_a_token_absent_from_S_is_surfaced():
    """The rmt shape: carried by most missed-R, by no susceptible isolate."""
    rows = [(f"R{i}", "R", ["rmtE1", "passenger"]) for i in range(6)]
    rows += [(f"S{i}", "S", ["passenger"]) for i in range(6)]
    labels, dets = _cell(rows)
    out = screen_drug(labels, dets, "tet", _rule(["aac(3)"]))
    assert out["status"] == "screened" and out["n_missed_R"] == 6
    names = [c["determinant"] for c in out["candidates"]]
    assert names[0] == "rmtE1"                 # ranked first by gap
    assert "passenger" not in names            # in 6/6 missed AND 6/6 S -> filtered out


def test_a_token_common_in_S_is_filtered_even_at_100pct_of_missed_R():
    """THE point of the screen. `mtrR` is in 23/23 missed-R and 26/26 S on the real gono tetracycline
    cell; a naive most-common-token screen would report it as the answer."""
    rows = [(f"R{i}", "R", ["mtrR"]) for i in range(8)]
    rows += [(f"S{i}", "S", ["mtrR"]) for i in range(8)]
    labels, dets = _cell(rows)
    out = screen_drug(labels, dets, "tet", _rule(["tet(M)"]))
    assert out["n_missed_R"] == 8              # every R is missed
    assert out["candidates"] == []             # ...and nothing is reported


def test_a_cell_with_no_countable_mechanism_reports_nothing():
    """The gono-azithromycin validation case: 110/110 missed, zero candidates. A screen that emits a
    confident answer where no gene can explain the resistance is worse than useless."""
    rows = [(f"R{i}", "R", ["efflux_marker", "housekeeping"]) for i in range(10)]
    rows += [(f"S{i}", "S", ["efflux_marker", "housekeeping"]) for i in range(10)]
    labels, dets = _cell(rows)
    assert screen_drug(labels, dets, "tet", _rule(["nothing"]))["candidates"] == []


def test_underpowered_and_perfect_cells_are_labelled_not_screened():
    rows = [("R1", "R", ["x"])] + [(f"S{i}", "S", []) for i in range(5)]
    labels, dets = _cell(rows)
    assert screen_drug(labels, dets, "tet", _rule(["nope"]))["status"] == "underpowered"

    rows = [(f"R{i}", "R", ["tet(M)"]) for i in range(5)]
    labels, dets = _cell(rows)
    out = screen_drug(labels, dets, "tet", _rule(["tet(M)"]))
    assert out["status"] == "no_missed_R" and out["n_missed_R"] == 0


def test_a_token_the_rule_already_counts_is_flagged_as_such():
    """`also_in_called_R` keeps the reader from reading a partially-counted token as a pure gap."""
    rows = [(f"R{i}", "R", ["dual"]) for i in range(4)]                  # missed (rule sees nothing here)
    rows += [(f"C{i}", "R", ["dual", "tet(M)"]) for i in range(3)]       # called R via tet(M)
    rows += [(f"S{i}", "S", []) for i in range(5)]
    labels, dets = _cell(rows)
    out = screen_drug(labels, dets, "tet", _rule(["tet(M)"]))
    dual = next(c for c in out["candidates"] if c["determinant"] == "dual")
    assert dual["also_in_called_R"] is True


def test_min_missed_threshold_is_a_real_guard():
    assert MIN_MISSED >= 3      # rates on 1-2 isolates are not interpretable


def test_the_committed_screen_artifact_matches_its_own_verdict():
    """The shipped result must still say what the memo says: no NEW actionable gap."""
    import json
    root = Path(__file__).resolve().parent.parent
    p = root / "wiki" / "rule_gap_screen_2026-08-25.json"
    if not p.exists():
        pytest.skip("screen artifact not present")
    d = json.loads(p.read_text(encoding="utf-8"))
    cells = {(c["organism"], c["drug"]): c for c in d["cells"]}
    # the validation case: every R missed, nothing reported
    azi = cells[("Neisseria gonorrhoeae", "azithromycin")]
    assert azi["n_missed_R"] == azi["n_R"] and azi["candidates"] == []
    # exactly two cells produced candidates, both gonococcal porin
    withc = [c for c in d["cells"] if c["candidates"]]
    assert len(withc) == 2
    assert all(c["candidates"][0]["determinant"].startswith("porB1b") for c in withc)
    assert "hypothesis-generating" in d["honest_scope"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
