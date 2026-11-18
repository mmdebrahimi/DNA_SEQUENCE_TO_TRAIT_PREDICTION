"""Tests for the decay-curve confound decomposition (`scripts/decay_curve_confound.py`).

All offline and substrate-free: every test drives the pure `decompose`.

The properties worth pinning are the ones that would let a confounded design read as a clean one:

  * the floor must actually gate separability, or "separable" becomes decoration;
  * a single group must NOT report as collinear, because there is nothing to be collinear with;
  * misaligned inputs must raise rather than zip short, since a short zip silently attributes one
    variant's distance to another variant's group;
  * the verdict must stay inside its closed enum.
"""
from __future__ import annotations

import pytest

from scripts.decay_curve_confound import VERDICTS, decompose

BINS = ((0, 1), (2, 3), (4, 7))
FLOOR = 30


def _mk(spec: list[tuple[str, int, int]]):
    """spec rows are (group, distance, count) -> parallel (distances, labels)."""
    d: list[int] = []
    g: list[str] = []
    for grp, dist, n in spec:
        d.extend([dist] * n)
        g.extend([grp] * n)
    return d, g


# ------------------------------------------------------------------ the core verdict


def test_perfect_collinearity_is_reported_as_confounded():
    d, g = _mk([("A", 1, 40), ("B", 2, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert r["verdict"] == "confounded:collinear"
    assert r["max_dominant_pct_across_scored_bins"] == 100.0


def test_a_group_spanning_two_distances_at_the_floor_is_separable():
    d, g = _mk([("A", 1, 40), ("A", 2, 40), ("B", 1, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert r["verdict"] == "separable:within-group-possible"
    assert r["within_group_separation"]["A"]["separable"] is True
    assert r["within_group_separation"]["B"]["separable"] is False


def test_the_floor_actually_gates_separability():
    """Same data, higher floor -> no longer separable. Otherwise the floor is decoration."""
    d, g = _mk([("A", 1, 40), ("A", 2, 40), ("B", 1, 40)])
    assert decompose(d, g, bins=BINS, min_bin_n=30)["verdict"] == "separable:within-group-possible"
    assert decompose(d, g, bins=BINS, min_bin_n=50)["verdict"] == "confounded:collinear"


def test_a_single_group_is_indeterminate_not_confounded():
    """With one group there is nothing for distance to be collinear WITH."""
    d, g = _mk([("A", 1, 40), ("A", 2, 40)])
    assert decompose(d, g, bins=BINS, min_bin_n=FLOOR)["verdict"] == "indeterminate:single-group"


def test_nothing_clearing_the_floor_is_indeterminate():
    d, g = _mk([("A", 1, 5), ("B", 2, 5)])
    assert decompose(d, g, bins=BINS, min_bin_n=FLOOR)["verdict"] == "indeterminate:no-scored-bins"


def test_verdict_always_inside_the_closed_enum():
    for spec in ([("A", 1, 40), ("B", 2, 40)],
                 [("A", 1, 40), ("A", 2, 40), ("B", 1, 40)],
                 [("A", 1, 40)],
                 [("A", 1, 1)]):
        d, g = _mk(spec)
        assert decompose(d, g, bins=BINS, min_bin_n=FLOOR)["verdict"] in VERDICTS


# -------------------------------------------------------------------- input safety


def test_misaligned_inputs_raise_rather_than_zip_short():
    with pytest.raises(ValueError):
        decompose([1, 2, 3], ["A", "B"], bins=BINS, min_bin_n=1)


def test_empty_input_is_indeterminate_not_a_crash():
    r = decompose([], [], bins=BINS, min_bin_n=FLOOR)
    assert r["verdict"] == "indeterminate:no-scored-bins"
    assert r["n_total"] == 0


# ------------------------------------------------------------------- the reported detail


def test_every_preregistered_bin_appears_even_when_empty():
    d, g = _mk([("A", 1, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert [row["bin_index"] for row in r["bins"]] == [0, 1, 2]
    assert r["bins"][1]["n"] == 0 and r["bins"][1]["scored"] is False


def test_composition_percentages_sum_to_about_a_hundred_per_bin():
    d, g = _mk([("A", 1, 30), ("B", 1, 10)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    comp = r["bins"][0]["composition"]
    assert round(sum(v["pct_of_bin"] for v in comp.values()), 1) == 100.0
    assert comp["A"]["n"] == 30 and comp["B"]["n"] == 10


def test_dominant_group_is_the_actual_majority():
    d, g = _mk([("A", 1, 10), ("B", 1, 30)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert r["bins"][0]["dominant_group"] == "B"
    assert r["bins"][0]["dominant_pct"] == 75.0


def test_the_supplied_floor_and_bins_are_echoed_into_the_result():
    """A later run must not be able to quietly substitute a different pre-registration."""
    d, g = _mk([("A", 1, 40), ("B", 2, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=37)
    assert r["min_bin_n_supplied_by_caller"] == 37
    assert r["bins_supplied_by_caller"] == [[0, 1], [2, 3], [4, 7]]


def test_within_group_cells_count_every_distance_not_only_scorable_ones():
    """The tiny off-mode cells are the evidence the confound is unbreakable; they must be visible."""
    d, g = _mk([("A", 1, 40), ("A", 2, 3), ("B", 2, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert r["within_group_separation"]["A"]["cells"] == {"d=1": 40, "d=2": 3}
    assert r["within_group_separation"]["A"]["distances_at_or_above_floor"] == [1]


def test_reading_text_names_the_worst_dominance_when_confounded():
    d, g = _mk([("A", 1, 40), ("B", 2, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert "100.0%" in r["reading"]


def test_confounded_result_states_what_it_does_not_establish():
    """A confounded design is uninformative, NOT evidence for the null."""
    d, g = _mk([("A", 1, 40), ("B", 2, 40)])
    r = decompose(d, g, bins=BINS, min_bin_n=FLOOR)
    assert "not evidence" in r["what_this_does_not_establish"]


# ------------------------------------------------- the real case this module was built for


def test_reproduces_the_shape_of_the_measured_gse108535_confound():
    """The committed artifact records bin 0-1 at 100.0% single and bin 2-3 at 98.6% double.

    Rebuilt here from the recorded per-cell counts so the shape is pinned without the substrate.
    """
    d, g = _mk([
        ("12T->G", 1, 1331), ("12T->G", 2, 16), ("12T->G", 3, 3),
        ("12T->G/7T->C", 2, 1371), ("12T->G/7T->C", 3, 3), ("12T->G/7T->C", 4, 1),
    ])
    r = decompose(d, g, bins=((0, 1), (2, 3), (4, 7), (8, 15)), min_bin_n=30)

    assert r["verdict"] == "confounded:collinear"
    b01, b23 = r["bins"][0], r["bins"][1]
    assert b01["n"] == 1331 and b01["dominant_pct"] == 100.0
    assert b01["dominant_group"] == "12T->G"
    assert b23["n"] == 1393 and b23["dominant_pct"] == 98.6
    assert b23["dominant_group"] == "12T->G/7T->C"
    # Neither group spans two distances at the floor -- the confound is unbreakable here.
    assert not any(w["separable"] for w in r["within_group_separation"].values())
