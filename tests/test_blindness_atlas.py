"""Tests for the determinant-blindness atlas pure core (scripts/build_blindness_atlas.blindness_row)."""
from scripts.build_blindness_atlas import blindness_row


def test_no_measured_r_is_unscorable():
    b = blindness_row([("S", "S", False), ("S", "S", True)])
    assert b["n_R"] == 0
    assert b["invisible_fraction"] is None


def test_fully_visible_cell():
    # every measured-R is called R -> zero invisible
    b = blindness_row([("R", "R", True), ("R", "R", True), ("S", "S", False)])
    assert b["n_R"] == 2
    assert b["n_invisible"] == 0
    assert b["invisible_fraction"] == 0.0


def test_fully_invisible_cell_truly_invisible():
    # every measured-R called S AND carries no determinant token (the azithro-mtr shape)
    b = blindness_row([("R", "S", False), ("R", "S", False)])
    assert b["n_R"] == 2
    assert b["n_invisible"] == 2
    assert b["invisible_fraction"] == 1.0
    assert b["n_truly_invisible"] == 2
    assert b["n_rule_limited"] == 0


def test_split_truly_vs_rule_limited():
    # 4 R: 1 visible, 1 invisible-with-no-determinant, 1 invisible-with-a-determinant, 1 visible
    b = blindness_row([
        ("R", "R", True),    # visible
        ("R", "S", False),   # truly invisible (zero determinant)
        ("R", "S", True),    # rule-limited (has a determinant, rule didn't fire)
        ("R", "R", True),    # visible
    ])
    assert b["n_R"] == 4
    assert b["n_invisible"] == 2
    assert b["invisible_fraction"] == 0.5
    assert b["n_truly_invisible"] == 1
    assert b["n_rule_limited"] == 1


def test_prediction_case_insensitive():
    # a lowercase 'r' prediction still counts as visible
    b = blindness_row([("R", "r", True)])
    assert b["n_invisible"] == 0
