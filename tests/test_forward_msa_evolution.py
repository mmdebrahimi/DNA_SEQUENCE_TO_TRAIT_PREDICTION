"""Offline tests for the MSA->evolution-score pipeline (dna_decode/forward/msa_evolution.py).

Pins a2m parsing (match-column extraction + insert handling), the query-pos->column map, the weighted
site-independent log-odds, and the pluggable-upgrade adapter -- all on a tiny synthetic MSA (no D:, no net).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.msa_evolution import (  # noqa: E402
    AA, evolution_table_from_scores, parse_a2m, query_pos_to_col, site_independent_table,
)

# focus MKtAY: M,K = match (pos1,2); t = insert (query pos3, NO column); A,Y = match (pos4,5)
SYNTH = """>focus
MKtAY
>s2
LK.AF
>s3
MR.GY
>s4
MK.AY
"""


def _write(tmp_path, text=SYNTH) -> Path:
    p = tmp_path / "synth.a2m"
    p.write_text(text, encoding="utf-8")
    return p


def test_parse_a2m_extracts_match_columns(tmp_path):
    name, focus_raw, match_cols = parse_a2m(_write(tmp_path))
    assert name == "focus" and focus_raw == "MKtAY"
    assert match_cols[0] == "MKAY"          # focus, insert 't' dropped
    assert match_cols[1] == "LKAF"          # s2, insert '.' dropped
    assert all(len(s) == 4 for s in match_cols)   # 4 match columns, fixed width


def test_query_pos_to_col_skips_inserts():
    # focus MKtAY: pos3 (t) is an insert -> absent; pos1,2 map to cols 0,1; pos4,5 -> cols 2,3
    m = query_pos_to_col("MKtAY")
    assert m == {1: 0, 2: 1, 4: 2, 5: 3}
    assert 3 not in m


def test_query_pos_to_col_focus_gap_consumes_column_without_residue():
    # a '-' in the focus is a match column but not a query residue
    m = query_pos_to_col("M-KA")
    assert m == {1: 0, 2: 2, 3: 3}          # col 1 consumed by '-', no query pos there


def test_site_independent_table_scores_only_match_positions(tmp_path):
    t = site_independent_table(_write(tmp_path), weights=[1.0, 1.0, 1.0, 1.0])
    # scorable positions are 1,2,4,5 (pos3 is an insert) -> no key mentions position 3
    assert not any(m[1:-1] == "3" for m in t)
    # every scorable position emits 19 alternatives
    assert sum(1 for m in t if m[1:-1] == "1") == 19
    # higher score = more preserved: at col0 (M/L/M/M) M is common -> M->M excluded, M->L modest, M->rare very negative
    assert t["M1L"] > t["M1W"]              # L is seen in the column; W never -> W more damaging


def test_site_independent_orientation_matches_blosum_sign(tmp_path):
    # column 2 (pos4): A / A / G / A -> A dominant. A->G (seen) should score higher than A->C (unseen)
    t = site_independent_table(_write(tmp_path), weights=[1.0, 1.0, 1.0, 1.0])
    assert t["A4G"] > t["A4C"]


def test_weights_length_mismatch_raises(tmp_path):
    with pytest.raises(ValueError, match="weights length"):
        site_independent_table(_write(tmp_path), weights=[1.0, 1.0])   # 2 != depth 4


def test_evolution_table_from_scores_is_a_passthrough_adapter():
    scores = {"M1L": 0.5, "M1W": -2.0}
    assert evolution_table_from_scores(scores) == scores
    with pytest.raises(ValueError, match="empty"):
        evolution_table_from_scores({})


def test_empty_msa_raises(tmp_path):
    p = tmp_path / "empty.a2m"
    p.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty MSA"):
        parse_a2m(p)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
