"""Pin the dataset-hunt scorecard logic (F1)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.dataset_candidate_scorecard import Candidate, GATE_KEYS, decoder_type, rank, score  # noqa: E402


def _gates(**over):
    g = {k: "pass" for k in GATE_KEYS}
    g.update(over)
    return g


def test_all_pass_learned_niche():
    c = Candidate("Deep", "yeast", "growth", _gates(), depth_estimate=1000)
    s = score(c)
    assert s["verdict"] == "PASS" and s["decoder_type"] == "learned-niche"


def test_catalog_relaxes_depth():
    # no depth, but a curated catalog -> deterministic path, G6 satisfied -> PASS
    c = Candidate("Cat", "human", "variant->disease", _gates(G6_depth_or_catalog="fail"),
                  has_curated_catalog=True)
    assert score(c)["verdict"] == "PASS"
    assert decoder_type(c) == "deterministic"


def test_any_fail_rejects():
    c = Candidate("Circular", "x", "y", _gates(G2_non_circular="fail"), depth_estimate=500)
    assert score(c)["verdict"] == "REJECT"


def test_unknown_is_verify_not_pass():
    c = Candidate("Gap", "fly", "trait", _gates(G3_sampling_independent="unknown"), depth_estimate=200)
    assert score(c)["verdict"] == "VERIFY"


def test_all_pass_but_no_paradigm_rejects():
    # shallow + no catalog -> neither paradigm -> REJECT even with all gates pass
    c = Candidate("Shallow", "x", "y", _gates(), depth_estimate=20)
    assert decoder_type(c) == "neither"
    assert score(c)["verdict"] == "REJECT"


def test_bad_gate_value_raises():
    with pytest.raises(ValueError):
        Candidate("Bad", "x", "y", _gates(G1_accessible="maybe"))


def test_rank_orders_pass_first():
    a = Candidate("A", "x", "y", _gates(G2_non_circular="unknown"), depth_estimate=200)   # VERIFY
    b = Candidate("B", "x", "y", _gates(), depth_estimate=1000)                            # PASS
    out = rank([a, b])
    assert [r["name"] for r in out] == ["B", "A"]
