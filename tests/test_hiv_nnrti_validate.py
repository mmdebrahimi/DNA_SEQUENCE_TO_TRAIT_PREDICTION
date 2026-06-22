"""Pins the load-bearing pure logic of the HIV NNRTI validation harness (no dataset / no network).

The headline AUC is only trustworthy if the rank-based Mann-Whitney AUC + the mutation parse are correct.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.hiv_nnrti_validate import (  # noqa: E402
    _auc_rank, _observed_rt_mutations, _parse_fold,
)


def test_auc_perfect_and_reversed_and_tie():
    assert _auc_rank([3, 2, 1], [0, -1, -2]) == 1.0        # all pos > all neg
    assert _auc_rank([0], [1]) == 0.0                       # pos < neg
    assert _auc_rank([1, 2], [1, 2]) == 0.5                 # identical distributions -> all ties = 0.5
    assert _auc_rank([2, 2], [1, 1]) == 1.0
    assert _auc_rank([], [1, 2]) is None                    # empty group -> None


def test_auc_partial_separation_half_credit_on_ties():
    # pos={2,1}, neg={1,0}: pairs (2,1)>,(2,0)>,(1,1)tie,(1,0)> -> (3 + 0.5)/4 = 0.875
    assert _auc_rank([2, 1], [1, 0]) == 0.875


def test_parse_fold():
    assert _parse_fold("25.0") == 25.0
    assert _parse_fold("NA") is None
    assert _parse_fold("") is None
    assert _parse_fold(".") is None
    assert _parse_fold(">100") == 100.0    # censored kept at the numeric bound (v0)
    assert _parse_fold("<0.5") == 0.5


def test_observed_rt_mutations_forms_wt_pos_mut():
    # P103=N -> K103N (wt K); '-' -> nothing; P181 mixture 'CY' -> Y181C only (Y == wt, skipped)
    row = {f"P{p}": "-" for p in (100, 101, 103, 106, 181, 188, 190, 230)}
    row["P103"] = "N"
    row["P181"] = "CY"
    row["P188"] = "L"
    muts = _observed_rt_mutations(row)
    assert muts == {"K103N", "Y181C", "Y188L"}


def test_observed_rt_mutations_consensus_only_is_empty():
    row = {f"P{p}": "-" for p in (100, 101, 103, 106, 181, 188, 190, 230)}
    assert _observed_rt_mutations(row) == set()


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
