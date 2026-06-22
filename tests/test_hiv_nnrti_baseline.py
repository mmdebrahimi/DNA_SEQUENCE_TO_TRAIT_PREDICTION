"""Pins the load-bearing pure logic of the HIV NNRTI OLS-baseline comparison (no dataset / no network)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.hiv_nnrti_baseline import (  # noqa: E402
    MIN_MUTS, _auc, _build_design_matrix, _confusion, _position_columns,
)


def test_position_columns_filter():
    header = ["SeqID", "EFV", "NVP", "P1", "P103", "P230", "CompMutList", "Px"]
    assert _position_columns(header) == ["P1", "P103", "P230"]  # P<digits> only


def test_build_design_matrix_min_muts_and_presence():
    # 12 rows with K103N; 3 rows with Y181C; consensus elsewhere. MIN_MUTS=10 keeps 103N, drops 181C.
    rows = []
    for i in range(12):
        rows.append({"P103": "N", "P181": "-"})
    for i in range(3):
        rows.append({"P103": "-", "P181": "C"})
    X, feats = _build_design_matrix(rows, ["P103", "P181"])
    assert "103N" in feats           # present in 12 >= MIN_MUTS
    assert "181C" not in feats       # present in 3 < MIN_MUTS -> dropped
    assert MIN_MUTS == 10
    j = feats.index("103N")
    assert X[:12, j].sum() == 12 and X[12:, j].sum() == 0


def test_confusion_metrics():
    pred = np.array([True, True, False, False])
    actual = np.array([True, False, False, True])   # tp=1, fp=1, tn=1, fn=1
    c = _confusion(pred, actual)
    assert (c["tp"], c["fp"], c["tn"], c["fn"]) == (1, 1, 1, 1)
    assert c["sens"] == 0.5 and c["spec"] == 0.5 and c["balanced_accuracy"] == 0.5


def test_auc_perfect_reversed_tie():
    label = np.array([True, True, False, False])
    assert _auc(np.array([3.0, 2.0, 1.0, 0.0]), label) == 1.0   # pos scores all higher
    assert _auc(np.array([0.0, 1.0, 2.0, 3.0]), label) == 0.0   # reversed
    assert _auc(np.array([1.0, 1.0, 1.0, 1.0]), label) == 0.5   # all ties


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
