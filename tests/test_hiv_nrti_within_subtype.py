"""Pins the subtype-grouping logic of the NRTI within-subtype transfer check (no dataset / no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.hiv_nrti_within_subtype import _group  # noqa: E402


def test_group_b_vs_nonb_vs_unknown():
    assert _group("B") == "B"
    assert _group("C") == "non-B"
    assert _group("CRF02_AG") == "non-B"
    assert _group("D") == "non-B"
    assert _group("Unknown") == "unknown"
    assert _group("") == "unknown"
    assert _group(None) == "unknown"


if __name__ == "__main__":
    test_group_b_vs_nonb_vs_unknown()
    print("PASS")
