"""Offline pins for the 1000G LD annotation (validates the committed artifact + the valid-population logic).

The capture itself hits Ensembl REST (network); this test asserts the committed annotation's shape + the
load-bearing ancestry finding, and the valid-population threshold logic on a synthetic dict."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.capture_1000g_ld import VALID_R2  # noqa: E402


def test_valid_population_threshold_logic():
    by = {"EUR": {"r2": 0.97}, "AFR": {"r2": 0.33}, "EAS": {"r2": 0.92}, "SAS": {"r2": None}}
    valid = [sp for sp, v in by.items() if isinstance(v.get("r2"), float) and v["r2"] >= VALID_R2]
    assert valid == ["EUR", "EAS"] and VALID_R2 == 0.90


def test_committed_abo_1000g_annotation():
    p = Path(__file__).resolve().parent.parent / "data" / "imputation" / "rs8176719_from_rs657152_1000g_ld.json"
    if not p.exists():
        import pytest
        pytest.skip("1000G LD annotation not captured yet")
    d = json.loads(p.read_text(encoding="utf-8"))
    assert d["tag"] == "rs657152" and d["target"] == "rs8176719"
    bp = d["by_superpopulation"]
    # the load-bearing finding: EUR-valid, AFR-invalid (the imputation map's ancestry limit)
    assert bp["EUR"]["r2"] >= 0.90 and "EUR" in d["valid_populations"]
    assert bp["AFR"]["r2"] < 0.90 and "AFR" not in d["valid_populations"]
