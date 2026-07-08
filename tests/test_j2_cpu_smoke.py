"""Offline test for the J2 CPU-smoke assay selector (no torch/network at import)."""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _load():
    spec = importlib.util.spec_from_file_location("j2cpu", ROOT / "scripts" / "j2_cpu_smoke.py")
    m = importlib.util.module_from_spec(spec); sys.modules["j2cpu"] = m; spec.loader.exec_module(m)
    return m


def test_smallest_assays_filters_and_orders():
    M = _load()
    ref = [
        {"DMS_id": "big", "target_seq": "M" * 300, "DMS_number_single_mutants": "999"},
        {"DMS_id": "tiny_ok", "target_seq": "M" * 40, "DMS_number_single_mutants": "50"},
        {"DMS_id": "tiny_toofew", "target_seq": "M" * 30, "DMS_number_single_mutants": "5"},   # <20 -> dropped
        {"DMS_id": "mid", "target_seq": "M" * 100, "DMS_number_single_mutants": "40"},
        {"DMS_id": "noseq", "target_seq": "", "DMS_number_single_mutants": "40"},               # no seq -> dropped
    ]
    got = M.smallest_assays(ref, k=5)
    assert [r["DMS_id"] for r in got] == ["tiny_ok", "mid", "big"]   # sorted by len, filtered, k-capped
    assert M.smallest_assays(ref, k=1)[0]["DMS_id"] == "tiny_ok"
