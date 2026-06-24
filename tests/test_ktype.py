"""Pins the Klebsiella K-antigen (wzi) caller -- the serotype sibling. Pure-logic tests always run; the
real blastn control (committed wzi_1 fixture -> KL1) runs only when blastn resolves.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.ktype.runner import _wzi_num, call_ktype, load_wzi_kl_map  # noqa: E402

_FIX = REPO / "tests" / "fixtures" / "ktype"


def _blastn():
    try:
        from dna_decode.pathotype.vf_runner import find_blastn
        return find_blastn()
    except Exception:
        return None


# --- pure logic (always run) ---

def test_wzi_kl_map_parses_primary_kl():
    m = load_wzi_kl_map(_FIX / "wzi_db" / "wzi.txt")
    assert m["1"] == "KL1"
    assert m["2"] == "KL2"   # 'KL2 (KL30)' -> primary 'KL2' (alternatives dropped)


def test_wzi_num():
    assert _wzi_num("wzi_1") == "1" and _wzi_num("wzi_137") == "137"
    assert _wzi_num("notawzi") is None


def test_offline_safe_missing_db():
    r = call_ktype(_FIX / "genome.fna", REPO / "no" / "such" / "dir")
    assert r["status"] == "unavailable" and r["kl_type"] is None


# --- real blastn control (skip if blastn absent) ---

@pytest.mark.skipif(_blastn() is None, reason="blastn not installed")
def test_control_wzi1_calls_kl1():
    r = call_ktype(_FIX / "genome.fna", _FIX / "wzi_db")
    assert r["status"] == "ok"
    assert r["wzi_allele"] == "wzi_1" and r["kl_type"] == "KL1" and r["predicted_k"] == "KL1"
    assert r["percent_identity"] >= 99.0


@pytest.mark.skipif(_blastn() is None, reason="blastn not installed")
def test_cli_emits_ktype_call_v0(capsys):
    from dna_decode.ktype.cli import main as ktype_main
    rc = ktype_main([str(_FIX / "genome.fna"), "--db-dir", str(_FIX / "wzi_db"),
                     "--sample-id", "ctl", "--json-only"])
    import json
    rec = json.loads(capsys.readouterr().out)
    assert rc == 0 and rec["schema"] == "ktype-call-v0"
    assert rec["predicted_k"] == "KL1" and rec["organism"] == "Klebsiella"
    assert rec["caller"]["caller_is_independent_baseline"] is False   # faithful-to-tool, honest


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
