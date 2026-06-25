"""Pins the S. pneumoniae capsular serotype caller -- the serotype sibling. Pure-logic tests always run;
the real blastn control (synthetic 19F cps fixture -> 19F) runs only when blastn resolves.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.pneumoserotype.runner import call_pneumo_serotype, serotype_of  # noqa: E402

_FIX = REPO / "tests" / "fixtures" / "pneumoserotype"


def _blastn():
    try:
        from dna_decode.pathotype.vf_runner import find_blastn
        return find_blastn()
    except Exception:
        return None


# --- pure logic (always run) ---

def test_serotype_of():
    assert serotype_of("serotype__19F__01") == "19F"
    assert serotype_of("serotype__6B__01") == "6B"
    assert serotype_of("junk") is None


def test_offline_safe_missing_db():
    r = call_pneumo_serotype(_FIX / "genome.fna", REPO / "no" / "such" / "dir")
    assert r["status"] == "unavailable" and r["serotype"] is None


# --- real blastn control (skip if blastn absent) ---

@pytest.mark.skipif(_blastn() is None, reason="blastn not installed")
def test_control_19f():
    r = call_pneumo_serotype(_FIX / "genome.fna", _FIX / "db")
    assert r["status"] == "ok"
    assert r["serotype"] == "19F" and r["best_reference"] == "serotype__19F__01"
    assert r["percent_identity"] >= 99.0


@pytest.mark.skipif(_blastn() is None, reason="blastn not installed")
def test_cli_emits_pneumo_serotype_call_v0(capsys):
    import json

    from dna_decode.pneumoserotype.cli import main as ps_main
    rc = ps_main([str(_FIX / "genome.fna"), "--db-dir", str(_FIX / "db"),
                  "--sample-id", "ctl", "--json-only"])
    rec = json.loads(capsys.readouterr().out)
    assert rc == 0 and rec["schema"] == "pneumo-serotype-call-v0"
    assert rec["serotype"] == "19F" and rec["organism"] == "Streptococcus pneumoniae"
    assert rec["caller"]["caller_is_independent_baseline"] is False   # faithful-to-tool, honest


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
