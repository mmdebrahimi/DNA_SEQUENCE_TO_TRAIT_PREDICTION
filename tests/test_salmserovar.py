"""Pins the Salmonella serovar caller -- the serotype sibling. Pure-logic tests always run; the real
blastn control (synthetic Typhimurium fixture O=4/H1=i/H2=1,2 -> Typhimurium) runs only when blastn resolves.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.salmserovar.runner import (  # noqa: E402
    call_serovar, load_serovar_table, parse_axis_antigen,
)

_FIX = REPO / "tests" / "fixtures" / "salmserovar"


def _blastn():
    try:
        from dna_decode.pathotype.vf_runner import find_blastn
        return find_blastn()
    except Exception:
        return None


# --- pure logic (always run) ---

def test_parse_axis_antigen():
    assert parse_axis_antigen("O__9__01") == ("O", "9")
    assert parse_axis_antigen("H1__d__01") == ("H1", "d")
    assert parse_axis_antigen("H2__1,2__01") == ("H2", "1,2")
    assert parse_axis_antigen("junk") is None


def test_load_serovar_table():
    t = load_serovar_table(_FIX / "db" / "serovar_table.tsv")
    assert t[("4", "i", "1,2")] == "Typhimurium"
    assert t[("9", "d", "-")] == "Typhi"
    assert ("O", "H1", "-") not in t   # header skipped


def test_offline_safe_missing_db():
    r = call_serovar(_FIX / "genome.fna", REPO / "no" / "such" / "dir")
    assert r["status"] == "unavailable" and r["serovar"] is None and r["antigenic_formula"] is None


def test_best_per_axis_is_identity_primary():
    # Regression for the shipped-0.5.2 bug: flagellin (fliC/fljB) alleles cross-hybridize at near-full
    # COVERAGE across antigen types, so coverage-only selection picked the WRONG H antigen (Typhimurium
    # LT2 gave 4:r:1,5,7 instead of 4:i:1,2). The true antigen is the highest-IDENTITY hit.
    from dna_decode.salmserovar.runner import _best_per_axis
    per_allele = {
        "H1__r__8": {"percent_identity": 92.1, "percent_coverage": 101.1, "called": True},  # cross-hybridizer
        "H1__i__5": {"percent_identity": 100.0, "percent_coverage": 100.0, "called": True},  # true antigen
    }
    best = _best_per_axis(per_allele)
    assert best["H1"]["antigen"] == "i"   # identity-primary picks i despite r's higher coverage


# --- real blastn control (skip if blastn absent) ---

@pytest.mark.skipif(_blastn() is None, reason="blastn not installed")
def test_control_typhimurium():
    r = call_serovar(_FIX / "genome.fna", _FIX / "db")
    assert r["status"] == "ok"
    assert r["o_antigen"] == "4" and r["h1_antigen"] == "i" and r["h2_antigen"] == "1,2"
    assert r["antigenic_formula"] == "4:i:1,2"
    assert r["serovar"] == "Typhimurium"


@pytest.mark.skipif(_blastn() is None, reason="blastn not installed")
def test_cli_emits_serovar_call_v0(capsys):
    import json

    from dna_decode.salmserovar.cli import main as sv_main
    rc = sv_main([str(_FIX / "genome.fna"), "--db-dir", str(_FIX / "db"),
                  "--sample-id", "ctl", "--json-only"])
    rec = json.loads(capsys.readouterr().out)
    assert rc == 0 and rec["schema"] == "serovar-call-v0"
    assert rec["serovar"] == "Typhimurium" and rec["organism"] == "Salmonella enterica"
    assert rec["caller"]["caller_is_independent_baseline"] is False   # faithful-to-tool, honest


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
