"""Pins the S. pneumoniae beta-lactam PBP->MIC engine (non-frozen). Pure-logic on a synthetic table
(always runs) + a real-DB smoke that skips when the gitignored CDC table is absent."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.organism_rules.pneumo_betalactam import (  # noqa: E402
    RULE_STATUS, apt_key, predict_mic, predict_rs,
)

_REAL_DB = REPO / "data" / "pneumo_betalactam_db" / "Ref_PBPtype_MIC.csv"

# synthetic mini lookup: APT -> {drug_col: mic}
_TABLE = {"1-2-2": {"PEN": 0.03, "MER": 0.06}, "13-11-16": {"PEN": 4.0, "MER": 1.0}}


def test_branding_non_frozen():
    assert RULE_STATUS == "KNOWLEDGE_BASELINE"


def test_apt_key():
    assert apt_key("1", "2", "2") == "1-2-2"
    assert apt_key("209", "253", "NEW") == "209-253-NEW"


def test_predict_mic_lookup():
    assert predict_mic(_TABLE, "1", "2", "2", "penicillin") == 0.03
    assert predict_mic(_TABLE, "13", "11", "16", "meropenem") == 1.0


def test_novel_pbp_type_is_nocall():
    # 'NEW' allele -> no fabricated MIC
    assert predict_mic(_TABLE, "209", "253", "NEW", "penicillin") is None


def test_absent_combo_is_nocall():
    assert predict_mic(_TABLE, "99", "99", "99", "penicillin") is None


def test_predict_rs_breakpoint_context():
    # PEN MIC 4.0 -> R at meningitis (S<=0.06) but R at non-meningitis too (R>=8? no -> I). Check meningitis R.
    r = predict_rs(_TABLE, "13", "11", "16", "penicillin", "meningitis")
    assert r["prediction"] == "R" and r["mic"] == 4.0 and r["apt"] == "13-11-16"
    # susceptible isolate
    s = predict_rs(_TABLE, "1", "2", "2", "penicillin", "meningitis")
    assert s["prediction"] == "S" and s["mic"] == 0.03


@pytest.mark.skipif(not _REAL_DB.exists(), reason="CDC PBP-MIC table not present (gitignored)")
def test_real_db_loads():
    from dna_decode.organism_rules.pneumo_betalactam import load_pbp_mic_table
    t = load_pbp_mic_table(_REAL_DB)
    assert len(t) > 100 and "0-0-0" in t and "PEN" in t["0-0-0"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
