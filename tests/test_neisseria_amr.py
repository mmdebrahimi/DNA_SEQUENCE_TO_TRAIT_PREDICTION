"""Pin the N. gonorrhoeae cipro curated rule (gyrA QRDR Ser91/Asp95 -> R; parC accessory-only)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.neisseria_amr import (  # noqa: E402
    call_ng_ciprofloxacin, call_ng_tetracycline,
)


def test_ng_tet_tetM_only():
    # tet(M) (high-level) -> R; rpsJ V57M is accessory-only (low-level, over-calls -> excluded from the call)
    assert call_ng_tetracycline(["tet(M)"])["prediction"] == "R"
    assert call_ng_tetracycline(["tet(M)_1"])["prediction"] == "R"       # startswith tet(M)
    r = call_ng_tetracycline(["rpsJ_V57M"])
    assert r["prediction"] == "S" and r["accessory_rpsJ_V57M"] == ["rpsJ_V57M"]
    assert call_ng_tetracycline([])["prediction"] == "S"
    assert call_ng_tetracycline(["tet(M)", "rpsJ_V57M"])["prediction"] == "R"


def test_gyrA_qrdr_confers_R():
    for sym in ("gyrA_S91F", "gyrA_S91Y", "gyrA_D95N", "gyrA_D95A", "gyrA_D95G"):
        r = call_ng_ciprofloxacin([sym])
        assert r["prediction"] == "R", sym
        assert r["matched_gyrA_qrdr"] == [sym]


def test_no_qrdr_is_S():
    assert call_ng_ciprofloxacin([])["prediction"] == "S"
    assert call_ng_ciprofloxacin(["mtrR", "penA_1", "porB"])["prediction"] == "S"       # irrelevant genes
    assert call_ng_ciprofloxacin(["gyrA_A75S"])["prediction"] == "S"                    # gyrA but NOT a QRDR codon


def test_parC_is_accessory_only():
    r = call_ng_ciprofloxacin(["parC_S87R"])
    assert r["prediction"] == "S"                          # parC alone does NOT flip the binary call
    assert r["accessory_parC_qrdr"] == ["parC_S87R"]
    r2 = call_ng_ciprofloxacin(["gyrA_S91F", "parC_S87R"])
    assert r2["prediction"] == "R" and r2["accessory_parC_qrdr"] == ["parC_S87R"]


def test_rule_is_nonfrozen_scoped():
    r = call_ng_ciprofloxacin(["gyrA_S91F"])
    assert r["rule_status"] == "CURATED_NONFROZEN" and r["rule_scope"] == "scorer_local"
