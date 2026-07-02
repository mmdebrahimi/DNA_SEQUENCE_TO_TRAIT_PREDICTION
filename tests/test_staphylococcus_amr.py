"""Pin the S. aureus cipro curated rule (gyrA 84/85/88 OR grlA/parC 80/84 QRDR -> R; BOTH primary)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.staphylococcus_amr import call_sa_ciprofloxacin  # noqa: E402


def test_gyrA_qrdr_confers_R():
    for sym in ("gyrA_S84L", "gyrA_S85P", "gyrA_E88K"):
        assert call_sa_ciprofloxacin([sym])["prediction"] == "R", sym


def test_parC_grlA_is_PRIMARY_in_staph():
    # unlike N. gonorrhoeae, parC/grlA alone DOES confer resistance in S. aureus (first-step target)
    r = call_sa_ciprofloxacin(["parC_S80F"])
    assert r["prediction"] == "R" and r["matched_parC_qrdr"] == ["parC_S80F"]
    assert call_sa_ciprofloxacin(["parC_E84K"])["prediction"] == "R"


def test_no_qrdr_is_S():
    assert call_sa_ciprofloxacin([])["prediction"] == "S"
    assert call_sa_ciprofloxacin(["mecA", "blaZ"])["prediction"] == "S"        # irrelevant genes
    assert call_sa_ciprofloxacin(["gyrA_A75S"])["prediction"] == "S"           # gyrA but NOT a QRDR codon


def test_both_targets_recorded():
    r = call_sa_ciprofloxacin(["gyrA_S84L", "parC_S80F"])
    assert r["prediction"] == "R"
    assert r["matched_gyrA_qrdr"] == ["gyrA_S84L"] and r["matched_parC_qrdr"] == ["parC_S80F"]
    assert r["rule_status"] == "CURATED_NONFROZEN"
