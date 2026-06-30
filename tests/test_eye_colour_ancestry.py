"""Pin the rs12913832 ancestry-confound quantification (D:-free M2 first cut).

Asserts the sourced 1000G blue-allele frequencies + that the confound summary correctly flags rs12913832
as strongly ancestry-informative (EUR >> all other super-pops) WITHOUT overclaiming (it's structural, not a
proof of inflation)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.eye_colour_ancestry import (  # noqa: E402
    BLUE_ALLELE, RS12913832_BLUE_FREQ_1000G, confound_summary,
)


def test_blue_allele_is_G():
    assert BLUE_ALLELE == "G"


def test_european_concentrated_frequencies():
    f = RS12913832_BLUE_FREQ_1000G
    assert f["EUR"] > 0.6
    assert f["EAS"] < 0.01           # essentially absent in East Asians
    assert f["AFR"] < 0.05
    assert f["EUR"] > 10 * f["AFR"]  # strongly ancestry-informative


def test_confound_summary_flags_ancestry_informative():
    s = confound_summary()
    assert s["is_ancestry_informative"] is True
    assert s["min_blue_freq_pop"] == "EAS"
    assert s["eur_over_min_ratio"] > 100
    # honesty: the summary must NOT claim the accuracy is inflated -- only structural
    assert "NOT a proof" in s["interpretation"]
    assert "within-European" in s["interpretation"]
