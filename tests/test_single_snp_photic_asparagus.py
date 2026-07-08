"""Photic-sneeze + asparagus-anosmia single-locus falsification cells (sourced directions + circularity guard)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.data import single_snp_traits as T  # noqa: E402


def test_photic_call_direction_strand_agnostic():
    # C risk allele -> sneezer (Eriksson 2010); homozygous T -> non-sneezer. C==G, T==A across strands.
    assert T.call_photic("CC") == "sneezer"
    assert T.call_photic("CT") == "sneezer"
    assert T.call_photic("TT") == "non-sneezer"
    assert T.call_photic("GG") == "sneezer"        # minus strand
    assert T.call_photic("AA") == "non-sneezer"
    assert T.call_photic("") == T.INDETERMINATE


def test_asparagus_call_direction_strand_agnostic():
    # A allele -> can-smell (dominant); homozygous G -> anosmic. A==T, G==C across strands.
    assert T.call_asparagus("AA") == "can-smell"
    assert T.call_asparagus("AG") == "can-smell"
    assert T.call_asparagus("GG") == "anosmic"
    assert T.call_asparagus("TT") == "can-smell"   # minus strand
    assert T.call_asparagus("CC") == "anosmic"
    assert T.call_asparagus("") == T.INDETERMINATE


def test_circularity_guard_drops_genotype_referencing_reports():
    # load-bearing: a self-report that names the SNP/genotype is contaminated -> must NOT be scored.
    assert T.bin_photic("photic sneezer with the snp") is None
    assert T.bin_photic("i have the snp but no sneezing") is None
    assert T.bin_photic("rs10427255") is None
    assert T.bin_photic("cc") is None
    assert T.bin_asparagus("gg - but... i can smell it") is None
    assert T.bin_asparagus("ag") is None


def test_photic_binner_clean_reports():
    assert T.bin_photic("no sneezing") == "non-sneezer"
    assert T.bin_photic("no") == "non-sneezer"
    assert T.bin_photic("photic sneezer") == "sneezer"
    assert T.bin_photic("yes") == "sneezer"
    assert T.bin_photic("not in light but do have sneezing fits on occasion") is None  # not photic-specific
    assert T.bin_photic("don't know") is None


def test_asparagus_binner_clean_reports():
    assert T.bin_asparagus("yes") == "can-smell"
    assert T.bin_asparagus("can't smell") == "anosmic"
    assert T.bin_asparagus("don't know") is None


def test_traits_registered_with_sources():
    for k in ("photic", "asparagus"):
        assert k in T.TRAITS
        assert "Eriksson 2010" in T.TRAITS[k].source
        assert T.TRAITS[k].tier == "WEAK_ASSOCIATION_CONTRAST"


def test_verdict_vs_null_thresholds():
    import sys as _s
    _s.path.insert(0, str(REPO / "scripts"))
    from single_snp_opensnp_validate import _verdict_vs_null
    assert _verdict_vs_null(0.90, 0.60) == "ABOVE_NULL"
    assert _verdict_vs_null(0.602, 0.574) == "NEAR_CHANCE"
    assert _verdict_vs_null(0.574, 0.852) == "BELOW_NULL"
