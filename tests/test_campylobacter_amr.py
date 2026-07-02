"""Pin the Campylobacter tet/gent curated rules — esp. the gent non-gent-marker exclusion (aad9/spw -> S)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.campylobacter_amr import (  # noqa: E402
    call_cj_gentamicin, call_cj_tetracycline,
)


def test_tet_tetO_family():
    for s in ("tet(O)", "tet(O/M/O)", "tet(O/W/32/O)"):
        assert call_cj_tetracycline([s])["prediction"] == "R", s
    assert call_cj_tetracycline([])["prediction"] == "S"
    assert call_cj_tetracycline(["tet(A)"])["prediction"] == "S"      # efflux, not the Campylobacter mechanism


def test_gent_true_marker_only():
    assert call_cj_gentamicin(["aph(2'')-Ia"])["prediction"] == "R"
    assert call_cj_gentamicin(["aac(3)-IV"])["prediction"] == "R"
    assert call_cj_gentamicin(["aac(6')-Ie-aph(2'')-Ia"])["prediction"] == "R"
    # aad9 (spectinomycin/streptomycin) + spw are NOT gentamicin -> must NOT flip the call
    assert call_cj_gentamicin(["aad9"])["prediction"] == "S"          # non-gent -> S (load-bearing)
    assert call_cj_gentamicin(["spw"])["prediction"] == "S"
    assert call_cj_gentamicin([])["prediction"] == "S"
