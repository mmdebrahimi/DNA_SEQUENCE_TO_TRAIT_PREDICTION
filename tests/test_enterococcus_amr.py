"""Pin the E. faecium curated rules — esp. the gentamicin intrinsic-exclusion (aac(6')-Ii must NOT confer R)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.enterococcus_amr import (  # noqa: E402
    call_efm_ciprofloxacin, call_efm_gentamicin, call_efm_tetracycline,
)


def test_cipro_qrdr():
    for s in ("gyrA_S83Y", "gyrA_E87G", "parC_S80I", "parC_E84K"):
        assert call_efm_ciprofloxacin([s])["prediction"] == "R", s
    assert call_efm_ciprofloxacin([])["prediction"] == "S"
    assert call_efm_ciprofloxacin(["gyrA_A75S"])["prediction"] == "S"     # non-QRDR codon


def test_tet_gene_presence():
    for s in ("tet(M)", "tet(L)", "tet(S)", "tet(O)"):
        assert call_efm_tetracycline([s])["prediction"] == "R", s
    assert call_efm_tetracycline(["vanA"])["prediction"] == "S"


def test_gent_intrinsic_exclusion():
    # HIGH-level gent needs aph(2''); the INTRINSIC aac(6')-Ii must NOT flip the call (the over-call trap)
    assert call_efm_gentamicin(["aac(6')-Ie-aph(2'')-Ia"])["prediction"] == "R"
    assert call_efm_gentamicin(["aph(2'')-Ia"])["prediction"] == "R"
    assert call_efm_gentamicin(["aac(6')-Ii"])["prediction"] == "S"       # intrinsic -> S (load-bearing)
    assert call_efm_gentamicin([])["prediction"] == "S"
