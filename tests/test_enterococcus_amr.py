"""Pin the E. faecium curated rules — esp. the gentamicin intrinsic-exclusion (aac(6')-Ii must NOT confer R)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.enterococcus_amr import (  # noqa: E402
    call_efm_ciprofloxacin, call_efm_gentamicin, call_efm_tetracycline,
    call_efm_vancomycin, call_efm_teicoplanin,
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


def test_vancomycin_ligase_presence():
    # any van ligase -> R; accessory-only / no van -> S
    for lig in ("vanA", "vanB", "vanD", "vanM"):
        assert call_efm_vancomycin([lig])["prediction"] == "R", lig
    assert call_efm_vancomycin([])["prediction"] == "S"
    assert call_efm_vancomycin(["tet(M)"])["prediction"] == "S"
    # accessory genes WITHOUT the bare ligase must NOT trigger R (the ligase is the determinant)
    assert call_efm_vancomycin(["vanH-A", "vanX-A", "vanR-A"])["prediction"] == "S"


def test_teicoplanin_vanA_vs_vanB_split():
    # vanA/vanD/vanM -> teico-R; vanB/vanC alone -> teico-S (the classic glycopeptide split)
    assert call_efm_teicoplanin(["vanA"])["prediction"] == "R"
    assert call_efm_teicoplanin(["vanD"])["prediction"] == "R"
    assert call_efm_teicoplanin(["vanB"])["prediction"] == "S"       # vanB = teicoplanin-SUSCEPTIBLE
    assert call_efm_teicoplanin(["vanC"])["prediction"] == "S"
    assert call_efm_teicoplanin([])["prediction"] == "S"


def test_van_vve_expression_disclosure():
    # structural van genes present but vanR/vanS regulator ABSENT -> VVE non-expression flag (disclosed,
    # NOT a call change: the genotype stays resistance-capable -> R). This is the SAMN11953784 case.
    vve = ["vanA", "vanH-A", "vanX-A", "vanY-A", "vanZ-A"]   # full structural, NO vanR-A/vanS-A
    v = call_efm_vancomycin(vve)
    assert v["prediction"] == "R"                            # genotype resistance-capable -> still R
    assert v["expression_note"] and "VVE" in v["expression_note"]
    # with the regulators present, no VVE flag
    full = vve + ["vanR-A", "vanS-A"]
    assert call_efm_vancomycin(full)["expression_note"] is None
