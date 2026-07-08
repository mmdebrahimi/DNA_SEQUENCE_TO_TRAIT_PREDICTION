"""Offline pins for the J3-ABO substrate parsers (no zip / no network)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import abo_opensnp_ingest as m  # noqa: E402


def test_parse_abo_group_common_forms():
    assert m.parse_abo_group("O+") == "O"
    assert m.parse_abo_group("0-") == "O"          # people write zero for O
    assert m.parse_abo_group("A+") == "A"
    assert m.parse_abo_group("A-") == "A"
    assert m.parse_abo_group("B+") == "B"
    assert m.parse_abo_group("Ab+") == "AB"
    assert m.parse_abo_group("AB-") == "AB"
    assert m.parse_abo_group("b+/-") == "B"


def test_parse_abo_group_genotype_style():
    assert m.parse_abo_group("AO - -") == "A"      # AO genotype -> A phenotype
    assert m.parse_abo_group("BO") == "B"
    assert m.parse_abo_group("OO") == "O"
    assert m.parse_abo_group("AA") == "A"


def test_parse_abo_group_rejects_non_abo():
    assert m.parse_abo_group("Rh negative") is None   # pure Rh, no ABO letter
    assert m.parse_abo_group("rather not say") is None
    assert m.parse_abo_group("") is None
    assert m.parse_abo_group("-") is None
    assert m.parse_abo_group("Can't smell") is None   # decoy-column junk


def test_deterministic_o_call():
    assert m.deterministic_o_call("DD") == "O"        # homozygous deletion -> type O
    assert m.deterministic_o_call("DI") == "non_O"
    assert m.deterministic_o_call("II") == "non_O"
    assert m.deterministic_o_call("ID") == "non_O"
    assert m.deterministic_o_call(None) is None
    assert m.deterministic_o_call("GG") is None        # non-D/I form -> abstain


def test_abo_snp_set():
    assert m.ABO_SNPS == {"rs8176719", "rs8176746", "rs8176747"}
