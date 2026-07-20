"""Pin the AR-Bank organism registry: uniform rule_fn dispatch + drug-map shape."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.organism_rules import ar_bank_registry as reg  # noqa: E402


def test_enterococcus_faecium_dispatch():
    fn = reg.rule_fn_for("enterococcus_faecium")
    # Levofloxacin <- cipro gyrA/parC QRDR rule (E. faecium codons Ser83/Glu87 gyrA, Ser80/Glu84 parC)
    assert fn("levofloxacin", ["gyrA_S83I"])["prediction"] == "R"
    assert fn("levofloxacin", [])["prediction"] == "S"
    # Doxycycline <- acquired tet gene rule
    assert fn("doxycycline", ["tet(M)"])["prediction"] == "R"
    assert fn("doxycycline", ["gyrA_S83I"])["prediction"] == "S"    # FQ marker, not tet
    # unsupported drug -> abstain
    assert fn("daptomycin", [])["prediction"] == "INDETERMINATE"


def test_gono_surfaced_through_registry():
    fn = reg.rule_fn_for("gono")
    assert fn("ciprofloxacin", ["gyrA_S91F"])["prediction"] == "R"
    assert fn("ceftriaxone", ["penA_A501P"])["prediction"] == "R"   # v0.1 A501-class
    assert fn("ceftriaxone", ["penA_G545S"])["prediction"] == "S"   # v0.1 non-A501 mosaic


def test_staph_only_levofloxacin():
    cfg = reg.config_for("staphylococcus_aureus")
    assert set(cfg["drug_map"]) == {"levofloxacin"}                  # no Rifampin/cipro label on the panel
    fn = reg.rule_fn_for("staphylococcus_aureus")
    assert fn("rifampicin", [])["prediction"] == "INDETERMINATE"     # not registered -> abstain


def test_shown_name_map_and_unknown_key():
    import pytest
    assert reg.shown_name_map("enterococcus_faecium")["levofloxacin"] == "Levofloxacin"
    with pytest.raises(KeyError):
        reg.config_for("nonexistent_organism")
