"""Metabolic carbon-utilization cell — deterministic, uptake-gated. Pins measured E. coli K-12 phenotypes
incl. the citrate anchor (the case a naive has-the-genes rule mis-calls). Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.metabolic.carbon_catalog import (  # noqa: E402
    CARBON_SOURCES, MetabolicInputError, call_carbon_utilization, genes_for, resolve_substrate,
    reference_integrity_ok)


def test_reference_integrity_guard():
    assert reference_integrity_ok() is True


def test_lactose_full_operon_positive():
    c = call_carbon_utilization("lactose", ["lacZ", "lacY"])
    assert c.capability == "utilizes" and c.confidence == "high"
    assert c.transporter_present and c.transporter_expressed


def test_lactose_enzyme_knockout_negative():
    # lacY present (can import) but lacZ absent (cannot cleave) -> classic Lac-
    c = call_carbon_utilization("lactose", ["lacY"])
    assert c.capability == "cannot_utilize" and c.confidence == "high"
    assert "lacZ" in c.enzymes_missing


def test_lactose_transporter_knockout_negative():
    # lacZ present but lacY absent (cannot import) -> Lac-
    c = call_carbon_utilization("lactose", ["lacZ"])
    assert c.capability == "cannot_utilize"
    assert c.transporter_present is False


def test_citrate_aerobic_anchor_negative():
    # THE ANCHOR: full cit operon present, aerobic -> Cit- because citT is anaerobic-only.
    c = call_carbon_utilization("citrate", ["citD", "citE", "citF", "citT"], condition="aerobic")
    assert c.capability == "cannot_utilize"
    assert c.transporter_present is True          # the gene IS there
    assert c.transporter_expressed is False       # but not expressed aerobically
    assert c.confidence == "medium"
    assert any("naive" in n for n in c.notes)     # names the trap explicitly


def test_citrate_anaerobic_positive():
    c = call_carbon_utilization("citrate", ["citD", "citE", "citF", "citT"], condition="anaerobic")
    assert c.capability == "utilizes"
    assert c.transporter_expressed is True


def test_naive_has_genes_rule_would_disagree_on_citrate():
    # A naive "all pathway genes present -> utilizes" rule sees the full cit operon and says +.
    genes = ["citD", "citE", "citF", "citT"]
    naive_positive = all(g in genes for g in genes_for("citrate"))
    real = call_carbon_utilization("citrate", genes, condition="aerobic")
    assert naive_positive is True                        # naive rule: +
    assert real.capability == "cannot_utilize"           # truth: - (the uptake gate)


def test_arabinose_either_transporter():
    # araE OR araFGH suffices
    c_low = call_carbon_utilization("L-arabinose", ["araA", "araB", "araD", "araE"])
    c_abc = call_carbon_utilization("L-arabinose", ["araA", "araB", "araD", "araF", "araG", "araH"])
    assert c_low.capability == "utilizes" and c_abc.capability == "utilizes"
    # partial ABC transporter -> not present
    c_bad = call_carbon_utilization("L-arabinose", ["araA", "araB", "araD", "araF"])
    assert c_bad.capability == "cannot_utilize" and c_bad.transporter_present is False


def test_glucose_positive():
    assert call_carbon_utilization("D-glucose", ["pgi", "ptsG"]).capability == "utilizes"


def test_alias_and_case_resolution():
    assert resolve_substrate("lac") == "lactose"
    assert resolve_substrate("ARABINOSE") == "L-arabinose"
    assert resolve_substrate("cit") == "citrate"
    # gene matching is case-insensitive
    assert call_carbon_utilization("lactose", ["LACZ", "LACY"]).capability == "utilizes"


def test_unknown_substrate_raises():
    try:
        call_carbon_utilization("plutonium", ["x"])
    except MetabolicInputError:
        return
    raise AssertionError("expected MetabolicInputError on unknown substrate")


def test_unknown_condition_raises():
    try:
        call_carbon_utilization("lactose", ["lacZ", "lacY"], condition="martian")
    except MetabolicInputError:
        return
    raise AssertionError("expected MetabolicInputError on unknown condition")


def test_as_dict_shape():
    d = call_carbon_utilization("citrate", ["citD", "citE", "citF", "citT"], condition="aerobic").as_dict()
    assert d["organism"] == "Escherichia_coli" and d["trait"] == "carbon_utilization"
    assert d["capability"] == "cannot_utilize" and "undetectable_mechanisms" in d


def test_catalog_transporters_nonempty():
    for name, ent in CARBON_SOURCES.items():
        assert ent.enzymes, name
        assert ent.transporters, name
        assert ent.transporter_expressed, name


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
