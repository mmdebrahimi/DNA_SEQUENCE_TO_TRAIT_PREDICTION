"""Pins the inline trust-surface: every decoder cell resolves to its HONEST validation tier, no tier is
fabricated, and the independence flag is consistent. Card-dependent tiers skip if a card JSON is absent;
the structural invariants (always-a-dict / no-fabrication / genus-normalization / UNKNOWN) always run.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.data import trust_surface as ts  # noqa: E402

_WIKI = REPO / "wiki"
_REQUIRED_KEYS = {"tier", "independent", "headline", "metric", "n", "cell", "source_card", "caveat"}


def _has(card: str) -> bool:
    return (_WIKI / card).exists()


# --- structural invariants (always run) ---

def test_always_returns_full_dict():
    for drug, org in [("ciprofloxacin", "Escherichia"), ("zzz_nonsense", None), ("rifampicin", None)]:
        b = ts.trust_block(drug, org)
        assert _REQUIRED_KEYS <= set(b), f"missing keys for {drug}"
        assert isinstance(b["caveat"], str) and b["caveat"]


def test_unknown_drug_is_unknown_and_never_fabricates():
    b = ts.lookup_trust("not_a_real_drug_xyz", "Nowhere")
    assert b["tier"] == ts.UNKNOWN
    assert b["metric"] is None and b["independent"] is False


def test_independent_flag_is_consistent():
    # only the two independent tiers may carry independent=True
    for drug, org in [("efavirenz", None), ("rifampicin", None), ("ciprofloxacin", "Escherichia"),
                      ("nirmatrelvir", None), ("fluconazole", "Candida_auris"), ("oxacillin", "Staphylococcus_aureus")]:
        b = ts.trust_block(drug, org)
        assert b["independent"] == (b["tier"] in (ts.INDEPENDENT_WETLAB, ts.INDEPENDENT_MEASURED))


def test_genus_normalization_collapses_organism_spellings():
    a = ts.trust_block("ciprofloxacin", "Escherichia")
    b = ts.trust_block("ciprofloxacin", "Escherichia_coli_Shigella")
    assert a["tier"] == b["tier"] and a["cell"] == b["cell"]


def test_one_line_is_ascii_safe():
    s = ts.one_line(ts.trust_block("rifampicin"))
    assert isinstance(s, str) and "validation:" in s
    assert "—" not in s  # no em-dash (cp1252-console trap)


# --- per-tier pins (skip if the backing card is absent) ---

@pytest.mark.skipif(not _has("hiv_decoder_report_card.json"), reason="hiv card absent")
def test_hiv_is_free_wetlab_independent():
    b = ts.trust_block("efavirenz")
    assert b["tier"] == ts.INDEPENDENT_WETLAB and b["independent"] is True
    assert b["metric"] is not None and "wetlab" in b["caveat"].lower() or "wet-lab" in b["caveat"].lower()


@pytest.mark.skipif(not _has("tb_report_card.json"), reason="tb card absent")
def test_tb_is_independent_measured():
    b = ts.trust_block("rifampicin", "Mycobacterium_tuberculosis")
    assert b["tier"] == ts.INDEPENDENT_MEASURED and b["independent"] is True
    assert b["source_card"].endswith("tb_report_card.md")


@pytest.mark.skipif(not _has("amr_portal_independent_report_card.json"), reason="amr portal card absent")
def test_ecoli_cipro_is_independent_measured():
    b = ts.trust_block("ciprofloxacin", "Escherichia")
    assert b["tier"] == ts.INDEPENDENT_MEASURED and b["independent"] is True
    assert b["metric"] is not None


def test_sarscov2_is_in_distribution():
    b = ts.trust_block("nirmatrelvir")
    assert b["tier"] == ts.IN_DISTRIBUTION and b["independent"] is False


def test_fungal_is_no_free_source():
    b = ts.trust_block("fluconazole", "Candida_auris")
    assert b["tier"] == ts.NO_FREE_PHENOTYPE_SOURCE and b["independent"] is False


def test_oxacillin_is_label_confounded():
    b = ts.trust_block("oxacillin", "Staphylococcus_aureus")
    assert b["tier"] == ts.LABEL_CONFOUNDED


def test_meropenem_acinetobacter_abstains():
    b = ts.trust_block("meropenem", "Acinetobacter")
    assert b["tier"] == ts.ABSTAINS_BY_DESIGN


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        try:
            fn(); print(f"PASS {fn.__name__}")
        except Exception as e:  # pragma: no cover
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns)} tests")
