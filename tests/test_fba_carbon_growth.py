"""FBA carbon-source growth: pure name/exchange helpers + real-model medium-swap smoke (slow)."""
from __future__ import annotations

import pytest

from dna_decode.fba.carbon_growth import (
    match_carbon_exchange,
    normalize_carbon_name,
)


def test_normalize_carbon_name_alias():
    assert normalize_carbon_name("D-Maltose") == "maltose"
    assert normalize_carbon_name("a-Ketoglutaric") == "2-oxoglutarate"
    # unaliased names pass through lowercased
    assert normalize_carbon_name("D-Glucose") == "d-glucose"


def test_match_carbon_exchange():
    idx = {"d-glucose": "EX_glc__D_e", "maltose": "EX_malt_e"}
    assert match_carbon_exchange("D-Glucose", idx) == "EX_glc__D_e"       # exact (lowercased)
    assert match_carbon_exchange("D-Maltose", idx) == "EX_malt_e"          # via alias -> maltose
    assert match_carbon_exchange("Unobtainium", idx) is None               # unresolved


@pytest.mark.slow
def test_predict_growth_real_iml1515():
    pytest.importorskip("cobra")
    from dna_decode.fba.carbon_growth import build_exchange_name_index, predict_growth
    from dna_decode.fba.model import load_model

    m = load_model()
    idx = build_exchange_name_index(m)
    # glucose grows (~0.877); a pentose (arabinose) grows but slower than glucose
    glc = predict_growth(m, idx["d-glucose"])
    ara = predict_growth(m, idx["l-arabinose"])
    assert glc > 0.5
    assert ara > 1e-4 and ara < glc                    # pentoses grow slower -> quantitative + ordered
    # a carbon source with no exchange -> 0.0 (no crash)
    assert predict_growth(m, "EX_not_a_real_metabolite_e") == 0.0


# ---- salt / hydrate / formula normalisation (Fitness Browser join, 2026-08-12) ----

def test_strip_salt_and_hydrate_handles_the_real_assay_labels():
    """Measured on the Fitness Browser's 28 Keio carbon sources: exact+alias matching resolved only 14,
    and all 13 misses were NAME-gaps whose exchanges exist in iML1515."""
    from dna_decode.fba.carbon_growth import strip_salt_and_hydrate as s

    assert s("Sodium succinate dibasic hexahydrate") == "succinate"
    assert s("Potassium acetate") == "acetate"
    assert s("Sodium pyruvate") == "pyruvate"
    assert s("D-Gluconic Acid sodium salt") == "d-gluconate"
    assert s("Glycolic Acid") == "glycolate"


def test_longest_token_first_or_disodium_leaves_a_stray_di():
    """THE bug this ordering fixes: replacing 'sodium salt' before 'disodium salt' left 'l-malate di',
    which matched nothing. Tokens are sorted longest-first."""
    from dna_decode.fba.carbon_growth import strip_salt_and_hydrate as s

    got = s("L-Malic acid disodium salt monohydrate")
    assert got == "l-malate", got
    assert "di" not in got.split()


def test_greek_letter_prefix_is_dropped():
    from dna_decode.fba.carbon_growth import strip_salt_and_hydrate as s

    assert s("a-Ketoglutaric acid disodium salt hydrate") == "2-oxoglutarate"


def test_stereochemistry_is_NEVER_stripped():
    """Load-bearing safety rule. D-lactate and L-lactate are DIFFERENT metabolites with different
    exchanges (EX_lac__D_e vs EX_lac__L_e); collapsing them would silently mis-assign the carbon source.
    This is why 'D-Trehalose dihydrate' stays unmapped rather than being force-matched to 'trehalose'."""
    from dna_decode.fba.carbon_growth import strip_salt_and_hydrate as s

    assert s("Sodium D-Lactate").startswith("d-")
    assert s("D-Trehalose dihydrate") == "d-trehalose"


def test_formula_suffix_stripping_is_conservative():
    """BiGG writes 'maltose c12h22o11'. Strip the formula token, but never a real word."""
    from dna_decode.fba.carbon_growth import _strip_formula_suffix as f

    assert f("maltose c12h22o11") == "maltose"
    assert f("glycolate c2h3o3") == "glycolate"
    assert f("d-glucose 6-phosphate") == "d-glucose 6-phosphate"   # 'phosphate' is a word, not a formula
    assert f("trehalose") == "trehalose"


def test_an_unmappable_label_returns_None_rather_than_guessing():
    """A racemic mixture and an amino-acid mixture are not single exchanges. Returning None keeps them
    VISIBLY unmapped instead of silently binding to the wrong metabolite."""
    from dna_decode.fba.carbon_growth import match_carbon_exchange

    idx = {"succinate": "EX_succ_e", "l-lactate": "EX_lac__L_e", "d-lactate": "EX_lac__D_e"}
    assert match_carbon_exchange("Sodium D,L-Lactate", idx) is None
    assert match_carbon_exchange("casamino acids", idx) is None
    assert match_carbon_exchange("Sodium succinate dibasic hexahydrate", idx) == "EX_succ_e"
