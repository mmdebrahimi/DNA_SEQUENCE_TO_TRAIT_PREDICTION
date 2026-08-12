"""Label-matched growth media (wheel-only; no cobra, no network)."""
from __future__ import annotations

import pytest

from dna_decode.fba.essentiality_labels import ESSENTIALITY_LABEL_CONDITION
from dna_decode.fba.medium import RICH_MEDIUM_EXCHANGES, apply_rich_medium, rich_medium


class _Rxn:
    def __init__(self, rid):
        self.id = rid


class _Model:
    """Minimal cobra-shaped stand-in: `.exchanges` and a settable `.medium` dict."""

    def __init__(self, exchange_ids, medium):
        self.exchanges = [_Rxn(i) for i in exchange_ids]
        self.medium = dict(medium)


def test_rich_medium_is_additive_and_keeps_the_carbon_source():
    """Load-bearing: supplements are ADDED to the reconstruction's own medium. Replacing it would drop the
    carbon source / oxygen bounds and change what is being measured."""
    m = _Model(["EX_glc__D_e", "EX_o2_e", "EX_leu__L_e"], {"EX_glc__D_e": 10.0, "EX_o2_e": 2.0})
    med = rich_medium(m)
    assert med["EX_glc__D_e"] == 10.0        # carbon source untouched
    assert med["EX_o2_e"] == 2.0             # oxygen bound untouched
    assert med["EX_leu__L_e"] == 10.0        # supplement opened


def test_supplements_absent_from_the_model_are_skipped_not_raised():
    """Portability across reconstructions: a model lacking an amino-acid exchange must not crash."""
    m = _Model(["EX_glc__D_e"], {"EX_glc__D_e": 10.0})
    assert rich_medium(m) == {"EX_glc__D_e": 10.0}
    assert apply_rich_medium(m) == []


def test_apply_returns_only_the_newly_opened_exchanges():
    m = _Model(["EX_glc__D_e", "EX_leu__L_e", "EX_ade_e"], {"EX_glc__D_e": 10.0, "EX_leu__L_e": 10.0})
    opened = apply_rich_medium(m)
    assert opened == ["EX_ade_e"]            # leu was ALREADY open -> not reported as newly opened


def test_rich_medium_covers_all_twenty_proteinogenic_amino_acids():
    """A partial supplement set would leave some biosynthesis genes spuriously essential -- the exact
    failure this module exists to remove."""
    aa = {"ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile",
          "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val"}
    present = {e.split("_")[1] for e in RICH_MEDIUM_EXCHANGES if e.startswith("EX_")}
    assert aa <= present


def test_uptake_rate_is_configurable():
    m = _Model(["EX_leu__L_e"], {})
    assert rich_medium(m, uptake=3.5)["EX_leu__L_e"] == 3.5


def test_sgd_yeast_labels_are_registered_as_rich():
    """SGD's inviable-null set comes from the deletion collection on YPD. Scoring it against a
    minimal-medium model charged the model for biology (MCC 0.2524 vs 0.3773)."""
    for alias in ("yeast", "scerevisiae", "saccharomyces_cerevisiae"):
        assert ESSENTIALITY_LABEL_CONDITION[alias] == "rich"


def test_an_unregistered_organism_falls_back_to_the_model_default():
    """No condition recorded => do NOT guess a medium; use the reconstruction's own."""
    assert ESSENTIALITY_LABEL_CONDITION.get("saureus", "default") == "default"
