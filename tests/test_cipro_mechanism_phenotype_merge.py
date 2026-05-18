"""Tests for scripts/cipro_mechanism_phenotype_merge.py — pure-logic functions.

Pins the per-row noise classification + strict primary mechanism definition +
opacity flag separation that drives the EP1 SUSPEND_CONDITION_4 gate. The
JSON-load + merge + packet-write paths are orchestration (skipped).
"""
from __future__ import annotations

import pytest

from scripts.cipro_mechanism_phenotype_merge import (
    CO_RESISTANCE_MECHANISMS,
    PRIMARY_CIPRO_MECHANISMS,
    _classify_noise,
)


def test_primary_cipro_mechanisms_strict_definition():
    # Strict primary = QRDR OR plasmid only. Efflux/regulatory/porin are
    # co-resistance modifiers, NOT primary.
    assert PRIMARY_CIPRO_MECHANISMS == {"QRDR_target_alteration", "plasmid_protect_modify"}


def test_co_resistance_mechanisms_separated_from_primary():
    # Efflux + regulatory + porin_loss are co-resistance only
    assert "efflux" in CO_RESISTANCE_MECHANISMS
    assert "regulatory" in CO_RESISTANCE_MECHANISMS
    assert "porin_loss" in CO_RESISTANCE_MECHANISMS
    assert PRIMARY_CIPRO_MECHANISMS.isdisjoint(CO_RESISTANCE_MECHANISMS)


# ---- R-strain classification ------------------------------------------------


def test_classify_clean_r_primary_mechanism():
    # HIGH_R + QRDR -> CLEAN_R_primary_mechanism
    noise, opacity, co_res = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "HIGH_R",
        "mechanisms_present": ["QRDR_target_alteration"],
    })
    assert noise == "CLEAN_R_primary_mechanism"
    assert opacity is False
    assert co_res == []


def test_classify_opaque_r_co_resistance_only():
    # HIGH_R + only efflux/regulatory -> OPAQUE_R_co_resistance_only + opacity flag
    noise, opacity, co_res = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "HIGH_R",
        "mechanisms_present": ["efflux"],
    })
    assert noise == "OPAQUE_R_co_resistance_only"
    assert opacity is True  # tool-incomplete signal
    assert co_res == ["efflux"]


def test_classify_opaque_r_no_mechanism():
    # HIGH_R + no AMRFinder hit -> OPAQUE_R_no_mechanism + opacity flag
    noise, opacity, _ = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "HIGH_R",
        "mechanisms_present": [],
    })
    assert noise == "OPAQUE_R_no_mechanism"
    assert opacity is True


def test_classify_noisy_r_borderline():
    noise, opacity, _ = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "BORDERLINE",
        "mechanisms_present": ["QRDR_target_alteration"],
    })
    assert noise == "NOISY_R_borderline"
    assert opacity is False  # noisy label, not tool gap


def test_classify_noisy_r_no_mic():
    noise, _, _ = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "NO_MIC",
        "mechanisms_present": [],
    })
    assert noise == "NOISY_R_no_mic"


# ---- S-strain classification ------------------------------------------------


def test_classify_clean_s_no_primary_mechanism():
    # HIGH_S + no primary mech -> CLEAN_S
    noise, opacity, _ = _classify_noise({
        "cohort_binary_label": 0,
        "mic_tier": "HIGH_S",
        "mechanisms_present": [],
    })
    assert noise == "CLEAN_S_no_primary_mechanism"
    assert opacity is False


def test_classify_suspect_s_silent_primary_mechanism():
    # HIGH_S + QRDR present -> SUSPECT (silent mechanism)
    noise, _, _ = _classify_noise({
        "cohort_binary_label": 0,
        "mic_tier": "HIGH_S",
        "mechanisms_present": ["QRDR_target_alteration"],
    })
    assert noise == "SUSPECT_S_silent_primary_mechanism"


def test_classify_suspect_s_borderline_primary_mechanism():
    # Borderline MIC + primary mech present -> likely mislabeled to S
    noise, _, _ = _classify_noise({
        "cohort_binary_label": 0,
        "mic_tier": "BORDERLINE",
        "mechanisms_present": ["plasmid_protect_modify"],
    })
    assert noise == "SUSPECT_S_borderline_primary_mechanism"


def test_classify_noisy_s_borderline():
    # Borderline MIC + no primary mech -> NOISY_S_borderline
    noise, _, _ = _classify_noise({
        "cohort_binary_label": 0,
        "mic_tier": "BORDERLINE",
        "mechanisms_present": ["efflux"],  # co-resistance, not primary
    })
    assert noise == "NOISY_S_borderline"


def test_classify_noisy_s_no_mic():
    noise, _, _ = _classify_noise({
        "cohort_binary_label": 0,
        "mic_tier": "NO_MIC",
        "mechanisms_present": [],
    })
    assert noise == "NOISY_S_no_mic"


# ---- Co-resistance modifier tracking ----------------------------------------


def test_co_resistance_modifiers_returned_separately():
    # Strain with QRDR + efflux + regulatory: primary = QRDR; co_res = efflux + regulatory
    _, _, co_res = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "HIGH_R",
        "mechanisms_present": ["QRDR_target_alteration", "efflux", "regulatory"],
    })
    assert "efflux" in co_res
    assert "regulatory" in co_res
    # Primary mechanism NOT in co_res
    assert "QRDR_target_alteration" not in co_res


def test_co_resistance_empty_when_only_primary():
    _, _, co_res = _classify_noise({
        "cohort_binary_label": 1,
        "mic_tier": "HIGH_R",
        "mechanisms_present": ["QRDR_target_alteration"],
    })
    assert co_res == []


# ---- Opacity flag separates tool-incomplete from labels-noisy ---------------


def test_opacity_flag_true_only_when_high_r_no_primary():
    # The opacity flag is specifically for "HIGH_R but no primary mech" cases —
    # AMRFinder might be incomplete, NOT that the label is wrong
    _, opacity1, _ = _classify_noise({
        "cohort_binary_label": 1, "mic_tier": "HIGH_R", "mechanisms_present": [],
    })
    assert opacity1 is True

    _, opacity2, _ = _classify_noise({
        "cohort_binary_label": 1, "mic_tier": "HIGH_R", "mechanisms_present": ["efflux"],
    })
    assert opacity2 is True  # co-resistance only -> tool may have missed primary

    _, opacity3, _ = _classify_noise({
        "cohort_binary_label": 1, "mic_tier": "BORDERLINE", "mechanisms_present": [],
    })
    assert opacity3 is False  # borderline label = labels noisy, NOT tool incomplete

    _, opacity4, _ = _classify_noise({
        "cohort_binary_label": 1, "mic_tier": "HIGH_R",
        "mechanisms_present": ["QRDR_target_alteration"],
    })
    assert opacity4 is False  # CLEAN; no tool-incomplete signal
