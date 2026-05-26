from scripts.drug_mechanism_phenotype_merge import classify_noise


def test_cef_clean_r_primary_mechanism():
    noise, opacity, co_res = classify_noise(
        {
            "cohort_binary_label": 1,
            "mic_tier": "HIGH_R",
            "mechanisms_present": ["acquired_beta_lactamase"],
        },
        "ceftriaxone",
    )
    assert noise == "CLEAN_R_primary_mechanism"
    assert opacity is False
    assert co_res == []


def test_cef_opaque_r_co_resistance_only():
    noise, opacity, co_res = classify_noise(
        {
            "cohort_binary_label": 1,
            "mic_tier": "HIGH_R",
            "mechanisms_present": ["porin_loss"],
        },
        "ceftriaxone",
    )
    assert noise == "OPAQUE_R_co_resistance_only"
    assert opacity is True
    assert co_res == ["porin_loss"]


def test_cef_suspect_s_silent_primary_mechanism():
    noise, opacity, _ = classify_noise(
        {
            "cohort_binary_label": 0,
            "mic_tier": "HIGH_S",
            "mechanisms_present": ["ampC_hyperproduction"],
        },
        "ceftriaxone",
    )
    assert noise == "SUSPECT_S_silent_primary_mechanism"
    assert opacity is False
