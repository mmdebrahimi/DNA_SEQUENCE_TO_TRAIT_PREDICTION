"""Smoke + structure tests for the HIV report-card roll-up (Rec 3)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.build_hiv_report_card import _latest, build  # noqa: E402


def test_latest_missing_prefix_returns_none():
    assert _latest("definitely_no_such_artifact_prefix_zzz_") is None


def test_build_structure_and_modality_separation():
    rc = build()
    assert rc["artifact"] == "hiv_decoder_report_card"
    assert isinstance(rc["cells"], list) and rc["n_cells"] == len(rc["cells"])
    # the load-bearing honesty: explicitly NOT conflated with the bacterial provenance-disjoint card
    assert "provenance-disjoint" in rc["modality"] and "NOT conflated" in rc["modality"]
    assert "Sierra" in rc["label_independence"]  # circularity-safe label note carried
    # every cell carries class + drug + a label-independent metric slot
    for c in rc["cells"]:
        assert c["drug_class"] in ("NNRTI", "NRTI", "PI", "INSTI", "CAI")
        assert "auc_call_separates_fold" in c and "delta_ols_minus_catalog" in c


def test_pi_insti_cai_cells_carry_ols_baseline():
    # The wrapper-vs-tool discipline gap closed 2026-06-22: PI/INSTI/CAI cells must carry the OLS
    # underlying-tool baseline (not just the catalog AUC). Relies on the committed baseline JSONs.
    rc = build()
    new_cells = [c for c in rc["cells"] if c["drug_class"] in ("PI", "INSTI", "CAI")]
    assert new_cells, "expected PI/INSTI/CAI cells (committed validation JSONs present)"
    scored = [c for c in new_cells if c["ols_baseline_balacc"] is not None]
    # at least the well-powered PI/INSTI/CAI drugs carry a numeric OLS baseline + catalog balacc + delta
    assert scored, "PI/INSTI/CAI cells must carry the OLS baseline once baseline JSONs exist"
    for c in scored:
        assert c["catalog_balacc"] is not None and c["delta_ols_minus_catalog"] is not None


def test_insti_v0_1_gain_wired_when_artifact_present():
    """INSTI v0.1 (2026-06-27) deconfounded gains must surface in the card (not silently render '-')."""
    rc = build()
    insti = [c for c in rc["cells"] if c["drug_class"] == "INSTI"]
    if not insti:
        return  # the committed hiv_insti_v0_validation_ artifact is absent on this checkout
    gains = [c.get("v0_1_mutant_gain") for c in insti]
    assert any(isinstance(g, (int, float)) for g in gains), (
        "INSTI v0.1 gains not wired — check hiv_insti_v0.1_validation_ ingestion in build_hiv_report_card"
    )


def test_the_censoring_of_the_underlying_label_is_disclosed():
    """The PhenoSense fold tables are right-censored at 100 (3TC 45.2% / NVP 33.0% AT the cap, both
    failing this repo's own assay_degeneracy bar). The `v0.1 gain` column is DOWNSTREAM of that: those
    catalogs were selected by thresholding an OLS coefficient fit on the censored response. A reader of
    that column must be able to see its provenance from the card itself."""
    rc = build()
    cav = " ".join(rc["honest_caveats"])
    assert "RIGHT-CENSORED AT 100" in cav
    assert "v0.1 gain" in cav


def test_the_censoring_caveat_says_what_is_NOT_affected():
    """Over-reading this as 'the HIV numbers are wrong' would be the opposite error. The AUC and both
    balacc columns are computed at a fold cutoff of 3; the ceiling is 100, so every censored observation
    is R under either reading and the R/S labels are untouched. The caveat must scope itself."""
    rc = build()
    cav = " ".join(rc["honest_caveats"])
    assert "are NOT" in cav and "cutoff of" in cav


def test_the_censoring_caveat_does_not_claim_the_unclaimed_refit():
    """The Tobit refit's own frozen verdict is IMPLEMENTATION_SUSPECT, so the card may CITE it but must
    not present it as an established correction, and must restate no number."""
    rc = build()
    cav = " ".join(rc["honest_caveats"])
    assert "IMPLEMENTATION_SUSPECT" in cav and "NOT claimed" in cav


def test_the_disclosure_is_augment_only_and_changes_no_cell():
    """Every disclosure layer in this project AUGMENTS and never restates a measurement. Pinned so a
    future edit cannot quietly move a published number under cover of adding a caveat."""
    rc = build()
    assert len(rc["cells"]) == 25
    for c in rc["cells"]:
        assert "censor" not in json.dumps(c).lower(), c.get("drug")


if __name__ == "__main__":
    test_latest_missing_prefix_returns_none()
    test_build_structure_and_modality_separation()
    test_pi_insti_cai_cells_carry_ols_baseline()
    test_insti_v0_1_gain_wired_when_artifact_present()
    print("PASS")
