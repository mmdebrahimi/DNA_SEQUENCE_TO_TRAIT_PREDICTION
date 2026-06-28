"""Smoke + structure tests for the HIV report-card roll-up (Rec 3)."""
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


if __name__ == "__main__":
    test_latest_missing_prefix_returns_none()
    test_build_structure_and_modality_separation()
    test_pi_insti_cai_cells_carry_ols_baseline()
    test_insti_v0_1_gain_wired_when_artifact_present()
    print("PASS")
