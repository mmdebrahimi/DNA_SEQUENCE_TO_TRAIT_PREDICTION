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
        assert c["drug_class"] in ("NNRTI", "NRTI")
        assert "auc_call_separates_fold" in c and "delta_ols_minus_catalog" in c


if __name__ == "__main__":
    test_latest_missing_prefix_returns_none()
    test_build_structure_and_modality_separation()
    print("PASS")
