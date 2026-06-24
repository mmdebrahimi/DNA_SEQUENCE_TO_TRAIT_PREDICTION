"""Offline test for the TB report-card builder pure logic (tier rows + headline-rule)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.build_tb_report_card import build, render_md  # noqa: E402


def test_build_has_namespace_separation_and_headline_rule():
    card = build()
    assert card["schema"] == "tb-report-card-v1"
    assert "homoplasic" in card["headline_rule"] or "homoplas" in card["headline_rule"].lower()
    # independent tier present + raw is the headline (each indep row carries raw_sens + a lineage DISCLOSURE)
    if card["independent"]:
        r = card["independent"][0]
        assert "raw_sens" in r and "lineage_sens_disclosure" in r
        assert r["tier"] == "PROVENANCE_DISJOINT_INDEPENDENT"
    md = render_md(card)
    assert "namespace-separate" in md.lower() and "RAW per-isolate" in md


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
