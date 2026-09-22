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


def test_the_in_distribution_rows_disclose_their_quality_filter():
    """`n` on the IN-DISTRIBUTION rows is a HIGH-quality SUBSET (rifampicin 8,955 of 12,097 labelled),
    because score_tb_cryptic_parquet.py filters PHENOTYPE_QUALITY == HIGH. Undisclosed, a reader takes
    `n` for the compendium. The card's numbers are unchanged -- the filter was always applied; what was
    missing was saying so."""
    md = render_md(build())
    assert "HIGH-QUALITY-ONLY" in md
    assert "8,955" in md and "12,097" in md


def test_the_quality_disclosure_says_the_excluded_tiers_are_not_exchangeable():
    """The load-bearing half. 'We dropped the noisy ones' implies a random sample; R-prevalence differs
    sharply by tier (RIF HIGH 0.385 / MEDIUM 0.590 / LOW 0.244), so the discarded isolates are a
    DIFFERENT population rather than a noisier copy of the same one."""
    md = render_md(build())
    assert "not " in md.lower() and "exchangeable" in md
    assert "0.385" in md and "0.590" in md and "0.244" in md
