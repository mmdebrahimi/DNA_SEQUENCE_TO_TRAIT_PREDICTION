"""Offline test for the cross-kingdom summary builder (preserves distinct tiers; no aggregate flatten)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.build_cross_kingdom_summary import build, render_md  # noqa: E402


def test_no_aggregate_headline_and_distinct_tiers():
    card = build()
    assert card["schema"] == "cross-kingdom-validation-summary-v1"
    assert card["no_aggregate_headline"] is True
    # each surface carries its OWN independence string (not a shared/averaged one)
    indep = [r["independence"] for r in card["surfaces"]]
    assert len(set(indep)) == len(indep) or len(indep) <= 1   # distinct per surface
    md = render_md(card)
    assert "No aggregate headline" in md and "category error" in md.lower()


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
