"""Guards for pointing the project's own source-concentration bar at the AR Bank arm."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import ar_bank_source_concentration as m  # noqa: E402

ARTIFACT = ROOT / "wiki" / "ar_bank_source_concentration_2026-09-10.json"


def test_the_bar_is_imported_not_restated():
    """Restating 0.60 here would make this a DIFFERENT standard wearing the same name -- the point is
    to judge our own cohorts by the number the frozen arm is judged by."""
    from scripts.source_diverse_validate import MAX_SOURCE_SHARE
    assert m.MAX_SOURCE_SHARE == MAX_SOURCE_SHARE


def test_concentration_is_computed_over_resolved_rows_only():
    panels = {"A": 1, "B": 1, "C": 2}
    rows = [{"biosample": x} for x in ("A", "B", "C", "UNKNOWN")]
    c = m.concentration(rows, panels)
    assert c["n_rows"] == 4 and c["n_resolved"] == 3, "an unresolvable row must not enter the denominator"
    # the reporter rounds to 4dp, so compare against the rounded value rather than the exact ratio
    assert c["n_panels"] == 2 and c["largest_share"] == pytest.approx(round(2 / 3, 4))


def test_a_cohort_with_no_resolvable_panel_returns_none_not_a_clean_pass():
    """Returning a passing verdict for a cohort we could not group would be the worst outcome: it would
    certify diversity that was never measured."""
    assert m.concentration([{"biosample": "X"}], {}) is None


def test_a_single_panel_cohort_fails_the_bar():
    c = m.concentration([{"biosample": "A"}, {"biosample": "B"}], {"A": 7, "B": 7})
    assert c["largest_share"] == 1.0 and c["passes_bar"] is False


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not generated")
def test_the_finding_is_confined_not_blanket():
    """The hypothesis was that EVERY AR Bank cell would fail. Measurement refutes it: the gonorrhoeae
    arm is single-panel while Klebsiella and E. coli draw on 5-6 panels and pass. A test that only
    asserted 'some fail' would have let the over-claim stand."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    cells = {c["cohort"]: c for c in d["cells"]}
    gono = [c for k, c in cells.items() if "_gono_" in k]
    assert gono and all(c["largest_share"] == 1.0 and not c["passes_bar"] for c in gono)
    others = [c for k, c in cells.items() if "_kleb_" in k or "_ecoli_" in k]
    assert others and all(c["passes_bar"] for c in others), (
        "if these now fail too, the finding is no longer 'confined' and the headline must be rewritten")
    assert 0 < d["n_failing_bar"] < d["n_cells"], "a blanket verdict in either direction is the wrong shape"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not generated")
def test_the_proxy_caveat_ships_with_the_numbers():
    """Panel is not BioProject and the bar was calibrated on BioProjects. A reader who sees the verdict
    without that sentence has been told a like-for-like number that this is not."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["grouping"] == "cdc_panel_id"
    assert "not a like-for-like BioProject number" in d["grouping_is_a_proxy"]
    assert d["honest_limits"]
