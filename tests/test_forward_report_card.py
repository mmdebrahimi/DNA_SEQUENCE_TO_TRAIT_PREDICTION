"""Guard: the forward report card builder runs read-only + renders every capability with no None-leak."""
import json, subprocess, sys
from pathlib import Path


def test_report_card_builds_and_is_honest(tmp_path=None):
    r = subprocess.run([sys.executable, "scripts/build_forward_report_card.py"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    card = json.loads(Path("wiki/forward_validation_report_card.json").read_text(encoding="utf-8"))
    assert card["schema"] == "forward-validation-report-card-v1"
    caps = card["capabilities"]
    assert len(caps) >= 7
    # no aggregate headline field; every capability has a tier + scope
    assert "aggregate" not in card and "headline" not in card
    for c in caps:
        assert c["tier"] and c["scope"] and c["metric"]
        assert "None" not in c["metric"]           # no unpopulated source leaked into a claim
    # the multi-mutant + epistasis rows (tonight's work) are present
    names = " ".join(c["capability"] for c in caps)
    assert "additive-null" in names and "epistasis" in names and "regime router" in names


if __name__ == "__main__":
    test_report_card_builds_and_is_honest(); print("PASS")
