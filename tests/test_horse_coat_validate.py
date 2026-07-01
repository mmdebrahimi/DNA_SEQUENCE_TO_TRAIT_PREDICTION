"""Pin the horse-coat validation harness on synthetic OBSERVED data (no network) + the data-wall honesty."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.horse_coat_validate import run, score_rows  # noqa: E402


def test_scores_observed_rows_and_surfaces_discordant():
    rows = [
        {"mc1r": "ee", "asip": "AA", "observed_colour": "chestnut"},  # correct
        {"mc1r": "EE", "asip": "aa", "observed_colour": "black"},     # correct
        {"mc1r": "Ee", "asip": "Aa", "observed_colour": "bay"},       # correct
        {"mc1r": "EE", "asip": "aa", "observed_colour": "bay"},       # DISCORDANT (rule-breaker)
        {"mc1r": "Ee", "asip": "Aa", "observed_colour": "grey"},      # excluded (non-base)
    ]
    r = score_rows(rows)
    assert r["n_scored"] == 4 and r["n_correct"] == 3
    assert r["concordance"] == 0.75
    assert r["n_excluded_nonbase"] == 1
    assert r["confusion_observed_to_predicted"]["black"]["black"] == 1
    assert r["confusion_observed_to_predicted"]["bay"]["black"] == 1   # the rule-breaker is visible


def test_data_wall_when_absent(tmp_path):
    res = run(tmp_path / "nope.tsv")
    assert res["status"] == "VALIDATION_DATA_WALL"
    assert "circular" in res["note"]
