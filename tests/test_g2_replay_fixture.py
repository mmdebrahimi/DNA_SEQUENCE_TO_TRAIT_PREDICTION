"""Verify the G2 Arabidopsis replay fixture against the markdown it was transcribed from.

WHY THIS EXISTS. The F1 falsifier's Arm C replays the committed Arabidopsis per-seed scores through
the verdict function. Those numbers lived ONLY inside a markdown table, so the alternative was to
hard-code twelve literals in the falsifier -- an assertion that can drift from its source with nothing
able to notice. The fixture is HAND-TRANSCRIBED; this module RE-PARSES the markdown and fails loudly
on any mismatch. Verified, never trusted -- the same discipline as the repo's contract-number bindings.

Offline, pure, no fixtures beyond the two committed files.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

MD = ROOT / "wiki" / "phase2_arabidopsis_result_2026-06-12.md"
JSON_PATH = ROOT / "wiki" / "phase2_arabidopsis_result_2026-06-12.json"

#: U+2212 MINUS SIGN, not ASCII hyphen. float("−0.173") raises; normalise before parsing.
MINUS = "−"

#: markdown row label -> fixture metric key. The labels carry U+00B2, bold markers and a suffix.
_LABEL_TO_KEY = {
    "global r2": "global_r2",
    "spearman": "spearman",
    "within-group r2": "within_group_r2",
}


def _norm_label(cell: str) -> str:
    """Strip bold markers, the '(de-confounded)' suffix and superscript-2 -> '2'."""
    s = cell.replace("**", "").replace("²", "2")
    s = re.sub(r"\(.*?\)", "", s)
    return s.strip().lower()


def _num(tok: str) -> float:
    """Parse one cell value, handling U+2212, a leading '+', and bold markers."""
    return float(tok.replace("**", "").replace(MINUS, "-").replace("+", "").strip())


def parse_markdown_table() -> dict:
    """Re-derive {metric: {embedding: [...], structure_only: float}} from the markdown."""
    text = MD.read_text(encoding="utf-8")
    body = text.split("## Results", 1)[1].split("## Interpretation", 1)[0]
    out: dict[str, dict] = {}
    for line in body.splitlines():
        if not line.strip().startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 4:
            continue
        key = _LABEL_TO_KEY.get(_norm_label(cells[0]))
        if key is None:                      # header row / separator / unknown metric
            continue
        out[key] = {
            "embedding": [_num(t) for t in cells[1].split("/")],
            "structure_only": _num(cells[2]),
            "winner_cell": cells[3].replace("**", "").strip(),
        }
    return out


def _fixture() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


# --- the transcription is verified, not trusted ---------------------------------------------------

def test_both_committed_files_exist_and_the_fixture_names_its_source():
    assert MD.exists() and JSON_PATH.exists()
    fx = _fixture()
    assert fx["source"] == "wiki/phase2_arabidopsis_result_2026-06-12.md"
    assert (ROOT / fx["source"]).exists(), "the fixture cites a source that does not resolve"


def test_the_markdown_parse_is_non_vacuous():
    """A parser that silently found nothing would make every comparison below pass trivially."""
    parsed = parse_markdown_table()
    assert set(parsed) == {"global_r2", "spearman", "within_group_r2"}, parsed
    for k, v in parsed.items():
        assert len(v["embedding"]) == 3, f"{k}: expected 3 seeds, got {v['embedding']}"


def test_every_fixture_value_matches_the_markdown_table():
    """THE GUARD. A mistyped transcription fails here."""
    parsed = parse_markdown_table()
    fx = _fixture()["metrics"]
    assert set(fx) == set(parsed)
    for key, want in parsed.items():
        got = fx[key]
        assert got["embedding"] == want["embedding"], f"{key} embedding: {got} vs md {want}"
        assert got["structure_only"] == want["structure_only"], \
            f"{key} structure_only: {got['structure_only']} vs md {want['structure_only']}"


def test_the_declared_winner_matches_the_markdown_winner_column():
    parsed = parse_markdown_table()
    fx = _fixture()["metrics"]
    for key, want in parsed.items():
        declared = fx[key]["winner"]
        assert declared.replace("_", "-") in want["winner_cell"], \
            f"{key}: fixture says {declared!r}, markdown winner cell is {want['winner_cell']!r}"


def test_the_verbatim_rows_are_actually_verbatim():
    """The fixture carries the raw rows so a human can eyeball the source without opening the md."""
    md = MD.read_text(encoding="utf-8")
    for row in _fixture()["verbatim_table_rows"]:
        assert row in md, f"not verbatim in the markdown: {row!r}"


# --- the STRUCTURAL fact Arm C turns on (not a transcription check) -------------------------------

def test_the_embedding_WINS_global_r2_while_LOSING_the_other_two():
    """This is what makes Arm C a hard case rather than a formality: the candidate beats the baseline
    on one metric. The verdict function must still come out NEGATIVE, driven by the de-confounded
    cell rather than by whichever metric the candidate happens to win. If this shape is ever absent,
    Arm C is not testing what it claims and the falsifier should say so."""
    fx = _fixture()["metrics"]
    g = fx["global_r2"]
    assert max(g["embedding"]) > g["structure_only"], "embedding must WIN global r2"
    s = fx["spearman"]
    assert max(s["embedding"]) < s["structure_only"], "embedding must LOSE spearman"
    w = fx["within_group_r2"]
    assert max(w["embedding"]) < w["structure_only"], "embedding must LOSE within-group r2"


def test_within_group_r2_is_negative_across_all_three_seeds():
    """The signature of the recorded failure: not merely worse than the baseline, but below zero --
    worse than predicting the mean."""
    w = _fixture()["metrics"]["within_group_r2"]["embedding"]
    assert all(v < 0 for v in w), w
    assert len(w) == 3


def test_the_minus_sign_trap_is_recorded_and_real():
    """A naive float() on the markdown's minus sign raises. Pinned so a future parser author does not
    rediscover it, and so the fixture's parsing_notes cannot quietly become false."""
    try:
        float(f"{MINUS}0.173")
        raised = False
    except ValueError:
        raised = True
    assert raised, "U+2212 now parses as a number; the parsing note is stale"
    assert _num(f"{MINUS}0.173") == -0.173
    assert MINUS in MD.read_text(encoding="utf-8")
    assert any("U+2212" in n for n in _fixture()["parsing_notes"])
