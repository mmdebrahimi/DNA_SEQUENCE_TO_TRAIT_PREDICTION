"""The SeqSero2 comparator must fail LOUDLY on schema drift, never quietly as an abstention.

WHY THIS MATTERS MORE THAN AN ORDINARY PARSE CHECK. This parser produces the REFERENCE column the
serovar cell measures its own delta against. Every field is read `.get(...) or ""`, and `score()`
maps an empty prediction to `no_call` -- so a renamed header would yield five empty strings on a
perfectly non-empty row, and SeqSero2 would appear to ABSTAIN on every isolate. Our delta would
improve, for exactly the wrong reason. The pre-existing checks (container crash, missing file, zero
rows) cannot see this: the run succeeds, the file exists, the row is there.

These tests exercise the parse path directly with a stubbed container, so they need no Docker.
"""
from __future__ import annotations

import importlib.util as _u
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_spec = _u.spec_from_file_location("ss2v", ROOT / "scripts" / "salmserovar_seqsero2_validate.py")
ss2v = _u.module_from_spec(_spec)
_spec.loader.exec_module(ss2v)

GOOD_HEADER = "\t".join(ss2v._REQUIRED_COLUMNS)
GOOD_ROW = "Typhimurium\tI 4,[5],12:i:1,2\t4\ti\t1,2"


def _run_with_tsv(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tsv_text: str) -> dict:
    """Drive `run_seqsero2` with a stubbed docker call that writes `tsv_text` as the result."""
    stage = tmp_path / "stage"
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")

    class _P:
        returncode = 0

    def fake_run(cmd, **kw):
        out = stage / "out"
        out.mkdir(parents=True, exist_ok=True)
        (out / "SeqSero_result.tsv").write_text(tsv_text, encoding="utf-8")
        return _P()

    monkeypatch.setattr(ss2v.subprocess, "run", fake_run)
    return ss2v.run_seqsero2(fasta, stage, "img")


def test_a_well_formed_result_parses(tmp_path, monkeypatch):
    """Non-vacuity: if this failed, every test below would 'pass' for the wrong reason."""
    got = _run_with_tsv(tmp_path, monkeypatch, f"{GOOD_HEADER}\n{GOOD_ROW}\n")
    assert "error" not in got
    assert got == {"serotype": "Typhimurium", "profile": "I 4,[5],12:i:1,2",
                   "O": "4", "H1": "i", "H2": "1,2"}


@pytest.mark.parametrize("renamed", list(ss2v._REQUIRED_COLUMNS))
def test_any_renamed_column_is_reported_as_schema_drift(tmp_path, monkeypatch, renamed):
    """Each required header, dropped one at a time -- a gate that only catches SOME renames would
    still let the comparator deflate silently on the others."""
    header = "\t".join(c if c != renamed else f"{c} (v2)" for c in ss2v._REQUIRED_COLUMNS)
    got = _run_with_tsv(tmp_path, monkeypatch, f"{header}\n{GOOD_ROW}\n")
    assert "error" in got, f"a renamed {renamed!r} parsed silently instead of erroring"
    assert got["error"].startswith("schema_drift"), got["error"]
    assert renamed in got["error"], "the error must name the column that moved"


def test_the_old_behaviour_would_have_been_indistinguishable_from_abstention(tmp_path, monkeypatch):
    """Pins the FAILURE MODE, not just the fix: with the headers renamed, the values the old parser
    would have produced score as `no_call` for every field -- i.e. SeqSero2 'abstaining'."""
    idx: dict = {}
    assert ss2v.score("", "Typhimurium", idx) == "no_call"
    assert ss2v.score(None, "Typhimurium", idx) == "no_call"
    header = "\t".join(f"{c} (v2)" for c in ss2v._REQUIRED_COLUMNS)
    got = _run_with_tsv(tmp_path, monkeypatch, f"{header}\n{GOOD_ROW}\n")
    assert got.get("error", "").startswith("schema_drift")


def test_right_headers_but_an_entirely_empty_row_is_its_own_error(tmp_path, monkeypatch):
    """Distinct from schema drift on purpose: the headers are fine, so this is not a rename -- but
    SeqSero2 fills at least one field on any successful run, so an all-empty parse is still not an
    abstention and must not be recorded as one."""
    got = _run_with_tsv(tmp_path, monkeypatch, f"{GOOD_HEADER}\n\t\t\t\t\n")
    assert got.get("error") == "all_fields_empty"


def test_a_partially_filled_row_is_accepted(tmp_path, monkeypatch):
    """A genuine SeqSero2 abstention on some axes must still parse -- over-tightening the gate into
    'every field must be populated' would reject real results."""
    got = _run_with_tsv(tmp_path, monkeypatch, f"{GOOD_HEADER}\nTyphimurium\t\t4\t\t\n")
    assert "error" not in got
    assert got["serotype"] == "Typhimurium" and got["H1"] == ""


def test_the_pinned_columns_are_the_ones_the_parser_actually_reads():
    """The constant and the parse body must not drift apart -- a pin naming columns nobody reads
    would be ceremony that passes while the real reads stay unguarded."""
    src = (ROOT / "scripts" / "salmserovar_seqsero2_validate.py").read_text(encoding="utf-8")
    body = src.split("def run_seqsero2", 1)[1]
    for col in ss2v._REQUIRED_COLUMNS:
        assert f'r.get("{col}")' in body, f"{col!r} is pinned but never read"
