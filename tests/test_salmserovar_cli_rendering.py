"""The rendered CLI surface must not contradict the call it is reporting.

Since the O decision procedure landed, the best-scoring O allele is routinely NOT the O call -- plain O9
has no allele of its own in either database and is reached from a `wbaV` hit. The evidence list is still
best-hit-per-axis, so an unmarked `O 9,46` row printed beside a `9:g,m:-` formula reads as a
contradiction, and the rule that produced the call was not rendered at all.

`A block carried only in JSON is not a disclosure` is this repo's own standard, and it applies to the
thing that DECIDED the call just as much as to a caveat. These tests drive the real `main()` with a
synthetic result, so they need no blastn, no DB and no genome.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar import cli as sv_cli  # noqa: E402


def _result(o="9", rule="A_wbaV_default_O9_tyr_not_discriminable_by_blast", best="9,46",
            serovar="Enteritidis") -> dict:
    """A branch-A result: the best O hit is `9,46` (a wbaV allele) but the CALL is plain `9`."""
    return {
        "status": "ok", "tool": "blastn", "method": "seqsero2_o_procedure_port_v1",
        "parameters": {"identity_threshold": 90.0, "coverage_threshold": 40.0},
        "o_antigen": o, "h1_antigen": "g,m", "h2_antigen": None, "o_antigen_rule": rule,
        "o_antigen_best_hit": best,
        "antigenic_formula": f"{o}:g,m:-", "serovar": serovar,
        "antigens": [
            {"axis": "H1", "antigen": "g,m", "best_allele": "H1__g,m__0",
             "percent_identity": 100.0, "percent_coverage": 100.0, "is_call": True},
            {"axis": "O", "antigen": best, "best_allele": "O__9,46__239__O-9,46_wbaV__1002",
             "percent_identity": 99.7, "percent_coverage": 100.0, "is_call": best == o},
        ],
    }


FIXTURE_FASTA = ROOT / "tests" / "fixtures" / "salmserovar" / "genome.fna"


def _run(monkeypatch, capsys, result: dict, extra: list[str] | None = None):
    """Drive the REAL `main()`. `call_serovar` is stubbed so no blastn/DB is needed, but everything
    downstream -- argument parsing, record assembly, rendering -- is the shipped code path."""
    if not FIXTURE_FASTA.exists():
        pytest.skip("salmserovar fixture genome absent")
    monkeypatch.setattr(sv_cli, "call_serovar", lambda *a, **k: result)
    argv = [str(FIXTURE_FASTA), "--sample-id", "s"]
    rc = sv_cli.main(argv + (extra or []))
    return rc, capsys.readouterr().out


def test_an_evidence_row_that_is_not_the_call_says_so(monkeypatch, capsys):
    """THE DEFECT. Printing `O 9,46` unmarked beside a `9:g,m:-` formula tells a reader the caller
    contradicted itself, when in fact the wbaV hit is the EVIDENCE FOR calling plain 9."""
    rc, out = _run(monkeypatch, capsys, _result())
    assert rc == 0
    o_line = next(ln for ln in out.splitlines() if ln.strip().startswith("O "))
    assert "NOT the call" in o_line, o_line


def test_the_rule_that_decided_the_o_call_is_printed(monkeypatch, capsys):
    """A differential-marker decision, a branch default and an ordinary best-allele pick are three
    different kinds of evidence. Reporting only the antigen string collapses them."""
    _, out = _run(monkeypatch, capsys, _result())
    assert "O call: 9" in out
    assert "A_wbaV_default_O9_tyr_not_discriminable_by_blast" in out


def test_no_contradiction_marker_when_the_best_hit_IS_the_call(monkeypatch, capsys):
    """The marker must not cry wolf on the ordinary case, or readers learn to ignore it."""
    _, out = _run(monkeypatch, capsys, _result(o="8", rule="C_argmax", best="8", serovar="Kentucky"))
    o_line = next(ln for ln in out.splitlines() if ln.strip().startswith("O "))
    assert "NOT the call" not in o_line, o_line
    assert "[rule: C_argmax]" in out


def test_the_json_record_carries_the_rule_and_the_best_hit(monkeypatch, capsys):
    _, out = _run(monkeypatch, capsys, _result(), extra=["--json-only"])
    rec = json.loads(out)
    assert rec["o_antigen"] == "9"
    assert rec["o_antigen_rule"] == "A_wbaV_default_O9_tyr_not_discriminable_by_blast"
    assert rec["o_antigen_best_hit"] == "9,46"
    o_row = next(r for r in rec["antigen_detail"] if r["axis"] == "O")
    assert o_row["is_call"] is False


def test_an_unresolved_o_call_still_reports_its_rule(monkeypatch, capsys):
    """`C_guard_1319_marker_only` is an ABSTENTION with a reason -- the reason is the useful part, and
    dropping it would make a principled refusal look like a blank."""
    res = _result(o=None, rule="C_guard_1319_marker_only", best="1,3,19", serovar=None)
    res["antigenic_formula"] = "O?:g,m:-"
    _, out = _run(monkeypatch, capsys, res)
    assert "O call: unresolved" in out
    assert "C_guard_1319_marker_only" in out


def test_a_legacy_db_run_is_visible_in_the_output(monkeypatch, capsys):
    """The DB is gitignored, so another machine may hold a role-less build. That run does NOT execute
    the ported procedure, and the output has to show it rather than look like a normal call."""
    res = _result(o="9,46", rule="legacy_db_no_role_field", best="9,46", serovar="Hillingdon")
    _, out = _run(monkeypatch, capsys, res)
    assert "legacy_db_no_role_field" in out


def test_rendering_survives_a_result_without_the_new_fields(monkeypatch, capsys):
    """Backward compatibility: an older result dict has no `o_antigen_rule` and no `is_call`. It must
    render as it always did rather than raising or emitting a bare `[rule: None]`."""
    res = _result()
    del res["o_antigen_rule"]
    for row in res["antigens"]:
        row.pop("is_call", None)
    rc, out = _run(monkeypatch, capsys, res)
    assert rc == 0
    assert "O call:" not in out
    assert "NOT the call" not in out
