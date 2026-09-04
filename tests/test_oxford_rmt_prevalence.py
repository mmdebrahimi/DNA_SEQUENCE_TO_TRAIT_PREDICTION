"""The probe's zero must be absence, not a parse failure -- and the guard that says so must fire.

A negative result from a scan that silently read nothing is indistinguishable, in the artifact, from a
negative result from a real absence. These tests pin BOTH directions: the guard refuses on a dead scan,
and it does not refuse on a live one.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "oxford_rmt_prevalence_probe.py"
ARTIFACT = ROOT / "wiki" / "oxford_rmt_prevalence_2026-09-03.json"

HEADER = ("Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\t"
          "Scope\tElement type\tElement subtype\tClass\tSubclass\tMethod\n")


def _row(name: str, sym: str, cls: str = "AMINOGLYCOSIDE", sub: str = "GENTAMICIN") -> str:
    return f"{name}\tp\tc\t1\t2\t+\t{sym}\tseq\tcore\tAMR\tAMR\t{cls}\t{sub}\tEXACTX\n"


def _cohort(tmp: Path, symbols: list[tuple[str, str]]) -> Path:
    d = tmp / "oxford"
    d.mkdir()
    (d / "amrfinder.tsv").write_text(HEADER + "".join(_row(n, s) for n, s in symbols), encoding="utf-8")
    (d / "main_data.csv").write_text(
        "guuid,Gentamicin_lower,Gentamicin_upper\n"
        + "".join(f"{n},4,5\n" for n, _ in symbols), encoding="utf-8")
    return d


def _run(cohort: Path, out: Path):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--oxford-dir", str(cohort), "--out", str(out)],
        capture_output=True, text=True, cwd=str(ROOT))


def test_guard_refuses_when_the_scan_found_almost_nothing(tmp_path):
    """A near-empty scan must REFUSE (exit 3), not report a clean zero."""
    cohort = _cohort(tmp_path, [("iso1", "aac(3)-IId"), ("iso2", "aac(3)-IId")])
    out = tmp_path / "o.json"
    r = _run(cohort, out)
    assert r.returncode == 3, r.stdout + r.stderr
    assert "REFUSING to report" in r.stderr
    assert not out.exists(), "a refused run must not leave an artifact behind"


def test_guard_passes_on_a_live_scan_and_finds_a_planted_carrier(tmp_path):
    """The same guard must NOT block a real scan -- otherwise it could never report anything."""
    syms = [(f"iso{i}", s) for i, s in enumerate(
        ["aac(3)-IId", "aac(3)-IIe", "aadA1", "aadA5", "aph(3')-Ia", "aadA16", "rmtB"])]
    cohort = _cohort(tmp_path, syms)
    out = tmp_path / "o.json"
    r = _run(cohort, out)
    assert r.returncode == 0, r.stdout + r.stderr
    d = json.loads(out.read_text(encoding="utf-8"))
    assert d["n_rmt_npma_carriers"] == 1
    assert d["verdict"] == "CARRIERS_PRESENT_TESTABLE"
    assert d["non_vacuity"]["symbols_detected_once"], "singleton proof must be populated"


def test_arma_is_not_counted_as_a_rescue_case(tmp_path):
    """armA is already counted by the frozen rule; conflating it with rmt would overstate the gap."""
    syms = [(f"iso{i}", s) for i, s in enumerate(
        ["aac(3)-IId", "aac(3)-IIe", "aadA1", "aadA5", "aph(3')-Ia", "aadA16", "armA"])]
    out = tmp_path / "o.json"
    r = _run(_cohort(tmp_path, syms), out)
    assert r.returncode == 0, r.stdout + r.stderr
    d = json.loads(out.read_text(encoding="utf-8"))
    assert d["n_rmt_npma_carriers"] == 0 and d["n_arma_carriers"] == 1


@pytest.mark.skipif(not ARTIFACT.exists(), reason="committed Oxford probe artifact absent")
def test_committed_artifact_records_a_proven_absence_not_a_silent_one():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n_rmt_npma_carriers"] == 0
    assert d["verdict"] == "NO_CARRIERS_CANNOT_TEST"
    # The whole point: the zero is only meaningful because the scan demonstrably worked.
    assert d["non_vacuity"]["distinct_aminoglycoside_symbols"] >= 20
    assert len(d["non_vacuity"]["symbols_detected_once"]) >= 3
    # And it must not be mistaken for evidence about the rule.
    assert any("does NOT test" in s for s in d["what_this_does_not_show"])
