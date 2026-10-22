"""Tests for the read-only promcalc toolchain probe.

All offline: every external fact the probe consults is an injectable parameter, so a "drive" is a
tmp_path directory and the verdict never depends on which volumes happen to be attached.

THE LOAD-BEARING TESTS are the two that keep a broken probe from reading as a clean negative --
`test_a_probe_that_resolves_nothing_is_INDETERMINATE_not_INPUTS_ABSENT` and
`test_probe_issues_no_subprocess`. This repo has a recorded instance of a scorer erroring and printing
a confident verdict from zero resolved calls, and the whole value of this probe is that reporting an
install is possible is a DIFFERENT claim from having performed one.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import promcalc_toolchain_probe as P  # noqa: E402


# ------------------------------------------------------------------ the ladder (pure)
def test_every_verdict_in_the_vocabulary_is_reachable():
    """A vocabulary entry no input can produce is decoration, not a verdict."""
    produced = {
        P.decide(probe_ok=False, inputs_found=True, env_resolves=True,
                 interpreter_available=True, free_gb_max=100.0),
        P.decide(probe_ok=True, inputs_found=False, env_resolves=True,
                 interpreter_available=True, free_gb_max=100.0),
        P.decide(probe_ok=True, inputs_found=True, env_resolves=True,
                 interpreter_available=True, free_gb_max=100.0),
        P.decide(probe_ok=True, inputs_found=True, env_resolves=False,
                 interpreter_available=False, free_gb_max=100.0),
        P.decide(probe_ok=True, inputs_found=True, env_resolves=False,
                 interpreter_available=True, free_gb_max=100.0),
    }
    assert produced == set(P.VERDICTS)


def test_a_probe_that_resolves_nothing_is_INDETERMINATE_not_INPUTS_ABSENT():
    """"I could not look" and "I looked and it is not there" must stay different findings."""
    assert P.decide(probe_ok=False, inputs_found=False, env_resolves=False,
                    interpreter_available=False, free_gb_max=0.0) == "INDETERMINATE_PROBE_FAILED"


def test_INPUTS_ABSENT_outranks_ENV_PRESENT():
    """A toolchain with nothing to run on does not make regeneration possible."""
    assert P.decide(probe_ok=True, inputs_found=False, env_resolves=True,
                    interpreter_available=True, free_gb_max=900.0) == "INPUTS_ABSENT"


def test_too_little_free_space_is_UNSATISFIABLE():
    assert P.decide(probe_ok=True, inputs_found=True, env_resolves=False,
                    interpreter_available=True,
                    free_gb_max=P.MIN_FREE_GB - 0.1) == "ENV_ABSENT_DEPS_UNSATISFIABLE"


def test_enough_free_space_is_SATISFIABLE():
    assert P.decide(probe_ok=True, inputs_found=True, env_resolves=False,
                    interpreter_available=True,
                    free_gb_max=P.MIN_FREE_GB + 0.1) == "ENV_ABSENT_DEPS_SATISFIABLE_HERE"


# ------------------------------------------------------------------ import extraction (pure)
def test_module_level_imports_are_column_zero_only():
    """Stage B's third-party import sits INSIDE a function; calling it module-level misdescribes it."""
    src = "import argparse\nfrom pathlib import Path\n\ndef f():\n    from promoter_calculator import X\n"
    assert P._module_level_imports(src) == ["argparse", "pathlib"]
    assert P._nested_imports(src) == ["promoter_calculator"]


def test_future_import_is_not_reported_as_a_dependency():
    assert P._module_level_imports("from __future__ import annotations\nimport json\n") == ["json"]


# ------------------------------------------------------------------ input discovery
def test_find_inputs_reports_a_reachable_substrate(tmp_path):
    d = tmp_path / "mpra"
    d.mkdir()
    (d / P.PEAK_TILE_FILE).write_bytes(b"not really gzip, and never opened")
    assert P.find_inputs(roots=[tmp_path]) == [str(d / P.PEAK_TILE_FILE).replace(chr(92), "/")]


def test_find_inputs_reports_nothing_when_the_substrate_is_absent(tmp_path):
    (tmp_path / "mpra").mkdir()
    assert P.find_inputs(roots=[tmp_path]) == []


def test_find_inputs_does_not_open_the_file(tmp_path, monkeypatch):
    """`genomewide.load_peak_tiles` owns GSE144621; a second reader is what the ownership guard forbids."""
    d = tmp_path / "mpra"
    d.mkdir()
    (d / P.PEAK_TILE_FILE).write_bytes(b"\x00\x01\x02")

    import gzip
    def _refuse(*a, **k):
        raise AssertionError("the probe must not open the substrate")
    monkeypatch.setattr(gzip, "open", _refuse)
    monkeypatch.setattr(Path, "read_bytes", lambda self: _refuse())
    assert len(P.find_inputs(roots=[tmp_path])) == 1


# ------------------------------------------------------------------ the record
def _rec(tmp_path, **kw):
    kw.setdefault("env_python", tmp_path / "nope" / "python.exe")
    kw.setdefault("roots", [tmp_path])
    kw.setdefault("vols", [{"drive": "Z:", "free_gb": 900.0, "total_gb": 931.0}])
    kw.setdefault("stage_b_source", "import json\n\ndef g():\n    from promoter_calculator import C\n")
    kw.setdefault("stage_a_source", "def export(cache_dir=None):\n    return 1\n")
    kw.setdefault("requirement_importable", False)
    return P.probe_toolchain(**kw)


def test_probe_issues_no_subprocess(tmp_path, monkeypatch):
    """Reporting that an install is possible is not the same claim as having performed one."""
    import subprocess

    def _refuse(*a, **k):
        raise AssertionError("the probe must not shell out")

    for name in ("run", "call", "check_call", "check_output", "Popen"):
        monkeypatch.setattr(subprocess, name, _refuse)
    rec = _rec(tmp_path)
    assert rec["installs_performed"] == 0


def test_probe_writes_nothing_without_out(tmp_path):
    before = sorted(p.name for p in tmp_path.iterdir())
    _rec(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == before


def test_record_carries_the_frozen_vocabulary_and_a_reachable_verdict(tmp_path):
    rec = _rec(tmp_path)
    assert rec["verdict_vocabulary"] == list(P.VERDICTS)
    assert rec["verdict"] in P.VERDICTS
    assert rec["schema"] == P.SCHEMA


def test_satisfiable_carries_the_isolation_CLAIM_and_its_STATUS(tmp_path):
    """The claim must TRAVEL with the verdict, whichever way it resolved.

    This test used to assert `isolation_required is True` and the words "SEPARATE environment",
    because stage B's docstring says `promoter_calculator` pulls sklearn 1.9.1. That was MEASURED
    FALSE on 2026-10-22: the distribution is `promotercalculator` 1.2.4, a py3-none-any pure-Python
    wheel with unpinned deps, and a dry-run resolution would install exactly one package and leave
    sklearn where it is.

    So the guard is inverted into the durable shape: do not pin the ANSWER (which is environment-
    dependent and already changed once), pin that the CLAIM, its STATUS and its EVIDENCE are all
    present. A verdict that silently dropped the claim would let a future reader re-assume isolation
    is required -- or re-assume it is not -- with nothing recorded either way.
    """
    d = tmp_path / "mpra"
    d.mkdir()
    (d / P.PEAK_TILE_FILE).write_bytes(b"x")
    rec = _rec(tmp_path)
    assert rec["verdict"] == "ENV_ABSENT_DEPS_SATISFIABLE_HERE"

    claim = rec["env"]["isolation_claim"]
    assert claim["claimed_by"]                      # who said it
    assert "sklearn" in claim["claim"]              # what was claimed
    assert "REFUTED" in claim["status"]             # how it resolved
    assert "dry-run" in claim["evidence"]           # on what evidence
    assert isinstance(rec["env"]["isolation_required"], bool)

    # The probe must NOT pretend it measured this itself -- it issues no subprocess by contract.
    assert "no subprocess" in claim["not_measured_by_this_probe"].lower()

    # The trap that 404s a naive install must stay recorded next to the claim.
    assert "promotercalculator" in claim["install_name_differs_from_import_name"]


def test_absent_substrate_yields_INPUTS_ABSENT(tmp_path):
    (tmp_path / "mpra").mkdir()
    assert _rec(tmp_path)["verdict"] == "INPUTS_ABSENT"


def test_present_env_yields_ENV_PRESENT(tmp_path):
    d = tmp_path / "mpra"
    d.mkdir()
    (d / P.PEAK_TILE_FILE).write_bytes(b"x")
    env = tmp_path / "python.exe"
    env.write_text("", encoding="utf-8")
    assert _rec(tmp_path, env_python=env)["verdict"] == "ENV_PRESENT"


def test_stage_b_is_reported_as_not_importing_dna_decode(tmp_path):
    """Stage B lives in an isolated interpreter that CANNOT import this package."""
    rec = _rec(tmp_path)
    assert rec["stage_b"]["imports_dna_decode_at_module_level"] is False
    assert "promoter_calculator" in rec["stage_b"]["nested_imports"]


def test_honest_limits_name_the_anchor_question(tmp_path):
    """A SATISFIABLE verdict must not read as "regeneration would reproduce the frozen bar"."""
    joined = " ".join(_rec(tmp_path)["honest_limits"])
    assert "anchors" in joined


def test_main_writes_only_when_given_out(tmp_path, capsys):
    out = tmp_path / "sub" / "probe.json"
    assert P.main(["--env-python", str(tmp_path / "nope.exe")]) == 0
    assert not out.exists()
    assert P.main(["--env-python", str(tmp_path / "nope.exe"), "--out", str(out)]) == 0
    assert out.is_file()
    import json
    assert json.loads(out.read_text(encoding="utf-8"))["schema"] == P.SCHEMA


def test_volumes_sees_at_least_one_real_drive():
    """The one deliberately non-hermetic case: a host with zero resolving drives cannot exist."""
    assert volumes_nonempty()


def volumes_nonempty() -> bool:
    return len(P.volumes()) >= 1
