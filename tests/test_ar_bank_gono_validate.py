"""Tests for the AR Bank gono scorer wiring: main.tsv symbol parse + the rule_predictor -> call_ng_amr chain.

The heavy CLI main() (network genome-resolve + AMRFinder) is exercised by the real scoring run; here we
pin the NEW reusable pieces: parse_determinant_symbols (synthetic + the real cached AR#0165 main.tsv) and
rule_predictor (with _run_dir monkeypatched so no Docker/network is touched).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scripts.external_cohort_revalidate as ecr  # noqa: E402
from dna_decode.organism_rules.neisseria_amr import call_ng_amr  # noqa: E402

_MAIN_TSV_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\t"
                    "Scope\tType\tSubtype\tClass\tSubclass\n")
_REAL_AR0165 = Path("data/raw/ar_bank_gono_smoke/amrfinder_runs/GCA_042036815.1/main.tsv")


def _write_main_tsv(path: Path, symbols):
    rows = [_MAIN_TSV_HEADER]
    for s in symbols:
        rows.append(f"NA\tctg1\t1\t2\t+\t{s}\tname\tcore\tAMR\tPOINT\tCLASS\tSUBCLASS\n")
    path.write_text("".join(rows), encoding="utf-8")


def test_parse_determinant_symbols_synthetic(tmp_path):
    mt = tmp_path / "main.tsv"
    _write_main_tsv(mt, ["gyrA_S91F", "penA_G545S", "blaTEM-1", ""])   # blank skipped
    syms = ecr.parse_determinant_symbols(mt)
    assert syms == ["gyrA_S91F", "penA_G545S", "blaTEM-1"]
    assert ecr.parse_determinant_symbols(tmp_path / "absent.tsv") == []   # missing file -> []


def test_rule_predictor_applies_call_ng_amr(tmp_path, monkeypatch):
    mt_dir = tmp_path / "run"
    mt_dir.mkdir()
    _write_main_tsv(mt_dir / "main.tsv", ["gyrA_S91F", "penA_G545S"])
    # monkeypatch the run-dir resolver so no Docker/network fires; ensure_run must NOT be called
    monkeypatch.setattr("scripts.organism_drug_validate._run_dir", lambda gca, own, glob: mt_dir)
    def _boom(*a, **k):
        raise AssertionError("ensure_run should not run when _run_dir resolves")
    monkeypatch.setattr("scripts.organism_drug_validate.ensure_run", _boom)

    predict_cro = ecr.rule_predictor("ceftriaxone", tmp_path / "o", tmp_path / "g", "glob",
                                     call_ng_amr, "Neisseria_gonorrhoeae")
    assert predict_cro("GCA_x") == "R"        # penA mosaic -> cef R
    predict_cip = ecr.rule_predictor("ciprofloxacin", tmp_path / "o", tmp_path / "g", "glob",
                                     call_ng_amr, "Neisseria_gonorrhoeae")
    assert predict_cip("GCA_x") == "R"        # gyrA S91F -> cipro R
    predict_azm = ecr.rule_predictor("azithromycin", tmp_path / "o", tmp_path / "g", "glob",
                                     call_ng_amr, "Neisseria_gonorrhoeae")
    assert predict_azm("GCA_x") == "S"        # no 23S -> azithro S


def test_rule_predictor_missing_run_is_indeterminate(tmp_path, monkeypatch):
    monkeypatch.setattr("scripts.organism_drug_validate._run_dir", lambda gca, own, glob: None)
    monkeypatch.setattr("scripts.organism_drug_validate.ensure_run", lambda *a, **k: None)
    predict = ecr.rule_predictor("ceftriaxone", tmp_path / "o", tmp_path / "g", "glob",
                                 call_ng_amr, "Neisseria_gonorrhoeae")
    assert predict("GCA_x") == "INDETERMINATE"


def test_real_cached_ar0165_main_tsv():
    """R3 real-surface: the actual AMRFinder -O Neisseria_gonorrhoeae output for AR#0165 parses + scores.
    Skips when the gitignored cached run is absent (fresh checkout / CI)."""
    if not _REAL_AR0165.exists():
        import pytest
        pytest.skip("cached AR#0165 gono main.tsv absent (data/raw is gitignored)")
    syms = ecr.parse_determinant_symbols(_REAL_AR0165)
    assert "gyrA_S91F" in syms and any(s.startswith("penA_") for s in syms)
    assert call_ng_amr("ciprofloxacin", syms)["prediction"] == "R"
    assert call_ng_amr("ceftriaxone", syms)["prediction"] == "R"
    assert call_ng_amr("azithromycin", syms)["prediction"] == "S"
    assert call_ng_amr("gentamicin", syms)["prediction"] == "INDETERMINATE"


def test_scorer_scope_excludes_gentamicin():
    import scripts.ar_bank_gono_validate as g
    assert "gentamicin" not in g.SCORABLE_DRUGS
    assert set(g.SCORABLE_DRUGS) == {"azithromycin", "cefixime", "ceftriaxone", "ciprofloxacin",
                                     "penicillin", "tetracycline"}
