"""Pins the CYP2C19 validation harness (scripts/pgx_cyp2c19_validate.py).

Unit tests on the filename->expected-diplotype parser (incl. the `rs...` absorb regression) always run.
The real-fixture concordance test runs against the COMMITTED PharmCAT fixtures under tests/data/ -> it
asserts the 6/6 core number the report card publishes is reproducible, and that the two non-core
blind-spot cases behave exactly as documented (s1s35 -> *1/*1; s1s4b aliases to *1/*17 via the shared
rs12248560). FAITHFUL-TO-PHARMCAT tier (the reference tool's own fixtures), NOT the GeT-RM independent panel.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from pgx_cyp2c19_validate import (  # noqa: E402
    _norm_diplo,
    build_report,
    expected_from_filename,
    validate_dir,
)

_FIX = REPO / "tests" / "data" / "pgx_cyp2c19"


@pytest.mark.parametrize("stem,expected", [
    ("s1s1", ("*1", "*1")),
    ("s1s2", ("*1", "*2")),
    ("s2s3", ("*2", "*3")),
    ("s1s4b", ("*1", "*4b")),                       # single-letter allele suffix kept
    ("s1s1rs12248560missing", ("*1", "*1")),        # the rs-absorb regression: NOT *1/*1r
    ("s1s2rs58973490het", ("*1", "*2")),
    ("s1s2rs3758581missing", ("*1", "*2")),
    ("s15s28", ("*15", "*28")),
])
def test_expected_from_filename(stem, expected):
    assert expected_from_filename(stem) == expected


@pytest.mark.parametrize("stem", ["noCall", "rs12769205only", "sUnks17"])
def test_non_diplotype_fixtures_skipped(stem):
    assert expected_from_filename(stem) is None


def test_norm_diplotype_order_independent():
    assert _norm_diplo("*2", "*1") == "*1/*2"
    assert _norm_diplo("*17", "*1") == "*1/*17"
    assert _norm_diplo("*3", "*2") == "*2/*3"


@pytest.mark.skipif(not _FIX.exists() or not any(_FIX.glob("*.vcf")),
                    reason="PharmCAT fixtures not present")
def test_pharmcat_core_concordance_is_6_of_6():
    rows = validate_dir(_FIX)
    rep = build_report(rows, "pharmcat")
    assert rep["core_diplotype_hits"] == "6/6"
    assert rep["core_phenotype_hits"] == "6/6"
    assert rep["caller_is_independent_baseline"] is False  # faithful-to-tool, not independent


@pytest.mark.skipif(not _FIX.exists() or not any(_FIX.glob("*.vcf")),
                    reason="PharmCAT fixtures not present")
def test_documented_blindspots():
    rows = {r["fixture"]: r for r in validate_dir(_FIX)}
    # *35 is rs12769205-defined (not in the core SNP set) -> mis-called *1/*1
    assert rows["s1s35"]["predicted_diplotype"] == "*1/*1"
    assert rows["s1s35"]["is_core"] is False
    # *4b carries the *17 SNP (rs12248560) + rs28399504; single-SNP *17 proxy aliases it to *1/*17
    assert rows["s1s4b"]["predicted_diplotype"] == "*1/*17"
    assert rows["s1s4b"]["is_core"] is False


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
