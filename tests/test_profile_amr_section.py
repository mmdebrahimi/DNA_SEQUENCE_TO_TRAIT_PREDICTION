"""Pins the unified-genome-profile AMR section (productization move 2): the profile now folds the AMR R/S
decoder + the move-1 inline trust badges into one report, offline-safely. Uses a committed main.tsv fixture
(no Docker, no gitignored data/ dependency) so it runs in CI.
"""
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.profile import cli as profile_cli  # noqa: E402

_MAIN_TSV = REPO / "tests" / "fixtures" / "amr_mini" / "main.tsv"
_FASTA = REPO / "tests" / "fixtures" / "ecoli_mini" / "genome.fna"
_DRUGS = ["ciprofloxacin", "ceftriaxone", "tetracycline", "gentamicin"]


def test_amr_section_calls_with_trust_badges():
    r = profile_cli._amr(_MAIN_TSV, _DRUGS, "Escherichia", None)
    assert r["status"] == "ok" and r["organism"] == "Escherichia"
    # the fixture is built so cipro (2 QRDR points) + cef (CTX-M) are R; tet/gent S
    assert r["calls"]["ciprofloxacin"]["prediction"] == "R"
    assert r["calls"]["ceftriaxone"]["prediction"] == "R"
    assert r["calls"]["tetracycline"]["prediction"] == "S"
    # every call carries an inline validation badge (the move-1 trust surface, composed here)
    for drug in _DRUGS:
        v = r["calls"][drug]["validation"]
        assert v["tier"] == "INDEPENDENT_MEASURED" and v["independent"] is True
        assert v["source_card"].endswith("amr_portal_independent_report_card.md")


def test_amr_section_offline_safe():
    assert profile_cli._amr(None, _DRUGS, "Escherichia", None)["status"] == "unavailable"
    assert profile_cli._amr("does/not/exist/main.tsv", _DRUGS, "Escherichia", None)["status"] == "unavailable"


def test_amr_section_per_drug_error_does_not_sink_section():
    # a bogus drug must not crash the whole section; the section stays ok, others still resolve
    r = profile_cli._amr(_MAIN_TSV, ["ciprofloxacin", "not_a_real_drug"], "Escherichia", None)
    assert r["status"] == "ok"
    assert "ciprofloxacin" in r["calls"]


def test_profile_main_includes_amr_section(capsys):
    rc = profile_cli.main([str(_FASTA), "--amrfinder-run", str(_MAIN_TSV.parent),
                           "--amr-organism", "Escherichia", "--json-only"])
    assert rc == 0
    out = capsys.readouterr().out
    rec = json.loads(out)
    assert "amr" in rec["decoders"]
    amr = rec["decoders"]["amr"]
    assert amr["status"] == "ok"
    assert amr["calls"]["ciprofloxacin"]["validation"]["tier"] == "INDEPENDENT_MEASURED"
    # the AMR section counts toward the unified decoder tally
    assert rec["decoders_total"] == 6


def test_profile_main_amr_unavailable_without_source(capsys):
    rc = profile_cli.main([str(_FASTA), "--json-only"])   # no --amrfinder-run / --run-amrfinder
    assert rc == 0
    rec = json.loads(capsys.readouterr().out)
    assert rec["decoders"]["amr"]["status"] == "unavailable"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
