"""CLI tests for `dna-decode metabolic` (single-source + JSON + list + error paths). Offline, no data."""
import json
import dna_decode.metabolic.cli as mc
from dna_decode.cli import main as unified_main


def test_lactose_positive_text(capsys):
    rc = mc.main(["--source", "lactose", "--genes", "lacZ,lacY"])
    out = capsys.readouterr().out
    assert rc == 0 and "UTILIZES" in out


def test_citrate_aerobic_anchor_json(capsys):
    rc = mc.main(["--source", "citrate", "--genes", "citD,citE,citF,citT",
                  "--condition", "aerobic", "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["capability"] == "cannot_utilize" and d["tier"] == "KNOWLEDGE_BASELINE"
    assert d["transporter_present"] is True and d["transporter_expressed_under_condition"] is False


def test_list(capsys):
    rc = mc.main(["--list"])
    out = capsys.readouterr().out
    assert rc == 0 and "citrate" in out and "lactose" in out


def test_unknown_source_errors(capsys):
    rc = mc.main(["--source", "plutonium", "--genes", "x"])
    assert rc == 2 and "unknown substrate" in capsys.readouterr().err


def test_needs_input(capsys):
    rc = mc.main([])
    assert rc == 2 and "--list" in capsys.readouterr().err


def test_routes_through_unified_cli(capsys):
    rc = unified_main(["metabolic", "--source", "citrate", "--genes", "citD,citE,citF,citT",
                       "--condition", "anaerobic", "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["capability"] == "utilizes"
