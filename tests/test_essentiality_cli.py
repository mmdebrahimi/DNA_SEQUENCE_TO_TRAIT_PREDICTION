"""CLI tests for `dna-decode essentiality` (single-gene + JSON + error paths). Offline, no data."""
import json
import dna_decode.essentiality.cli as ec
from dna_decode.cli import main as unified_main


def test_single_gene_essential(capsys):
    rc = ec.main(["--gene", "ftsZ", "--product", "cell division protein FtsZ"])
    out = capsys.readouterr().out
    assert rc == 0 and "ESSENTIAL" in out


def test_single_gene_nonessential_json(capsys):
    rc = ec.main(["--gene", "lacZ", "--product", "beta-galactosidase lactose catabolism", "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["prediction"] == "non_essential" and d["tier"] == "KNOWLEDGE_BASELINE"


def test_needs_input(capsys):
    rc = ec.main([])
    assert rc == 2 and "feature-table" in capsys.readouterr().err


def test_routes_through_unified_cli(capsys):
    rc = unified_main(["essentiality", "--gene", "rpsA", "--product", "30S ribosomal protein S1", "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["prediction"] == "essential"
