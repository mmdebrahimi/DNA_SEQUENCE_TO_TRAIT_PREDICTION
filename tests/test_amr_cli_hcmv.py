"""CLI routing tests for the HCMV cell via `dna-amr --drug <hcmv> --observed GENE:SUB[,...]`.

Wheel-only observed mode (no BLAST/Docker). Pins: R/S dispatch, the amr-mechanism-call-v1 record shape,
the organism relabel (Escherichia -> HCMV), genome-mode-is-v0.1 error, and --amrfinder-run rejection.
"""
from __future__ import annotations

import json

from dna_decode.amr.cli import main


def _run(argv, capsys):
    rc = main(argv)
    out = capsys.readouterr().out
    return rc, out


def test_cli_gcv_ul97_resistant_json(capsys):
    rc, out = _run(["--drug", "ganciclovir", "--observed", "UL97:M460V", "--json-only"], capsys)
    assert rc == 0
    d = json.loads(out)
    assert d["prediction"] == "R"
    assert d["drug"] == "ganciclovir"
    assert d["provenance"]["organism"] == "HCMV"          # bacterial default relabelled
    assert d["determinants"][0]["symbol"] == "UL97"


def test_cli_benign_susceptible(capsys):
    rc, out = _run(["--drug", "ganciclovir", "--observed", "UL54:C304S", "--json-only"], capsys)
    assert rc == 0
    assert json.loads(out)["prediction"] == "S"


def test_cli_letermovir_ul56(capsys):
    rc, out = _run(["--drug", "letermovir", "--observed", "UL56:C325W", "--json-only"], capsys)
    assert rc == 0
    assert json.loads(out)["prediction"] == "R"


def test_cli_multi_gene_observed(capsys):
    rc, out = _run(["--drug", "ganciclovir", "--observed", "UL97:M460V,UL54:F412L", "--json-only"], capsys)
    assert rc == 0
    d = json.loads(out)
    assert d["prediction"] == "R"
    assert {x["symbol"] for x in d["determinants"]} == {"UL97", "UL54"}


def test_cli_genome_mode_is_v01(capsys, tmp_path):
    g = tmp_path / "x.fna"
    g.write_text(">c\nACGT\n")
    rc, _ = _run(["--drug", "letermovir", "--genome-fasta", str(g)], capsys)
    assert rc == 3          # HCMV genome-FASTA mode deferred to v0.1


def test_cli_amrfinder_run_rejected(capsys, tmp_path):
    rc, _ = _run(["--drug", "ganciclovir", "--amrfinder-run", str(tmp_path)], capsys)
    assert rc == 2          # --amrfinder-run is bacterial-only
