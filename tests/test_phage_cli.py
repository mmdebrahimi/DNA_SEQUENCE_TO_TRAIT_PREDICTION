"""Tests for `dna-decode phage` / `dna-phage` (dna_decode/phage/cli.py).

Offline: the wheel-only --lineage catalogue path (CALLED / abstain), JSON output, and the
genome-fasta degrade path (reference set absent -> INDETERMINATE, never a fabricated call).
"""
from __future__ import annotations

import json

import dna_decode.phage.cli as pcli
from dna_decode.cli import main as unified_main


def test_lineage_clade_conserved_called(capsys):
    rc = pcli.main(["--lineage", "Vequintavirus", "--json"])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["status"] == "CALLED" and out["receptor"] == "ECA"
    assert out["clade_conserved"] is True


def test_lineage_rbp_variable_abstains(capsys):
    # T-even is RBP-variable -> INDETERMINATE (no fabricated single receptor), still exit 0
    rc = pcli.main(["--lineage", "Tequatrovirus,Straboviridae", "--json"])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["status"] == "INDETERMINATE" and out["receptor"] is None
    assert "varies" in out["reason"].lower()


def test_lineage_uncatalogued_indeterminate(capsys):
    rc = pcli.main(["--lineage", "Novelgenus", "--json"])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["status"] == "INDETERMINATE" and out["receptor"] is None


def test_genome_fasta_missing_reference_degrades(tmp_path, capsys):
    q = tmp_path / "q.fna"
    q.write_text(">q\nACGT\n", encoding="utf-8")
    # point at a non-existent reference set -> actionable INDETERMINATE, not a crash / fabricated call
    rc = pcli.main(["--genome-fasta", str(q), "--reference-manifest", str(tmp_path / "none.tsv"),
                    "--reference-dir", str(tmp_path / "none"), "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["status"] == "INDETERMINATE"
    assert "fetch_basel_genomes" in out["reason"]


def test_routes_through_unified_cli(capsys):
    # `dna-decode phage --lineage ...` must dispatch to the phage CLI
    rc = unified_main(["phage", "--lineage", "Enquatrovirus"])
    assert rc == 0
    assert "NfrA" in capsys.readouterr().out


def test_scope_rail_present(capsys):
    pcli.main(["--lineage", "Vequintavirus"])
    assert "receptor-CLASS only" in capsys.readouterr().out


def test_rbp_fasta_path_calls_mixed_clade_receptor(tmp_path, capsys):
    # pull one reference RBP protein into a query fasta -> --rbp-fasta calls its receptor (offline, wheel-only)
    from dna_decode.phage.rbp_caller import DEFAULT_RBP_REFERENCE
    label = receptor = None; seq = []
    for line in open(DEFAULT_RBP_REFERENCE, encoding="utf-8"):
        if line.startswith(";"): continue
        if line.startswith(">"):
            if label: break
            label, _, receptor = line[1:].strip().partition("|")
        elif label: seq.append(line.strip())
    q = tmp_path / "rbp.faa"; q.write_text(f">q\n{''.join(seq)}\n", encoding="utf-8")
    rc = pcli.main(["--rbp-fasta", str(q), "--json"])
    import json
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["status"] == "CALLED" and out["receptor"] == receptor
