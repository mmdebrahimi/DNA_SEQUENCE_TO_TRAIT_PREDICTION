"""Tests for the `dna-identify` CLI surface + its registration (Step 5b).

Offline: no Docker, no sketch. The runner is stubbed; what is tested here is the ARGPARSE surface, the
derived defaults, the coverage-honest listing, and the registration contract.

The repo-wide guard `tests/test_advertised_commands.py` already covers this command's advertised
examples automatically -- its `_console_scripts()` helper DISCOVERS routes from
`pyproject.toml [project.scripts]`, so adding the console entry enrolled `dna-identify` with no edit
there. That was verified, not assumed, before writing this file.
"""
from __future__ import annotations

import json

import pytest

from dna_decode.identify import cli as icli
from dna_decode.identify.thresholds import AMBIGUITY_MARGIN, MAX_DISTANCE


# --------------------------------------------------------------------------- argparse surface

def test_the_advertised_examples_all_parse():
    """Every command this module's docstring tells a user to run must resolve through the REAL parser.
    An advertised command is a promise; a rendered string is not a runnable command."""
    p = icli.build_parser()
    p.parse_args(["--genome-fasta", "g.fna"])
    p.parse_args(["--genome-fasta", "g.fna", "--json"])
    p.parse_args(["--list-organisms"])


def test_genome_fasta_is_required_unless_listing(capsys):
    """argparse.error exits 2; the point is it REFUSES rather than running on nothing."""
    with pytest.raises(SystemExit) as e:
        icli.main([])
    assert e.value.code == 2


def test_argparse_defaults_DERIVE_from_the_frozen_constants():
    """THE SHIPPED-DEFAULT TRAP, which has already cost this repo real accuracy: a CLI that restates a
    validated constant keeps overriding it after the constant moves (a `--coverage default=80.0` ran a
    threshold replaced five days earlier, worth -0.225). Equality today is necessary but not
    sufficient -- the literal must not be hardcoded either."""
    import inspect
    a = icli.build_parser().parse_args(["--genome-fasta", "g.fna"])
    assert a.max_distance == MAX_DISTANCE
    assert a.ambiguity_margin == AMBIGUITY_MARGIN

    src = inspect.getsource(icli.build_parser)
    assert "default=MAX_DISTANCE" in src, "the default must reference the constant, not a literal"
    assert "default=AMBIGUITY_MARGIN" in src
    assert f"default={MAX_DISTANCE}" not in src, "hardcoded literal: equality today drifts tomorrow"


# --------------------------------------------------------------------------- coverage honesty

def test_list_organisms_separates_COVERED_from_merely_DECLARED(tmp_path, capsys, monkeypatch):
    """A REAL over-claim found by running the CLI, not by reading it.

    The vocabulary declares 14 organisms; the reference sketch covers 10. A genome of the other 4 can
    only ABSTAIN, so listing the vocabulary alone advertises coverage that does not exist -- the same
    class as `dna-hla` advertising two SNP tags its own measurement had demoted.
    """
    man = tmp_path / "identify_reference_manifest_2026-01-01.json"
    man.write_text(json.dumps({"accessions_by_organism": {
        "escherichia_coli": ["GCA_1.1"], "campylobacter": ["GCA_2.1"]}}), encoding="utf-8")

    assert icli.main(["--list-organisms", "--manifest", str(man)]) == 0
    out = capsys.readouterr().out
    assert "2 of" in out and "have reference coverage" in out
    assert "IDENTIFIABLE" in out and "escherichia_coli" in out
    assert "DECLARED but NOT in the reference" in out
    assert "ABSTAINS, it is NOT identified" in out
    # an organism with no reference genomes must appear in the DECLARED section, not the identifiable one
    ident, _, declared = out.partition("DECLARED but NOT in the reference")
    assert "candida_auris" in declared and "candida_auris" not in ident


def test_list_organisms_survives_a_missing_manifest(tmp_path, capsys):
    """Degrades to 'nothing is covered' rather than crashing -- and must NOT then claim coverage."""
    assert icli.main(["--list-organisms", "--manifest", str(tmp_path / "nope.json")]) == 0
    out = capsys.readouterr().out
    assert "0 of" in out and "DECLARED but NOT in the reference" in out


# --------------------------------------------------------------------------- fault vs abstention

def test_an_unavailable_reference_exits_3_and_emits_NO_abstention_record(tmp_path, capsys, monkeypatch):
    """An infrastructure fault must never read as 'organism unsupported'. Exit 3, stderr, and NOT a
    record with abstained=true -- otherwise a wedged Docker mount looks like a biological answer."""
    from dna_decode.identify.runner import ReferenceUnavailable

    def boom(*a, **k):
        raise ReferenceUnavailable("sketch missing")
    monkeypatch.setattr(icli, "query", boom)
    monkeypatch.setattr(icli, "_manifest_labels", lambda m: {"GCA_1.1": "escherichia_coli"})

    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    rc = icli.main(["--genome-fasta", str(g), "--json"])
    cap = capsys.readouterr()
    assert rc == 3
    assert "REFERENCE UNAVAILABLE" in cap.err
    assert "abstained" not in cap.out


# --------------------------------------------------------------------------- registration contract

def test_the_route_declares_a_regime():
    """`regime_for_route` RAISES on an unknown route, so a decoder cannot ship into a regime silently."""
    from dna_decode.data.cell_regime import regime_for_route
    assert regime_for_route("dna-identify") == "curated_catalog_exists"


def test_the_trait_is_reachable_and_contracted():
    from dna_decode.cli import TRAITS
    from dna_decode.data.cell_registry import cells
    assert "identify" in TRAITS
    got = [c for c in cells() if c.route == "dna-identify"]
    assert len(got) == 1, got
    c = got[0]
    # The registry's fields are STRUCTURED, which I had to learn by failing the existing guards:
    # `native_abstention` is a TERM string, `abstention_vocab` a single enum member, `falsifier_ref` a
    # repo-relative PATH that must exist, and `incoming_data_gate` a comma-separated G1..G10 list.
    # Prose in any of them is a defect, not richer documentation.
    from dna_decode.data.cell_registry import AbstentionVocab
    assert c.native_abstention == "ABSTAIN"
    assert c.abstention_vocab is AbstentionVocab.ABSTAIN_BY_DESIGN
    assert c.falsifier_ref == "scripts/identify_validate.py"
    assert c.incoming_data_gate == "G1", "the circular-label gate is the one that bites this cell"

    # the ABSTAIN REASON vocabulary lives in the decision core, not in the contract field
    from dna_decode.identify.core import AbstainReason
    assert {r.value for r in AbstainReason} >= {
        "no_hit", "above_max_distance", "ambiguous_top2", "reference_unavailable",
        "unlabelled_reference"}

    # the two measured limits must travel with the cell
    assert "69.2%" in c.validation_slice and "congener" in c.validation_slice.lower()
    assert "lineage" in c.demotion_rule.lower()
    assert "lineage-disjoint" in c.demotion_rule.lower(), "the open falsifier must be named"


def test_the_console_entry_exists_in_pyproject():
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent / "pyproject.toml").read_text(encoding="utf-8")
    assert 'dna-identify = "dna_decode.identify.cli:main"' in src
