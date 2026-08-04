"""Flagellar motility catalog: pure decision logic + literature anchors + CLI (wheel-only, no deps)."""
from __future__ import annotations

import pytest

from dna_decode.motility.flagellar_catalog import (
    MotilityInputError,
    call_motility,
    catalog_genes,
)

# a full motile gene set (E. coli K-12 MG1655-like: all 5 modules + chemotaxis)
_MOTILE = ["flhD", "flhC", "fliA", "fliC", "motA", "motB", "fliF", "fliG", "flhA", "fliI",
           "cheA", "cheW", "cheY", "cheZ"]
_MOTILE_NO_CHE = ["flhD", "flhC", "fliA", "fliC", "motA", "motB", "fliF", "fliG", "flhA", "fliI"]


def test_full_set_is_motile_and_chemotactic():
    c = call_motility(_MOTILE)
    assert c.motile is True and c.verdict == "MOTILE"
    assert c.chemotaxis_competent is True
    assert c.missing_modules == ()


def test_swims_without_chemotaxis():
    # a che-mutant STILL SWIMS -- chemotaxis must NOT gate motility (biology)
    c = call_motility(_MOTILE_NO_CHE)
    assert c.motile is True
    assert c.chemotaxis_competent is False


@pytest.mark.parametrize("drop,module", [
    (["flhD", "flhC"], "master_regulator"),   # no master -> nothing transcribed
    (["fliC", "fljB"], "flagellin"),          # no filament
    (["motA", "motB"], "motor"),              # no motor
    (["fliA"], "sigma28"),                    # no class-3 sigma
    (["fliF"], "basal_export"),               # basal body incomplete
])
def test_knockout_of_a_module_is_nonmotile(drop, module):
    genes = [g for g in _MOTILE if g not in drop]
    c = call_motility(genes)
    assert c.motile is False
    assert module in c.missing_modules


def test_shigella_pattern_nonmotile():
    # Shigella flexneri lost its flagellar genes (pseudogenes) -> only chemotaxis remnants present
    c = call_motility(["cheA", "cheW", "cheY"])
    assert c.motile is False
    assert "master_regulator" in c.missing_modules and "flagellin" in c.missing_modules


def test_flagellin_is_any_fljB_alternate():
    # fljB (phase-2 flagellin) satisfies the flagellin module in place of fliC
    genes = [g for g in _MOTILE_NO_CHE if g != "fliC"] + ["fljB"]
    assert call_motility(genes).motile is True


def test_empty_input_raises():
    with pytest.raises(MotilityInputError):
        call_motility([])


def test_catalog_genes_covers_modules():
    g = catalog_genes()
    assert {"flhD", "flhC", "fliC", "motA", "motB", "fliA", "fliF", "cheA"} <= g


def test_cli_via_unified_dispatch(capsys):
    from dna_decode.cli import _delegate
    rc = _delegate("motility", ["--genes", ",".join(_MOTILE)])
    assert rc == 0
    assert "MOTILE" in capsys.readouterr().out
