"""Offline tests for the Arabidopsis flowering-habit cell. Pure rule logic; no network/no D:/no GPU.

Pins the literature anchors — especially the one a NAIVE "FRI decides" rule gets WRONG (Da(1)-12 /
Shakhdara: functional FRI yet summer-annual, via a weak FLC). That contract is the anti-fabrication guard.
"""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.organism_rules.arabidopsis_flowering import (  # noqa: E402
    FloweringInputError, call_flowering_habit, reference_integrity_ok, status_for,
)
from dna_decode.organism_rules.flowering_cli import main as flowering_main  # noqa: E402
from dna_decode import cli as unified  # noqa: E402


def test_reference_integrity_biology_contract():
    # THE load-bearing guard: a corrupted catalog/rule fails this.
    assert reference_integrity_ok() is True


def test_col0_is_rapid_cycler():
    c = call_flowering_habit("Col", "Col")          # FRI-LoF + strong FLC
    assert c.habit == "summer_annual_early" and c.vernalization_required is False


def test_col_fri_sf2_is_winter_annual():
    c = call_flowering_habit("Sf-2", "Col")         # the classic late-flowering Col-FRI line
    assert c.habit == "winter_annual_late" and c.vernalization_required is True and c.confidence == "high"


def test_flc_route_beats_naive_fri_only_rule():
    # Da(1)-12 has a FUNCTIONAL FRI yet is a summer annual (weak FLC). A naive "FRI decides" rule says LATE.
    c = call_flowering_habit("Sf-2", "Da(1)-12")
    assert c.habit == "summer_annual_early" and c.fri_status == "functional" and c.flc_status == "weak"
    assert any("naive FRI-only rule" in n for n in c.notes)


def test_flc_null_early_despite_functional_fri():
    for flc in ("Van-0", "Bur-0"):
        c = call_flowering_habit("Sf-2", flc)
        assert c.habit == "summer_annual_early" and c.confidence == "high"


def test_fri_route_confidence_capped_by_lz0_counterexample():
    c = call_flowering_habit("Col", "Col")
    assert c.confidence == "medium"                 # FRI-LoF route is NOT high — Lz-0 can still be late
    assert any("Lz-0" in n for n in c.notes)


def test_ler_double_hit():
    c = call_flowering_habit("Ler", "Ler")          # FRI deletion + weak FLC
    assert c.habit == "summer_annual_early" and c.flc_status == "weak"


def test_unknown_locus_abstains():
    c = call_flowering_habit("unknown", "Col")
    assert c.habit == "ABSTAIN" and c.vernalization_required is None and c.confidence == "low"


def test_unknown_allele_name_raises():
    with pytest.raises(FloweringInputError):
        call_flowering_habit("NotARealAllele", "Col")


def test_status_passthrough():
    assert status_for("FRI", "functional") == "functional"
    assert status_for("FLC", "Ler") == "weak"


def test_undetectable_mechanisms_surfaced():
    d = call_flowering_habit("Col", "Col").as_dict()
    assert d["undetectable_mechanisms"] and "PARTIAL" in d["scope_limit"]


def test_cli_human_and_json(capsys):
    rc = flowering_main(["--fri", "Sf-2", "--flc", "Col"])
    out = capsys.readouterr().out
    assert rc == 0 and "WINTER_ANNUAL_LATE" in out
    rc = flowering_main(["--fri", "Col", "--flc", "Col", "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["habit"] == "summer_annual_early" and d["organism"] == "Arabidopsis_thaliana"


def test_cli_list_alleles(capsys):
    rc = flowering_main(["--list-alleles"])
    out = capsys.readouterr().out
    assert rc == 0 and "FRI:" in out and "FLC:" in out and "Sf-2" in out and "source:" in out


def test_cli_bad_allele_exits_2(capsys):
    rc = flowering_main(["--fri", "Bogus", "--flc", "Col"])
    assert rc == 2 and "error" in capsys.readouterr().err


def test_dispatch_through_unified_cli(capsys):
    rc = unified.main(["flowering", "--fri", "Sf-2", "--flc", "Col", "--json"])
    assert rc == 0
    assert json.loads(capsys.readouterr().out)["habit"] == "winter_annual_late"


def test_flowering_in_registry_and_list(capsys):
    assert "flowering" in unified.TRAITS
    rc = unified.main(["list"])
    out = capsys.readouterr().out
    assert rc == 0 and "flowering" in out and "arabidopsis" in out.lower()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
