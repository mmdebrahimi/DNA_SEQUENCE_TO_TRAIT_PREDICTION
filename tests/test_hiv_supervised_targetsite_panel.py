"""Tests for the PI/INSTI supervised panel (offline; real run skip-guarded)."""
from __future__ import annotations

import importlib.util, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))
pytest.importorskip("sklearn")

spec = importlib.util.spec_from_file_location("m", REPO / "scripts" / "hiv_supervised_targetsite_panel.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
import dna_decode.data.hiv_amr as H


def test_genes_wired_to_hiv_amr_classes():
    assert mod.GENES["PI"]["positions"] == set(H.PI_CLASS.positions)
    assert mod.GENES["INSTI"]["positions"] == set(H.INSTI_CLASS.positions)
    assert mod.GENES["PI"]["cols"] == ("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV")
    assert mod.GENES["INSTI"]["cols"] == ("RAL", "EVG", "DTG", "BIC", "CAB")


@pytest.mark.skipif(not (REPO / "data" / "hiv_ref" / "HIV1_PR_HXB2_cds.fna").exists(),
                    reason="needs committed PR reference")
def test_translate_protease_length():
    prot = mod.translate(REPO / "data" / "hiv_ref" / "HIV1_PR_HXB2_cds.fna")
    assert len(prot) == 99 and prot[29] == "D"      # protease is 99 aa; D30 wild-type (PI major position 30)


@pytest.mark.skipif(not (REPO / "data" / "hiv_ref" / "HIV1_IN_HXB2_cds.fna").exists(),
                    reason="needs committed IN reference")
def test_translate_integrase_length():
    prot = mod.translate(REPO / "data" / "hiv_ref" / "HIV1_IN_HXB2_cds.fna")
    assert len(prot) == 288 and prot[147] == "Q"    # integrase is 288 aa; Q148 wild-type (INSTI major position 148)
