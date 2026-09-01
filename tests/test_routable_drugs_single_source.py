"""The CLI's `--drug` choices and the registry's routable set must be ONE definition.

THE DEFECT THIS PINS. These two unions were spelled out separately. HCMV shipped 2026-07-23 with five
CLI-routable drugs, was added to the CLI copy, and missed in the registry copy. For weeks the registry
believed those drugs were not routable, so `test_every_cli_amr_drug_has_a_contract` -- which exists
precisely to stop a decoder shipping without an evidence contract -- could not see them, and five cells
rendered with no contract at all.

The guard was correct the whole time; its INPUT was wrong. A coverage test whose input set is
hand-maintained cannot catch "someone forgot to maintain the input set". Fifth instance of that bug class
in this repo.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.data.cell_registry import cells, cli_routable_manifest  # noqa: E402
from dna_decode.data.routable_drugs import all_routable_amr_drugs  # noqa: E402


def test_the_cli_parser_and_the_registry_read_the_same_union():
    """Not "they happen to agree" -- they must be the SAME call. Resolved through the real parser."""
    import argparse
    from unittest.mock import patch

    seen = {}
    real = argparse.ArgumentParser.add_argument

    def spy(self, *a, **kw):
        if "--drug" in a and kw.get("choices"):
            seen["choices"] = set(kw["choices"])
        return real(self, *a, **kw)

    from dna_decode.amr import cli as amr_cli
    with patch.object(argparse.ArgumentParser, "add_argument", spy):
        try:
            amr_cli.main([])          # builds the parser, then exits on missing args
        except SystemExit:
            pass
    assert seen.get("choices"), "could not observe the --drug choices through the real parser"
    assert seen["choices"] == all_routable_amr_drugs(), (
        "the CLI's --drug choices diverged from all_routable_amr_drugs()")


def test_every_kingdom_catalog_reaches_the_routable_set():
    """NON-VACUITY. A union that silently dropped a catalog is exactly the defect; check each arrives."""
    from dna_decode.data.antimalarial_amr import supported_antimalarial_drugs
    from dna_decode.data.antiviral_amr import supported_antiviral_drugs
    from dna_decode.data.fungal_amr import supported_fungal_drugs
    from dna_decode.data.hcmv_amr import all_supported_hcmv_drugs
    from dna_decode.data.hiv_amr import all_supported_hiv_drugs
    from dna_decode.data.mic_tiers import supported_drugs
    from dna_decode.data.sarscov2_amr import all_supported_sarscov2_drugs

    routable = all_routable_amr_drugs()
    for name, cat in [("bacterial", supported_drugs), ("fungal", supported_fungal_drugs),
                      ("antimalarial", supported_antimalarial_drugs),
                      ("influenza", supported_antiviral_drugs), ("hiv", all_supported_hiv_drugs),
                      ("sarscov2", all_supported_sarscov2_drugs), ("hcmv", all_supported_hcmv_drugs)]:
        drugs = set(cat())
        assert drugs, f"{name} catalog is empty — this check would be vacuous for it"
        assert drugs <= routable, f"{name} drugs missing from the routable union: {sorted(drugs - routable)}"


def test_hcmv_is_contracted_now_the_regression_that_started_this():
    """The specific cells that were invisible. All five, on the viral track, with a real tier."""
    hcmv = {c.target: c for c in cells() if c.organism == "HCMV"}
    from dna_decode.data.hcmv_amr import all_supported_hcmv_drugs
    assert set(hcmv) == set(all_supported_hcmv_drugs()), sorted(hcmv)
    for c in hcmv.values():
        assert c.track == "viral" and c.route == "dna-amr"
        assert c.evidence_tier.value != "not_censused", f"{c.target} still reads as never-validated"


def test_the_registry_manifest_matches_the_shared_union():
    assert cli_routable_manifest()["dna-amr"] == {d.lower() for d in all_routable_amr_drugs()}
