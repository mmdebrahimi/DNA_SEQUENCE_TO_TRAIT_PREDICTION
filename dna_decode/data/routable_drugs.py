"""The ONE definition of which drugs `dna-amr --drug` accepts.

WHY THIS MODULE EXISTS. This union was written out twice -- once as the CLI's argparse `choices`, once
inside `cell_registry.cli_routable_manifest()` -- and the two drifted. HCMV shipped 2026-07-23 with five
CLI-routable drugs and was added to the CLI copy but not the registry copy, so for weeks the registry
believed those drugs were not routable, the coverage test could not see them, and five shipped decoders
rendered with NO evidence contract at all. The registry's own comment two lines below that union boasts
that `traits` is "DERIVED from `per_target`, never hand-listed" -- while the drug union directly above it
was hand-listed.

That is the fifth instance in this repo of a hand-enumerated list drifting from the data that defines it
([[feedback_hardcoded_exclusion_list_undercovers]]). A guard whose input set is hand-maintained cannot
catch the case where someone forgets to maintain it. Both call sites now import from here, so a seventh
catalog cannot be added to one and missed by the other.

Pure; imports only the per-kingdom catalogs. No I/O.
"""
from __future__ import annotations


def all_routable_amr_drugs() -> set[str]:
    """Every drug `dna-amr --drug` routes, across all kingdoms. The CLI's choices ARE this set."""
    from dna_decode.data.antimalarial_amr import supported_antimalarial_drugs
    from dna_decode.data.antiviral_amr import supported_antiviral_drugs
    from dna_decode.data.fungal_amr import supported_fungal_drugs
    from dna_decode.data.hcmv_amr import all_supported_hcmv_drugs
    from dna_decode.data.hiv_amr import all_supported_hiv_drugs
    from dna_decode.data.mic_tiers import supported_drugs
    from dna_decode.data.sarscov2_amr import all_supported_sarscov2_drugs

    return (set(supported_drugs())
            | set(supported_fungal_drugs())
            | set(supported_antimalarial_drugs())
            | set(supported_antiviral_drugs())
            | set(all_supported_hiv_drugs())
            | set(all_supported_sarscov2_drugs())
            | set(all_supported_hcmv_drugs()))


def all_routable_tb_drugs() -> set[str]:
    """Every drug `dna-tb --drug` routes: the WHO catalogue's own long drug names.

    SEPARATE from `all_routable_amr_drugs` because TB is a SEPARATE console script (`dna-tb`) on a
    separate route -- it lives in the NON-frozen `organism_rules` package and takes a VCF against H37Rv
    rather than a genome FASTA or `--observed`. Folding it into the `dna-amr` union would claim those
    drugs are reachable through a command that rejects them.

    The catalogue's `DRUG_CATALOGUE_NAME` carries BOTH long names and short aliases (`inh`, `rif`), so
    the long form is selected by length. Short aliases stay accepted at the CLI but are not separate
    cells -- `inh` and `isoniazid` are one decoder, and registering both would double-count the surface.

    Reading the CATALOGUE rather than the contracts is deliberate: it is what the CLI validates against,
    so the registry's coverage test compares two independent sets and can actually fail.
    """
    from dna_decode.data.tb_who_catalogue import DRUG_CATALOGUE_NAME

    return {d for d in DRUG_CATALOGUE_NAME if len(d) > 4}
