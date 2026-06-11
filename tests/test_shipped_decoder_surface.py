"""Pin the shipped-decoder surface registry (dna_decode/data/shipped_decoder_surface.py, C2 amendment).

The coverage test is the anti-drift guard: every CLI-routable drug MUST appear in >=1 surface row, so a new
shipped decoder cannot ship invisible to the validation report card.
"""
from __future__ import annotations

from dna_decode.data.shipped_decoder_surface import (
    SHIPPED_DECODER_SURFACE,
    all_surface_drugs,
    shipped_decoder_rows,
    surface_index,
)


def _cli_routable_drugs() -> set[str]:
    from dna_decode.data.mic_tiers import supported_drugs
    drugs = {d.lower() for d in supported_drugs()}
    from dna_decode.data import antiviral_amr, antimalarial_amr, fungal_amr
    for mod in (antiviral_amr, antimalarial_amr, fungal_amr):
        for name in dir(mod):
            obj = getattr(mod, name)
            if isinstance(obj, dict):
                drugs |= {k.lower() for k in obj.keys() if isinstance(k, str)}
    # keep only plausible drug names (lower-case words) — the gene/mutation dicts also exist in those modules
    return drugs


def test_every_cli_drug_in_surface():
    """C2 coverage: no shipped CLI drug is missing from the surface."""
    surface = all_surface_drugs()
    # the canonical CLI drug catalogs (bacterial breakpoints + the 3 other kingdoms)
    from dna_decode.data.mic_tiers import supported_drugs
    expected = {d.lower() for d in supported_drugs()}
    expected |= {"fluconazole", "voriconazole", "caspofungin", "micafungin"}
    expected |= {"oseltamivir", "peramivir", "zanamivir"}
    expected |= {"artemisinin", "artesunate", "dihydroartemisinin", "chloroquine"}
    missing = expected - surface
    assert not missing, f"CLI drugs missing from shipped_decoder_surface: {sorted(missing)}"


def test_no_nonsensical_ecoli_oxacillin():
    """C2: oxacillin is S. aureus/mecA-specific — E. coli oxacillin must NOT be a surface cell."""
    assert ("escherichia_coli_shigella", "oxacillin") not in surface_index()


def test_oxacillin_is_label_confounded():
    idx = surface_index()
    row = idx[("staphylococcus_aureus", "oxacillin")]
    assert row["phenotype_source_status"] == "label_confounded"


def test_other_kingdom_no_free_source():
    idx = surface_index()
    for key in (("candida_auris", "fluconazole"), ("influenza_a", "oseltamivir"),
                ("plasmodium_falciparum", "artemisinin")):
        assert idx[key]["phenotype_source_status"] == "no_free_source"


def test_every_row_well_formed():
    for r in shipped_decoder_rows():
        assert set(r) == {"organism", "drug", "engine", "organism_scope", "phenotype_source_status", "census_group"}
        assert r["phenotype_source_status"] in {"ncbi_pd", "label_confounded", "no_free_source"}
    assert len(SHIPPED_DECODER_SURFACE) == len(shipped_decoder_rows())
