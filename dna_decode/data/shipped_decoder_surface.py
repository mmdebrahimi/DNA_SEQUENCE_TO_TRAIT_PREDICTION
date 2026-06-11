"""Checked-in SHIPPED-DECODER SURFACE — the authoritative deployed-claim set for the validation report card.

Each entry is a (organism, drug) cell the suite makes a deployed claim about (or intentionally abstains on).
The validation report card derives its rows from THIS registry UNION the observed scored/census keys, so a
shipped decoder that has never been censused still renders (as NOT_CENSUSED) and a new decoder cannot ship
invisibly. A coverage test (`tests/test_shipped_decoder_surface.py`) asserts every CLI-routable drug appears
in >=1 row here.

This is the DEPLOYED-CLAIM set, NOT a cross-product (pre-exec /brainstorm C2 / Q1): it lists exactly the
organism×drug pairs the suite claims or abstains on — E. coli × {cipro,cef,tet,gent}; the calibrated-registry
organisms; S. aureus × oxacillin (label-confounded); and the other-kingdom decoders. It deliberately does
NOT contain E. coli × oxacillin (oxacillin is S. aureus/mecA-specific) or any other unclaimed pair.

`phenotype_source_status`:
- `ncbi_pd`          — bacterial cell scoreable on free NCBI-PD AST_phenotypes (a real validation cell).
- `label_confounded` — the phenotype LABEL is an unreliable surrogate (oxacillin AST vs mecA; cefoxitin is
                       the CLSI surrogate) -> report-card state LABEL_CONFOUNDED, distinct from NOT_CENSUSED.
- `no_free_source`   — no free isolate-level phenotype source exists (fungal/antiviral/antimalarial) ->
                       report-card state NO_FREE_PHENOTYPE_SOURCE (structural non-cell).
"""
from __future__ import annotations

# (organism, drug, engine, organism_scope, phenotype_source_status, census_group)
SHIPPED_DECODER_SURFACE: list[tuple[str, str, str, str, str, str | None]] = [
    # --- E. coli default DRUG_RULE bacterial decoders (NCBI-PD census-able) ---
    ("Escherichia_coli_Shigella", "ciprofloxacin", "amrfinder_curated", "Escherichia_coli", "ncbi_pd", "Escherichia_coli_Shigella"),
    ("Escherichia_coli_Shigella", "ceftriaxone",   "amrfinder_curated", "Escherichia_coli", "ncbi_pd", "Escherichia_coli_Shigella"),
    ("Escherichia_coli_Shigella", "tetracycline",  "amrfinder_curated", "Escherichia_coli", "ncbi_pd", "Escherichia_coli_Shigella"),
    ("Escherichia_coli_Shigella", "gentamicin",    "amrfinder_curated", "Escherichia_coli", "ncbi_pd", "Escherichia_coli_Shigella"),
    # --- Klebsiella: calibrated cipro + DRUG_RULE cef/tet/gent/meropenem (all NCBI-PD census-able) ---
    ("Klebsiella", "ciprofloxacin", "amrfinder_calibrated", "Klebsiella", "ncbi_pd", "Klebsiella"),
    ("Klebsiella", "ceftriaxone",   "amrfinder_curated",    "Klebsiella", "ncbi_pd", "Klebsiella"),
    ("Klebsiella", "tetracycline",  "amrfinder_curated",    "Klebsiella", "ncbi_pd", "Klebsiella"),
    ("Klebsiella", "gentamicin",    "amrfinder_curated",    "Klebsiella", "ncbi_pd", "Klebsiella"),
    ("Klebsiella", "meropenem",     "amrfinder_curated",    "Klebsiella", "ncbi_pd", "Klebsiella"),
    # --- calibrated registry: cipro on Campylobacter + Salmonella ---
    ("Campylobacter", "ciprofloxacin", "amrfinder_calibrated", "Campylobacter", "ncbi_pd", "Campylobacter"),
    ("Salmonella",    "ciprofloxacin", "amrfinder_calibrated", "Salmonella",    "ncbi_pd", "Salmonella"),
    # --- calibrated registry: meropenem abstainers ---
    ("Acinetobacter",          "meropenem", "amrfinder_calibrated", "Acinetobacter",          "ncbi_pd", "Acinetobacter"),
    ("Pseudomonas_aeruginosa", "meropenem", "amrfinder_calibrated", "Pseudomonas_aeruginosa", "ncbi_pd", "Pseudomonas_aeruginosa"),
    # --- S. aureus oxacillin: LABEL-CONFOUNDED (oxacillin AST unreliable; cefoxitin is the CLSI surrogate) ---
    ("Staphylococcus_aureus", "oxacillin", "amrfinder_curated", "Staphylococcus_aureus", "label_confounded", "Staphylococcus_aureus"),
    # --- other-kingdom decoders: NO free isolate-level phenotype source ---
    ("Candida_auris", "fluconazole",  "fungal_erg11", "Candida_auris", "no_free_source", None),
    ("Candida_auris", "voriconazole", "fungal_erg11", "Candida_auris", "no_free_source", None),
    ("Candida_auris", "caspofungin",  "fungal_fks1",  "Candida_auris", "no_free_source", None),
    ("Candida_auris", "micafungin",   "fungal_fks1",  "Candida_auris", "no_free_source", None),
    ("Influenza_A", "oseltamivir", "influenza_na", "Influenza_A", "no_free_source", None),
    ("Influenza_A", "peramivir",   "influenza_na", "Influenza_A", "no_free_source", None),
    ("Influenza_A", "zanamivir",   "influenza_na", "Influenza_A", "no_free_source", None),
    ("Plasmodium_falciparum", "artemisinin",        "pf_kelch13", "Plasmodium_falciparum", "no_free_source", None),
    ("Plasmodium_falciparum", "artesunate",         "pf_kelch13", "Plasmodium_falciparum", "no_free_source", None),
    ("Plasmodium_falciparum", "dihydroartemisinin", "pf_kelch13", "Plasmodium_falciparum", "no_free_source", None),
    ("Plasmodium_falciparum", "chloroquine",        "pf_pfcrt",   "Plasmodium_falciparum", "no_free_source", None),
]

_FIELDS = ("organism", "drug", "engine", "organism_scope", "phenotype_source_status", "census_group")


def shipped_decoder_rows() -> list[dict]:
    """The deployed-claim cells as dicts keyed by `_FIELDS`."""
    return [dict(zip(_FIELDS, row)) for row in SHIPPED_DECODER_SURFACE]


def surface_index() -> dict[tuple[str, str], dict]:
    """(organism.lower(), drug.lower()) -> surface row dict, for report-card classification."""
    return {(r["organism"].lower(), r["drug"].lower()): r for r in shipped_decoder_rows()}


def all_surface_drugs() -> set[str]:
    """Every drug named in the surface (lowercased) — used by the coverage test."""
    return {r["drug"].lower() for r in shipped_decoder_rows()}
