"""AR Isolate Bank organism registry — the config that generalizes the gono validation harness to any
organism with an existing NON-FROZEN organism_rules cell.

Each entry maps an organism key -> {needle (AR Bank organism substring), amrfinder_organism (AMRFinder -O),
registry_organism (display), drug_map (canonical_drug -> (AR-Bank-page-shown-name, call_function)), note}.
The generalized scorer (`scripts/ar_bank_organism_validate.py`) + label builder
(`scripts/build_ar_bank_organism_labels.py`) read from here, so validating a new organism is a registry
edit, not new code. Frozen decoder surface is untouched — every call_function is a curated NON-FROZEN cell.
"""
from __future__ import annotations

from dna_decode.organism_rules import enterococcus_amr as _efm
from dna_decode.organism_rules import staphylococcus_amr as _sa
from dna_decode.organism_rules.neisseria_amr import call_ng_amr as _ng_amr

# organism_key -> config
AR_BANK: dict[str, dict] = {
    "enterococcus_faecium": {
        "needle": "Enterococcus faecium",
        "amrfinder_organism": "Enterococcus_faecium",
        "registry_organism": "Enterococcus faecium",
        # canonical_drug -> (AR-Bank INT column name, call_function)
        "drug_map": {
            # Levofloxacin (a fluoroquinolone) shares the enterococcal FQ mechanism (gyrA/parC QRDR) with
            # the cell's ciprofloxacin rule -- a documented cross-drug (same-mechanism) application.
            "levofloxacin": ("Levofloxacin", _efm.call_efm_ciprofloxacin),
            # Doxycycline (a tetracycline) <- the acquired-tet-gene rule (tet(M)/(L)/(S)/(O)).
            "doxycycline": ("Doxycycline", _efm.call_efm_tetracycline),
            # High-Level Gentamicin has NO S/I/R INT on the AR Bank page -> unscorable (dropped).
        },
        "note": ("efm FQ-QRDR + acquired-tet rules; Levofloxacin<-cipro-rule (same gyrA/parC QRDR), "
                 "Doxycycline<-tet-rule (tet-class). HLG unscorable (no INT)."),
    },
    "enterococcus_faecalis": {
        "needle": "Enterococcus faecalis",
        "amrfinder_organism": "Enterococcus_faecalis",
        "registry_organism": "Enterococcus faecalis",
        "drug_map": {
            # The efm FQ-QRDR + tet-gene mechanisms are pan-enterococcal (conserved QRDR codons, acquired
            # tet genes) -> applied to faecalis as a documented cross-species mechanistic transfer.
            "levofloxacin": ("Levofloxacin", _efm.call_efm_ciprofloxacin),
            "doxycycline": ("Doxycycline", _efm.call_efm_tetracycline),
        },
        "note": "pan-enterococcal FQ-QRDR + tet-gene rules applied to faecalis (mechanism transfer; efm cell).",
    },
    "staphylococcus_aureus": {
        "needle": "Staphylococcus aureus",
        "amrfinder_organism": "Staphylococcus_aureus",
        "registry_organism": "Staphylococcus aureus",
        "drug_map": {
            # Only Levofloxacin is scorable: the panel has NO Rifampin (so call_sa_rifampicin has no label)
            # and NO ciprofloxacin -- Levofloxacin (FQ) <- the cipro gyrA/grlA QRDR rule (same mechanism).
            "levofloxacin": ("Levofloxacin", _sa.call_sa_ciprofloxacin),
        },
        "note": ("only Levofloxacin scorable (panel has no Rifampin/ciprofloxacin); Levofloxacin<-cipro "
                 "gyrA/grlA QRDR rule. ~70 FREE genomes -> Kaggle-native AMRFinder (local Docker wedges >50)."),
    },
    "gono": {
        "needle": "gonorrhoeae",
        "amrfinder_organism": "Neisseria_gonorrhoeae",
        "registry_organism": "Neisseria gonorrhoeae",
        # gono has its own multi-drug dispatcher (call_ng_amr); expose it uniformly.
        "drug_map": {
            d: (d.capitalize(), (lambda dr: (lambda syms: _ng_amr(dr, syms)))(d))
            for d in ("azithromycin", "cefixime", "ceftriaxone", "ciprofloxacin", "penicillin", "tetracycline")
        },
        "note": "gono cell (call_ng_amr) surfaced through the uniform registry.",
    },
}


def config_for(organism_key: str) -> dict:
    if organism_key not in AR_BANK:
        raise KeyError(f"unknown AR-Bank organism {organism_key!r}; known: {sorted(AR_BANK)}")
    return AR_BANK[organism_key]


def rule_fn_for(organism_key: str):
    """Return a uniform `rule_fn(drug, symbols) -> {'prediction': ...}` for the organism, dispatching to
    the registered per-drug call functions. An unregistered drug ABSTAINS (INDETERMINATE)."""
    dm = {c.lower(): fn for c, (_shown, fn) in config_for(organism_key)["drug_map"].items()}

    def rule_fn(drug: str, symbols: list[str]) -> dict:
        fn = dm.get((drug or "").strip().lower())
        if fn is None:
            return {"prediction": "INDETERMINATE", "rule": f"no {organism_key} rule for {drug!r}",
                    "rule_status": "ABSTAIN_UNSUPPORTED_DRUG", "rule_scope": "scorer_local"}
        return fn(symbols)

    return rule_fn


def shown_name_map(organism_key: str) -> dict[str, str]:
    """canonical_drug -> the AR Bank page's INT column name (for the label builder)."""
    return {c: shown for c, (shown, _fn) in config_for(organism_key)["drug_map"].items()}
