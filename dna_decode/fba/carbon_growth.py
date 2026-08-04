"""FBA carbon-source growth: set a sole carbon source in the model's medium -> predicted growth rate.

The quantitative-trait complement to the binary essentiality axis: for a carbon source, swap it into the
model's (known-growing) default medium as the sole carbon and optimize biomass -> a growth RATE (/h). This
is the FBA analog of a Biolog carbon-utilization assay.

**Honest validation scope (see wiki/fba_carbon_growth_validation_2026-08-03.md):** a clean *measured
growth-rate* dataset across many carbon sources does NOT exist fetchably, and Biolog pos+neg data for the
iML1515 K-12 strain is SI-locked (the 190-source Biolog set is E. coli Nissle, a strain mismatch). So the
reachable validation is **RECALL** on measured-positive K-12 carbon sources (the Keio/Wetmore assay only ran
on sources E. coli grows on); full specificity + a growth-rate correlation are EXTERNAL-walled.

Pure helpers here (medium-swap logic, name normalization); the cobra optimize lives behind `predict_growth`.
"""
from __future__ import annotations

# common carbon-source-name -> BiGG metabolite-name aliases (exact-name matching misses these)
_CARBON_ALIASES: dict[str, str] = {
    "d-maltose": "maltose", "d-gluconic acid": "d-gluconate", "d-glucuronic acid": "d-glucuronate",
    "d-galacturonic acid": "d-galacturonate", "d-glucose-6-phosphate": "d-glucose 6-phosphate",
    "d-glucosamine hydrochloride": "d-glucosamine", "a-ketoglutaric": "2-oxoglutarate",
    "n-acetyl-d-glucosamine": "n-acetyl-d-glucosamine", "l-fucose": "l-fucose",
    "d-mannitol": "d-mannitol", "d-sorbitol": "d-sorbitol", "d-mannose": "d-mannose",
}


def normalize_carbon_name(name: str) -> str:
    """PURE: a Biolog/assay carbon-source label -> a candidate BiGG metabolite-name key (lowercased)."""
    key = name.strip().lower()
    return _CARBON_ALIASES.get(key, key)


def build_exchange_name_index(model) -> dict[str, str]:
    """PURE-ish (reads model): {exchanged-metabolite-name(lower) -> EX_ reaction id} for carbon matching."""
    idx: dict[str, str] = {}
    for r in model.reactions:
        if r.id.startswith("EX_") and r.id.endswith("_e"):
            mets = list(r.metabolites)
            if mets:
                idx[mets[0].name.lower()] = r.id
    return idx


def match_carbon_exchange(name: str, name_index: dict[str, str]) -> str | None:
    """PURE: resolve a carbon-source name to an EX_ id via exact name, then the alias, then None."""
    k = name.strip().lower()
    if k in name_index:
        return name_index[k]
    alias = normalize_carbon_name(name)
    return name_index.get(alias)


def predict_growth(model, carbon_exchange: str, base_medium: dict | None = None,
                   default_carbon: str = "EX_glc__D_e", uptake: float = 10.0) -> float:
    """Growth rate (/h) with `carbon_exchange` swapped in as the sole carbon source.

    Uses `model.medium` (the known-growing default) minus the default carbon plus the target -> keeps every
    inorganic component correct (the robust way; a zero-all-then-reopen approach silently cuts a nutrient).
    Returns 0.0 for no growth / an infeasible medium.
    """
    base = dict(base_medium if base_medium is not None else model.medium)
    med = dict(base)
    med.pop(default_carbon, None)
    med[carbon_exchange] = uptake
    with model:
        try:
            model.medium = med
        except Exception:
            return 0.0  # exchange not in this model
        val = model.slim_optimize()
    return 0.0 if (val is None or val != val or val < 0) else float(val)
