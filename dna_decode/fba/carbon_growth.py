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


# Counter-ions, hydration states and salt forms that assay catalogues put in the label but BiGG does not.
# "Sodium succinate dibasic hexahydrate" and BiGG's "succinate" are the SAME carbon source; the assay name
# describes the reagent bottle, the model name describes the metabolite.
_SALT_HYDRATE_TOKENS: tuple[str, ...] = (
    "monohydrate", "dihydrate", "hexahydrate", "trihydrate", "hydrate", "anhydrous",
    "sodium salt", "disodium salt", "potassium salt", "dipotassium salt", "hydrochloride",
    "disodium", "dibasic", "monobasic", "sodium", "potassium", "salt",
)
# acid <-> conjugate-base: assay labels say "acid", BiGG says the anion
_ACID_TO_ANION: tuple[tuple[str, str], ...] = (
    ("gluconic acid", "gluconate"), ("glucuronic acid", "glucuronate"),
    ("galacturonic acid", "galacturonate"), ("glycolic acid", "glycolate"),
    ("malic acid", "malate"), ("ketoglutaric acid", "2-oxoglutarate"),
    ("lactic acid", "lactate"), ("pyruvic acid", "pyruvate"),
    ("succinic acid", "succinate"), ("acetic acid", "acetate"), ("fumaric acid", "fumarate"),
)


def strip_salt_and_hydrate(name: str) -> str:
    """PURE: drop counter-ion / hydration / salt wording from an assay carbon-source label.

    Load-bearing for joining an assay catalogue to a genome-scale model. Measured 2026-08-12 on the
    Fitness Browser's 28 Keio carbon sources: exact+alias matching resolved 14; **all 13 of the misses were
    NAME-gaps, not model-gaps** -- every one of their exchanges (EX_succ_e, EX_pyr_e, EX_ac_e, EX_mal__L_e,
    EX_malt_e, EX_tre_e, EX_glcn_e, EX_g6p_e, EX_akg_e, EX_galur_e, EX_glyclt_e, EX_lac__D_e, EX_lac__L_e)
    is PRESENT in iML1515. Only "casamino acids" is genuinely unmappable (an amino-acid mixture, not one
    exchange).

    Deliberately does NOT touch stereochemistry (`D-` / `L-`): D-lactate and L-lactate are different
    metabolites with different exchanges, and collapsing them would silently mis-assign the carbon source.
    """
    k = " ".join(name.strip().lower().split())
    for acid, anion in _ACID_TO_ANION:
        if acid in k:
            k = k.replace(acid, anion)
    # LONGEST-FIRST is load-bearing: replacing "sodium salt" before "disodium salt" leaves a stray "di"
    # ("L-Malic acid disodium salt monohydrate" -> "l-malate di"), which then matches nothing.
    for tok in sorted(_SALT_HYDRATE_TOKENS, key=len, reverse=True):
        k = k.replace(tok, " ")
    k = " ".join(k.split()).strip(" ,-")
    # greek-letter prefixes: assay labels write "a-Ketoglutaric", BiGG writes "2-oxoglutarate"
    for pre in ("alpha-", "a-", "beta-", "b-", "gamma-"):
        if k.startswith(pre) and len(k) > len(pre) + 2:
            k = k[len(pre):]
            break
    return k


def normalize_carbon_name(name: str) -> str:
    """PURE: a Biolog/assay carbon-source label -> a candidate BiGG metabolite-name key (lowercased)."""
    key = name.strip().lower()
    return _CARBON_ALIASES.get(key, key)


def build_exchange_name_index(model) -> dict[str, str]:
    """PURE-ish (reads model): {exchanged-metabolite-name(lower) -> EX_ reaction id} for carbon matching."""
    idx: dict[str, str] = {}
    formula_free: dict[str, set[str]] = {}
    for r in model.reactions:
        if r.id.startswith("EX_") and r.id.endswith("_e"):
            mets = list(r.metabolites)
            if mets:
                nm = mets[0].name.lower()
                idx[nm] = r.id
                # Many BiGG names carry a trailing molecular formula ("maltose c12h22o11",
                # "glycolate c2h3o3"). Index the formula-free form too, but ONLY where it is
                # unambiguous -- if two metabolites collapse to the same bare name, neither is added,
                # so we never silently bind a carbon source to the wrong exchange.
                bare = _strip_formula_suffix(nm)
                if bare != nm:
                    formula_free.setdefault(bare, set()).add(r.id)
    for bare, ids in formula_free.items():
        if len(ids) == 1 and bare not in idx:
            idx[bare] = next(iter(ids))
    return idx


def _strip_formula_suffix(name: str) -> str:
    """PURE: drop a trailing molecular-formula token from a BiGG metabolite name.

    'maltose c12h22o11' -> 'maltose'. A formula token is letters+digits with at least one digit and no
    vowel-word shape; conservative by design -- a real word like 'phosphate' is never stripped.
    """
    parts = name.split()
    if len(parts) < 2:
        return name
    last = parts[-1]
    if any(ch.isdigit() for ch in last) and all(ch.isalnum() for ch in last) and len(last) >= 3:
        letters = "".join(ch for ch in last if ch.isalpha())
        if letters and len(letters) <= 6:          # c12h22o11 -> "choo"; a formula, not a word
            return " ".join(parts[:-1])
    return name


def match_carbon_exchange(name: str, name_index: dict[str, str]) -> str | None:
    """PURE: resolve a carbon-source name to an EX_ id.

    Four passes, cheapest and most specific first: exact name -> curated alias -> salt/hydrate-stripped ->
    alias of the stripped form. Returns None rather than guessing, so an unmappable label (e.g. the
    "casamino acids" mixture) stays visibly unmapped instead of silently binding to the wrong exchange.
    """
    k = name.strip().lower()
    if k in name_index:
        return name_index[k]
    alias = normalize_carbon_name(name)
    if alias in name_index:
        return name_index[alias]
    stripped = strip_salt_and_hydrate(name)
    if stripped in name_index:
        return name_index[stripped]
    return name_index.get(_CARBON_ALIASES.get(stripped, stripped))


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
