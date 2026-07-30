"""General N-category multinomial-logistic pigmentation predictor — the engine the HIrisPlex hair (4-cat)
and HIrisPlex-S skin (5-cat) models plug into, generalizing the eye-colour `predict_eye_color`.

The math is the SAME multinomial-logistic softmax the IrisPlex eye model uses (irisplex.py), lifted to an
arbitrary number of categories + an externally-supplied coefficient table:

    for each non-reference category c:  Z_c = intercept_c + sum_i beta[c][i] * x_i
    P(reference) = 1 / D ;  P(c) = exp(Z_c) / D ;  D = 1 + sum_c exp(Z_c)

where x_i in {0,1,2} = count of SNP i's counted allele. The reference category has all-zero coefficients by
construction (its Z=0 → exp(0)=1, the `1` in D).

WHY a separate engine (not just more constants in irisplex.py): hair/skin need DOZENS of SNPs across 4/5
categories, and their published coefficients are NOT yet transcribed (image-encoded across Walsh 2017
Table 2 + its erratum + Walsh 2013 Table 3 — see wiki/pigment_hirisplex_coefficient_sourcing_2026-07-30.md).
This engine is validated NOW by reproducing the shipped, 1000G-population-validated EYE model exactly
(test_pigment_multinomial.py), so the ONLY remaining step for hair/skin is filling their `PigmentModel`
coefficient tables from the located sources — an attended transcription, never fabricated here.

Pure-python (math only), wheel-only, offline, deterministic. Regime-A curated catalog. Benign visible-trait
genetics, NOT a forensic tool. The frozen AMR/forward surfaces + the shipped eye model are untouched.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

_VALID_BASES = set("ACGT")


class MissingGenotypeError(ValueError):
    """A required SNP genotype is absent (never a silent wrong call)."""


@dataclass(frozen=True)
class PigmentSNP:
    """One predictor SNP: its rsID, the counted allele (same strand as the coefficients), and the per-
    non-reference-category beta coefficients {category_name: beta}."""
    rsid: str
    counted_allele: str
    betas: dict            # {non_reference_category: beta_value}
    required: bool = False  # a dominant SNP that must be present even under allow_missing (e.g. HERC2 for eye)


@dataclass(frozen=True)
class PigmentModel:
    """A multinomial-logistic pigmentation model. `categories[0]` is the REFERENCE category (all-zero
    coefficients). `intercepts` maps each NON-reference category to its intercept. `snps` carries per-SNP
    betas for each non-reference category. `coefficients_pending` flags a structure-only stub (hair/skin
    before transcription) so a consumer never mistakes an empty model for a real predictor."""
    trait: str
    categories: tuple[str, ...]            # categories[0] = reference
    intercepts: dict                       # {non_reference_category: intercept}
    snps: tuple[PigmentSNP, ...]
    source: str
    coefficients_pending: bool = False

    @property
    def reference(self) -> str:
        return self.categories[0]

    @property
    def non_reference(self) -> tuple[str, ...]:
        return self.categories[1:]


@dataclass
class PigmentPrediction:
    trait: str
    probabilities: dict          # {category: prob}
    call: str
    confidence: str              # high / medium / low
    counted_alleles: dict
    notes: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"trait": self.trait, "model": "multinomial_logistic", "regime": "A_curated_catalog",
                "probabilities": self.probabilities, "call": self.call, "confidence": self.confidence,
                "counted_alleles": self.counted_alleles, "notes": self.notes}


def _count_allele(genotype: str, allele: str) -> int:
    g = "".join(c for c in genotype.upper() if c in _VALID_BASES)
    if len(g) != 2:
        raise MissingGenotypeError(f"genotype {genotype!r} is not a diploid A/C/G/T call")
    return g.count(allele)


def predict(model: PigmentModel, genotypes: dict, *, allow_missing: bool = False) -> PigmentPrediction:
    """Predict the categorical pigmentation phenotype from {rsID: genotype-string}.

    - Requires every SNP by default (missing → MissingGenotypeError); a SNP flagged `required` is ALWAYS
      required even under allow_missing. `allow_missing=True` imputes a non-required missing SNP as x=0
      (biased; caps confidence low).
    - Genotypes are counted on the SAME strand as each SNP's `counted_allele` (strand harmonization for
      real array data is the caller's responsibility — see the 1000G validator).
    """
    if model.coefficients_pending:
        raise ValueError(f"{model.trait} model coefficients are not yet populated (stub); see "
                         "wiki/pigment_hirisplex_coefficient_sourcing_2026-07-30.md — refusing to predict "
                         "from an empty coefficient table (never fabricated).")
    norm = {k.lower().strip(): v for k, v in genotypes.items()}
    notes: list = []
    counted: dict = {}
    z = {c: model.intercepts[c] for c in model.non_reference}
    for snp in model.snps:
        key = snp.rsid.lower()
        if key not in norm or norm[key] in (None, "", "--", "NN"):
            if snp.required:
                raise MissingGenotypeError(f"{snp.rsid} (a required predictor) is missing")
            if not allow_missing:
                raise MissingGenotypeError(f"missing genotype for {snp.rsid}; pass allow_missing=True to impute x=0")
            x = 0
            notes.append(f"{snp.rsid} missing → imputed x=0 (allow_missing); confidence capped low")
        else:
            x = _count_allele(norm[key], snp.counted_allele)
        counted[snp.rsid] = x
        for c in model.non_reference:
            z[c] += snp.betas.get(c, 0.0) * x

    exps = {c: math.exp(z[c]) for c in model.non_reference}
    d = 1.0 + sum(exps.values())
    probs = {model.reference: 1.0 / d}
    for c in model.non_reference:
        probs[c] = exps[c] / d
    call = max(probs, key=probs.get)
    top = probs[call]
    if any("imputed" in n for n in notes):
        conf = "low"
    elif top >= 0.70:
        conf = "high"
    elif top >= 0.50:
        conf = "medium"
    else:
        conf = "low"
    return PigmentPrediction(model.trait, {k: round(v, 6) for k, v in probs.items()}, call, conf,
                             counted, notes)


def reference_integrity_ok(model: PigmentModel, anchors: list[tuple[dict, str]]) -> bool:
    """Biology contract guard: every (genotype, expected_top_category) anchor must predict that category.
    A corrupted/fabricated coefficient table fails this loudly. `anchors` are known genotype→phenotype
    directions (e.g. for eye: HERC2 GG → blue)."""
    for genotype, expected in anchors:
        if predict(model, genotype).call != expected:
            return False
    return True
