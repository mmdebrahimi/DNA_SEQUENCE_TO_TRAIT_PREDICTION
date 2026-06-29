"""IrisPlex 6-SNP eye-colour model (v0.1) — the DEPLOYED forensic rule, applied to independent data.

This is the eye-colour decoder's v0.1: the published IrisPlex/HIrisPlex multinomial-logistic model
(Walsh et al., IrisPlex FSI:Genetics 2011 5(3):170-180; HIrisPlex 2013) over 6 pigmentation SNPs.
v0 (`eye_colour.py`) was rs12913832 alone and ABSTAINED on heterozygotes; v0.1 adds the other 5 SNPs so a
heterozygote's eye colour can be resolved. It is the DEPLOYED model (validate-the-published-rule-on-
independent-data), NOT a classifier trained here — consistent with the project's deterministic-decoder
philosophy.

COEFFICIENTS ARE SOURCED, NOT FABRICATED (load-bearing rail). Pulled from the brianbhsu/eye-color open
implementation's `input/input.txt` (github.com/brianbhsu/eye-color), which encodes the published Walsh
HIrisPlex eye model; values match the canonical published model (rs12913832 brown beta 5.41 / intermediate
beta 3.16 is its signature). Structure: blue = reference category; two equations (intermediate, brown);
each SNP coded additively as the count of its effect allele.

STRAND HANDLING: 5 of 6 SNPs are non-palindromic (A/G, C/T, G/T) → strand-agnostic effect-allele counting
({effect, complement(effect)}) is unambiguous, matching v0's proven rs12913832 strand-agnosticism. The
6th, rs16891982, is C/G PALINDROMIC → strand-agnostic counting is impossible (C and G are each other's
complement AND both real alleles); we count the LITERAL effect allele C (the DTC forward-strand
convention, e.g. 23andMe) and flag `palindromic_assumed_forward_strand` — a named v0.1 caveat. Its
coefficient is moderate (0.53/1.46 vs rs12913832's 3.16/5.41), so a strand mismatch is bounded.
"""
from __future__ import annotations

from math import exp

# (intermediate, brown) betas; row 0 is the constant. SOURCED — see module docstring.
IRISPLEX_COEFFICIENTS: dict[str, tuple[float, float]] = {
    "constant": (-2.3640093, -2.6415884),
    "rs12913832": (3.1627512, 5.412669),
    "rs1800407": (-0.3869865, -1.3480642),
    "rs12896399": (-0.5080515, -0.7537442),
    "rs16891982": (0.5304902, 1.464204),
    "rs1393350": (-0.2088037, -0.4246789),
    "rs12203592": (-0.0019755, -0.6515579),
}

# effect allele per SNP (from the sourced table) + palindromic flag (C/G or A/T = strand-ambiguous).
_EFFECT: dict[str, str] = {
    "rs12913832": "A", "rs1800407": "T", "rs12896399": "T",
    "rs16891982": "C", "rs1393350": "A", "rs12203592": "T",
}
_PALINDROMIC: frozenset[str] = frozenset({"rs16891982"})  # C/G
_SNP_ORDER = ["rs12913832", "rs1800407", "rs12896399", "rs16891982", "rs1393350", "rs12203592"]
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
DEFAULT_THRESHOLD = 0.7  # Walsh et al. recommended category threshold


def _normalise(genotype: str) -> str:
    return "".join(c for c in genotype.upper() if c in "ACGT") if genotype else ""


def effect_allele_count(rsid: str, genotype: str) -> int | None:
    """Additive count (0/1/2) of the effect allele for `rsid`. None if genotype is not 2 callable bases.

    Strand-agnostic ({effect, complement}) for non-palindromic SNPs; literal effect allele for the
    palindromic rs16891982 (forward-strand assumption — see module docstring)."""
    g = _normalise(genotype)
    if len(g) != 2:
        return None
    eff = _EFFECT[rsid]
    targets = {eff} if rsid in _PALINDROMIC else {eff, _COMPLEMENT[eff]}
    return sum(1 for a in g if a in targets)


def predict_irisplex(genotypes: dict[str, str]) -> dict:
    """Predict eye colour from the 6 IrisPlex SNP genotypes (each a 2-allele string e.g. 'AG').

    Complete-case: rs12913832 is MANDATORY (the model is undefined without it) and ALL 6 must be callable
    (no imputation — imputing absent SNPs to reference would bias). Returns probabilities + argmax
    prediction + the 0.7-threshold category (blue/intermediate/brown/undefined).
    """
    counts: dict[str, int] = {}
    for rsid in _SNP_ORDER:
        c = effect_allele_count(rsid, genotypes.get(rsid, ""))
        if c is None:
            if rsid == "rs12913832":
                return {"prediction": "INDETERMINATE", "status": "MISSING_RS12913832"}
            return {"prediction": "INDETERMINATE", "status": "INCOMPLETE_SNP_SET", "missing": rsid}
        counts[rsid] = c

    vec = [1.0] + [float(counts[r]) for r in _SNP_ORDER]
    keys = ["constant"] + _SNP_ORDER
    lin_int = sum(IRISPLEX_COEFFICIENTS[k][0] * v for k, v in zip(keys, vec))
    lin_brown = sum(IRISPLEX_COEFFICIENTS[k][1] * v for k, v in zip(keys, vec))
    e_int, e_brown = exp(lin_int), exp(lin_brown)
    denom = 1.0 + e_int + e_brown
    p_blue, p_int, p_brown = 1.0 / denom, e_int / denom, e_brown / denom

    probs = {"blue": p_blue, "intermediate": p_int, "brown": p_brown}
    argmax = max(probs, key=probs.__getitem__)
    top = probs[argmax]
    category = argmax if top >= DEFAULT_THRESHOLD else "undefined"
    return {
        "prediction": argmax,
        "category_at_0.7": category,
        "p_blue": round(p_blue, 4), "p_intermediate": round(p_int, 4), "p_brown": round(p_brown, 4),
        "effect_allele_counts": counts,
        "status": "PREDICTED",
        "palindromic_assumed_forward_strand": ["rs16891982"],
    }
