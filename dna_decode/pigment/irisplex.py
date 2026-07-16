"""IrisPlex eye-colour predictor — deterministic multinomial-logistic model over 6 curated SNPs.

The model (Walsh et al. 2011, "IrisPlex: A sensitive DNA tool for accurate prediction of blue and brown eye
colour", Forensic Sci Int Genet 5(3):170-180) classifies an individual as blue / intermediate / brown from
the number of a specified allele at each of 6 pigmentation SNPs, via multinomial logistic regression with
**blue as the reference category**:

    Z_intermediate = c_int + sum_i (beta_int[i]  * x_i)
    Z_brown        = c_brown + sum_i (beta_brown[i] * x_i)
    D = 1 + exp(Z_intermediate) + exp(Z_brown)
    P(blue)         = 1 / D
    P(intermediate) = exp(Z_intermediate) / D
    P(brown)        = exp(Z_brown) / D

where x_i in {0,1,2} = the count of the SNP's counted allele in the genotype.

COEFFICIENT PROVENANCE (reference-integrity — NOT fabricated): the coefficient table below is transcribed
VERBATIM from the machine-readable reimplementation `brianbhsu/eye-color` (`input/input.txt`,
https://github.com/brianbhsu/eye-color), which implements the Walsh 2011 IrisPlex model. The CANONICAL source
is Walsh 2011 Table / the FROG-kb tool (frog.med.yale.edu) / the HIrisPlex webtool
(hirisplex.erasmusmc.nl) — a v0.1 task is to re-verify these decimals against Walsh 2011 directly. The
biology contract (rs12913832 AA -> brown, GG -> blue) is asserted by `reference_integrity_ok()` + the tests,
so a corrupted/fabricated coefficient set fails loudly rather than mis-calling silently.

Pure-python (math only), wheel-only, offline. Deterministic. This is Regime-A/curated-catalog, NOT a learned
embedding. Scope: benign visible-trait genetics, NOT a forensic tool.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# (rsID, counted_allele, beta_intermediate, beta_brown) — verbatim from brianbhsu/eye-color input.txt.
# The "Constant" row is the per-category intercept.
_CONSTANT = (-2.3640093, -2.6415884)   # (c_intermediate, c_brown)
IRISPLEX_SNPS: tuple[tuple[str, str, float, float], ...] = (
    ("rs12913832", "A", 3.1627512, 5.4126690),   # HERC2 — dominant predictor
    ("rs1800407", "T", -0.3869865, -1.3480642),  # OCA2
    ("rs12896399", "T", -0.5080515, -0.7537442),  # SLC24A4
    ("rs16891982", "C", 0.5304902, 1.4642040),   # SLC45A2
    ("rs1393350", "A", -0.2088037, -0.4246789),  # TYR
    ("rs12203592", "T", -0.0019755, -0.6515579),  # IRF4
)
_VALID_BASES = set("ACGT")


class MissingGenotypeError(ValueError):
    """Raised when a required IrisPlex SNP genotype is absent (never a silent wrong call)."""


@dataclass
class IrisPlexResult:
    p_blue: float
    p_intermediate: float
    p_brown: float
    call: str                 # "blue" | "intermediate" | "brown"
    confidence: str           # "high" | "medium" | "low"
    counted_alleles: dict     # rsID -> x (0/1/2), for audit
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "trait": "eye_colour", "model": "IrisPlex", "regime": "A_curated_catalog",
            "p_blue": self.p_blue, "p_intermediate": self.p_intermediate, "p_brown": self.p_brown,
            "call": self.call, "confidence": self.confidence,
            "counted_alleles": self.counted_alleles, "notes": self.notes,
        }


def _count_allele(genotype: str, allele: str) -> int:
    """Count occurrences of `allele` in a diploid genotype string ('AA' / 'A/G' / 'A|G')."""
    g = "".join(c for c in genotype.upper() if c in _VALID_BASES)
    if len(g) != 2:
        raise MissingGenotypeError(f"genotype {genotype!r} is not a diploid A/C/G/T call")
    return g.count(allele)


def predict_eye_color(genotypes: dict, *, allow_missing: bool = False) -> IrisPlexResult:
    """Predict eye colour from a dict {rsID: genotype-string} for the 6 IrisPlex SNPs.

    - By default REQUIRES all 6 SNPs (missing any -> MissingGenotypeError; never a silent wrong call).
    - `allow_missing=True` imputes a missing SNP as x=0 (counted-allele absent) and flags low confidence —
      use only when you accept the bias (e.g. a DTC array missing one probe). rs12913832 (the dominant SNP)
      is ALWAYS required even under allow_missing.
    - Genotypes are counted on the SAME strand as the coefficient table's allele; strand harmonization for
      real DTC data (openSNP) is a documented v0.1 follow-on.
    """
    notes: list[str] = []
    norm = {k.lower().strip(): v for k, v in genotypes.items()}
    z_int = _CONSTANT[0]
    z_brown = _CONSTANT[1]
    counted: dict = {}
    for rsid, allele, b_int, b_brown in IRISPLEX_SNPS:
        if rsid not in norm or norm[rsid] in (None, "", "--", "NN"):
            if rsid == "rs12913832":
                raise MissingGenotypeError("rs12913832 (HERC2, the dominant IrisPlex SNP) is required")
            if not allow_missing:
                raise MissingGenotypeError(f"missing genotype for {rsid}; pass allow_missing=True to impute x=0")
            x = 0
            notes.append(f"{rsid} missing -> imputed x=0 (allow_missing); confidence capped low")
        else:
            x = _count_allele(norm[rsid], allele)
        counted[rsid] = x
        z_int += b_int * x
        z_brown += b_brown * x

    e_int = math.exp(z_int)
    e_brown = math.exp(z_brown)
    d = 1.0 + e_int + e_brown
    p_blue, p_int, p_brown = 1.0 / d, e_int / d, e_brown / d

    probs = {"blue": p_blue, "intermediate": p_int, "brown": p_brown}
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
    return IrisPlexResult(round(p_blue, 6), round(p_int, 6), round(p_brown, 6), call, conf, counted, notes)


def reference_integrity_ok() -> bool:
    """Biology contract guard: rs12913832 GG -> blue-dominant, AA -> brown-dominant. A corrupted/fabricated
    coefficient set fails this. Returns True iff the model reproduces the known HERC2 pigmentation direction."""
    def _g(gt):
        return {rsid: gt if rsid == "rs12913832" else _absent_hom(allele)
                for rsid, allele, _, _ in IRISPLEX_SNPS}

    def _absent_hom(allele):  # a homozygous genotype with 0 counted alleles
        other = next(b for b in "ACGT" if b != allele)
        return other + other

    gg = predict_eye_color(_g("GG"))   # 0 A at HERC2 (+ 0 counted elsewhere)
    aa = predict_eye_color(_g("AA"))   # 2 A at HERC2
    return gg.call == "blue" and gg.p_blue > 0.7 and aa.call == "brown" and aa.p_brown > 0.7
