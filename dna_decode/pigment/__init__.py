"""Visible-trait pigmentation decoder — the deterministic curated-catalog form of "DNA -> appearance".

v0 = IrisPlex eye colour (6 SNPs -> P(blue)/P(intermediate)/P(brown)) via a published multinomial-logistic
coefficient table. This is the project's WINNING paradigm (a curated causal-locus catalog + published
coefficients applied deterministically) pointed at a NEW trait class (visible traits), NOT a learned genomic
embedding (the 0-for-5 closed regime). Scope: benign visible-trait genetics (textbook eye-colour prediction) —
explicitly NOT a forensic/surveillance tool.
"""
from dna_decode.pigment.irisplex import (
    IRISPLEX_SNPS,
    IrisPlexResult,
    MissingGenotypeError,
    predict_eye_color,
    reference_integrity_ok,
)

__all__ = [
    "predict_eye_color",
    "IrisPlexResult",
    "IRISPLEX_SNPS",
    "MissingGenotypeError",
    "reference_integrity_ok",
]
