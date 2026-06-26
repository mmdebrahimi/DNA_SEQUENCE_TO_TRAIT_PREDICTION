"""Independent functional-evidence layer for the PGx cells (the circularity-break).

The PGx phenotype is "faithful-to-CPIC", and the per-allele FUNCTION assignments ARE CPIC's -- the one
circular link (caller + "truth" both from CPIC). This module attaches, per catalogued allele, an
INDEPENDENT (non-CPIC) functional signal and a verdict (AGREE / DISAGREE / FLAG / NO_SIGNAL) vs the CPIC
function. Agreement raises confidence; DISAGREE/FLAG surface where "faithful-to-CPIC" rests on clinical
evidence the sequence-level signals miss.

HONESTY (load-bearing): the independent signals are ORTHOGONAL to CPIC's curation, NOT ground truth:
  * missense -> Ensembl VEP variant-effect predictors (SIFT/PolyPhen) -- ML predictions, independent of CPIC.
  * stop_gained / splice -> the sequence CONSEQUENCE class itself (a fact, not a model).
  * regulatory promoter -> the documented cis-regulatory EXPRESSION effect from the primary functional
    literature (a MEASURED signal, separate from CPIC's function call). GTEx-liver-significant-eQTL was the
    intended confirmation but did not resolve via the GTEx v2 API in this run (true absence or threshold) ->
    recorded as a deferred enhancement, NOT asserted.
Values curated 2026-06-25 from a live Ensembl VEP REST fetch (rs1799853/rs1057910/rs4244285/rs4986893) +
cited primary literature. This is a small-N per-allele EVIDENCE ANNOTATION, not a concordance-% claim.
"""
from __future__ import annotations

from dataclasses import dataclass

SCHEMA = "pgx-functional-evidence-v0"


@dataclass(frozen=True)
class IndependentEvidence:
    gene: str
    allele: str
    rsid: str
    cpic_function: str        # CPIC's assignment (the thing being independently cross-checked)
    variant_class: str        # missense | stop_gained | synonymous_cryptic_splice | regulatory_promoter
    independent_signal: str   # human-readable independent signal
    signal_source: str        # VEP / consequence / primary-literature
    verdict: str              # AGREE | DISAGREE | FLAG | NO_SIGNAL
    note: str


def derive_missense_verdict(cpic_reduces_function: bool, predictor_damaging: bool | None) -> str:
    """Verdict for a MISSENSE allele: does the (independent) predictor agree with CPIC's direction?
    cpic_reduces_function = CPIC says decreased/no function. predictor_damaging = predictor leans damaging
    (None = no usable prediction)."""
    if predictor_damaging is None:
        return "NO_SIGNAL"
    if predictor_damaging == cpic_reduces_function:
        return "AGREE"
    return "DISAGREE"


# Per-allele independent evidence (curated from the grounded fetch + literature; see module docstring).
EVIDENCE: list[IndependentEvidence] = [
    IndependentEvidence(
        "CYP2C19", "*2", "rs4244285", "no function", "synonymous_cryptic_splice",
        "Ensembl VEP most-severe = synonymous_variant; the no-function mechanism is a cryptic splice site "
        "(aberrant splicing, de Morais 1994) NOT captured by the consequence class",
        "VEP + primary-literature", "FLAG",
        "A consequence-only predictor would UNDER-call this (synonymous); CPIC no-function rests on the "
        "documented splice defect -- a real flag, not a defect of the cell."),
    IndependentEvidence(
        "CYP2C19", "*3", "rs4986893", "no function", "stop_gained",
        "Ensembl VEP most-severe = stop_gained (p.W212X) -> unambiguous loss of function",
        "VEP-consequence", "AGREE",
        "Premature-stop consequence independently confirms CPIC no-function."),
    IndependentEvidence(
        "CYP2C19", "*17", "rs12248560", "increased function", "regulatory_promoter",
        "Documented cis-regulatory promoter variant INCREASING CYP2C19 expression via a new GATA site "
        "(Sim 2006); GTEx-liver-significant-eQTL not resolved via the v2 API this run",
        "primary-literature (Sim 2006)", "AGREE",
        "Independent expression-INCREASE direction matches CPIC increased-function; GTEx-eQTL confirmation "
        "deferred (not asserted)."),
    IndependentEvidence(
        "CYP2C9", "*2", "rs1799853", "decreased function", "missense",
        "Ensembl VEP: PolyPhen probably/possibly_damaging, SIFT deleterious(low-conf)/tolerated mixed "
        "(p.R144C) -> damaging-leaning",
        "VEP-predictor", "AGREE",
        "Damaging-leaning predictors match CPIC decreased-function."),
    IndependentEvidence(
        "CYP2C9", "*3", "rs1057910", "no function", "missense",
        "Ensembl VEP: PolyPhen BENIGN, SIFT mixed (p.I359L) -> predictors UNDER-call",
        "VEP-predictor", "DISAGREE",
        "In-silico predictors call this conservative I->L substitution benign, but CPIC assigns NO function "
        "on clinical/functional evidence -- a case where faithful-to-CPIC rests on MORE than sequence "
        "prediction. The honest value of this layer."),
    IndependentEvidence(
        "VKORC1", "-1639A", "rs9923231", "decreased expression (warfarin sensitivity)", "regulatory_promoter",
        "Documented cis-regulatory promoter variant LOWERING VKORC1 expression (Rieder 2005); "
        "GTEx-liver-significant-eQTL not resolved via the v2 API this run",
        "primary-literature (Rieder 2005)", "AGREE",
        "Independent expression-DECREASE direction matches the warfarin-sensitivity (low-dose) call."),
]


def summary() -> dict:
    from collections import Counter
    c = Counter(e.verdict for e in EVIDENCE)
    return {"n": len(EVIDENCE), "AGREE": c["AGREE"], "DISAGREE": c["DISAGREE"],
            "FLAG": c["FLAG"], "NO_SIGNAL": c["NO_SIGNAL"]}
