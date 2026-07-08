"""Integrated multi-gene PGx interpretation — combines per-gene calls into drug-context annotations.

This is the PAYOFF of completing three CPIC drug-gene groups: warfarin (VKORC1 + CYP2C9 + CYP4F2),
statins (SLCO1B1 + ABCG2), thiopurines (TPMT + NUDT15). Each function here takes the already-computed
per-gene phenotype/function strings and returns a QUALITATIVE, annotation-only drug-context summary.

HONESTY (load-bearing): these are DIRECTIONAL annotations, NOT clinical doses. A real warfarin dose needs
the IWPC/Gage algorithm + clinical inputs (age, weight, BSA, amiodarone, indication) that a VCF does not
carry — so we NEVER emit a mg/day number, only a genotype-derived dose-requirement DIRECTION. The thiopurine
combination IS a real CPIC rule (the more-deficient of TPMT/NUDT15 governs), so that one mirrors CPIC's own
two-gene logic. NOT a clinical tool.
"""
from __future__ import annotations

# Canonical phenotype/function strings emitted by the per-gene callers.
_METAB = {"Normal Metabolizer": 0, "Intermediate Metabolizer": 1, "Poor Metabolizer": 2}
_SENS = {"Normal sensitivity": 0, "Intermediate sensitivity": 1, "High sensitivity": 2}
_C4F2 = {"Normal Function": 0, "Intermediate": 1, "Reduced Function": 2}
_TRANSPORT = {"Normal Function": 0, "Decreased Function": 1, "Poor Function": 2}


def _direction(score: float) -> str:
    if score <= -1.5:
        return "much_lower_dose_requirement"
    if score <= -0.5:
        return "lower_dose_requirement"
    if score < 0.5:
        return "standard_dose_requirement"
    if score < 1.5:
        return "higher_dose_requirement"
    return "much_higher_dose_requirement"


def interpret_warfarin(vkorc1_sensitivity: str, cyp2c9_phenotype: str, cyp4f2_function: str) -> dict:
    """Combine the 3 CPIC warfarin genes into a qualitative dose-requirement DIRECTION (annotation only).

    VKORC1 A allele (higher sensitivity) -> LOWER dose; CYP2C9 reduced metabolism -> LOWER dose (slower
    clearance); CYP4F2 *3 reduced function -> HIGHER dose. Additive qualitative model (NOT the IWPC mg/day)."""
    v = _SENS.get(vkorc1_sensitivity)
    c = _METAB.get(cyp2c9_phenotype)
    f = _C4F2.get(cyp4f2_function)
    if v is None or c is None or f is None:
        return {"drug": "warfarin", "status": "indeterminate",
                "reason": "one or more of VKORC1/CYP2C9/CYP4F2 not callable"}
    # VKORC1 sensitivity lowers dose (-), CYP2C9 reduced metabolism lowers dose (-), CYP4F2 *3 raises (+).
    score = -(v) - (c) + (0.5 * f)
    return {
        "drug": "warfarin", "status": "ok",
        "dose_direction": _direction(score),
        "contributors": {
            "VKORC1": f"{vkorc1_sensitivity} (dose {'lower' if v else 'standard'})",
            "CYP2C9": f"{cyp2c9_phenotype} (dose {'lower' if c else 'standard'})",
            "CYP4F2": f"{cyp4f2_function} (dose {'higher' if f else 'standard'})"},
        "caveat": ("Qualitative genotype-derived warfarin dose DIRECTION combining VKORC1 + CYP2C9 + CYP4F2 "
                   "(CPIC Johnson 2017). NOT a mg/day dose — the IWPC/Gage algorithm needs clinical inputs "
                   "(age/weight/BSA/amiodarone) a VCF lacks. Annotation only. NOT a clinical tool."),
    }


def interpret_statins(slco1b1_function: str, abcg2_function: str) -> dict:
    """Combine SLCO1B1 (521T>C, simvastatin) + ABCG2 (Q141K, rosuvastatin) into a statin myopathy-risk band.

    Both are decreased-uptake/efflux -> higher systemic statin exposure -> higher myopathy risk. The two act
    on DIFFERENT statins (SLCO1B1 uptake = simvastatin-dominant; ABCG2 efflux = rosuvastatin), so we report
    per-statin risk, not a single merged number."""
    s = _TRANSPORT.get(slco1b1_function)
    a = _TRANSPORT.get(abcg2_function)
    if s is None or a is None:
        return {"drug": "statins", "status": "indeterminate",
                "reason": "SLCO1B1 or ABCG2 not callable"}
    risk = {0: "typical_risk", 1: "increased_risk", 2: "high_risk"}
    return {
        "drug": "statins", "status": "ok",
        "simvastatin_myopathy_risk": risk[s],     # SLCO1B1-driven
        "rosuvastatin_exposure_risk": risk[a],    # ABCG2-driven
        "contributors": {"SLCO1B1": slco1b1_function, "ABCG2": abcg2_function},
        "caveat": ("Per-statin myopathy risk: SLCO1B1 521T>C drives SIMVASTATIN risk (uptake), ABCG2 Q141K "
                   "drives ROSUVASTATIN exposure (efflux) — different statins, reported separately (CPIC "
                   "Cooper-DeHoff 2022). Annotation only. NOT a clinical tool."),
    }


def interpret_thiopurines(tpmt_phenotype: str, nudt15_phenotype: str) -> dict:
    """Combine TPMT + NUDT15 into thiopurine (azathioprine/mercaptopurine) toxicity risk.

    THIS IS A REAL CPIC RULE (Relling 2019): the MORE-deficient of the two genes governs — a poor
    metabolizer in EITHER gene => high myelosuppression risk (drastically reduce / avoid); an intermediate
    metabolizer in either => reduce dose. We take the max deficiency."""
    t = _METAB.get(tpmt_phenotype)
    n = _METAB.get(nudt15_phenotype)
    if t is None or n is None:
        return {"drug": "thiopurines", "status": "indeterminate",
                "reason": "TPMT or NUDT15 not callable"}
    worst = max(t, n)
    band = {0: ("normal_risk", "standard starting dose"),
            1: ("increased_risk", "start with reduced dose (CPIC IM)"),
            2: ("high_risk", "drastically reduce or select alternative (CPIC PM)")}[worst]
    governing = "TPMT" if t >= n else "NUDT15"
    if t == n:
        governing = "TPMT+NUDT15 (equal)"
    return {
        "drug": "thiopurines", "status": "ok",
        "toxicity_risk": band[0], "dosing_note": band[1],
        "governing_gene": governing,
        "contributors": {"TPMT": tpmt_phenotype, "NUDT15": nudt15_phenotype},
        "caveat": ("Combined thiopurine toxicity risk — CPIC (Relling 2019) uses the MORE-deficient of "
                   "TPMT / NUDT15 (real two-gene rule). Annotation only. NOT a clinical tool."),
    }


def interpret_all(results: dict) -> dict:
    """Build all available drug-context interpretations from a PGP-UK-style results dict.

    `results` maps gene -> {phenotype|function|sensitivity, ...} (the pgx_decode_pgp_uk realizer shape).
    Missing/indeterminate genes yield an indeterminate interpretation for that drug (never a crash)."""
    def _get(gene, key):
        return (results.get(gene) or {}).get(key)
    out = {}
    out["warfarin"] = interpret_warfarin(
        _get("vkorc1", "sensitivity"), _get("cyp2c9", "phenotype"), _get("cyp4f2", "function"))
    out["statins"] = interpret_statins(_get("slco1b1", "function"), _get("abcg2", "function"))
    out["thiopurines"] = interpret_thiopurines(_get("tpmt", "phenotype"), _get("nudt15", "phenotype"))
    return out
