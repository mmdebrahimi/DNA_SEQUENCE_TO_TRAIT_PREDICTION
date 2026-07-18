"""Forward edit router — the capstone that unifies the two working regimes behind one entry point.

Given an edit (gene + amino-acid mutation, optionally a drug), classify which G2P regime it falls in and
route to the RIGHT predictor (per feedback_g2p_decoder_regime_boundary):

  Regime A (edit hits a curated AMR determinant)  -> the determinant catalogue -> R/S call
  Regime B (edit changes protein molecular fitness) -> the DMS-validated predictor (BLOSUM/ESM2)
  Regime C (organism-level polygenic trait)         -> ABSTAIN (closed negative)
  unknown  (no predictor available)                 -> ABSTAIN (honest "no predictor")

The router NEVER sends a resistance edit to the likelihood predictor (which fails on the antagonistic
direction — the resistance-conservativeness finding) nor an organismal trait to any predictor. The
Regime-A membership test is catalogue-driven (inject the resistance-determinant key set; the demo builds it
from the real WHO TB catalogue), so the module stays offline-testable.
"""
from __future__ import annotations

from .variant_effect import predict_effect

REGIME_A = "A_determinant"
REGIME_B = "B_molecular"
REGIME_C = "C_organismal"
REGIME_UNKNOWN = "unknown"


def classify_edit(gene: str, *, drug: str | None = None, determinant_locus: bool = False,
                  molecular_predictor: bool = False, organismal: bool = False) -> tuple[str, str]:
    """Classify an edit into a regime + a one-line reason. Precedence: organismal > determinant > molecular.

    - organismal=True                      -> Regime C (closed negative; abstain)
    - determinant_locus=True AND drug      -> Regime A (a curated AMR determinant locus for `drug`)
    - molecular_predictor=True             -> Regime B (a protein with a DMS/ESM predictor)
    - otherwise                            -> unknown (no predictor)
    """
    if organismal:
        return REGIME_C, "organism-level polygenic trait — closed-negative regime"
    if determinant_locus and drug:
        return REGIME_A, f"{gene} is a curated AMR determinant locus for {drug}"
    if molecular_predictor:
        return REGIME_B, f"{gene} has a DMS/ESM molecular-fitness predictor"
    return REGIME_UNKNOWN, f"no forward predictor available for {gene} (edit not routable)"


def variant_key(gene: str, mutation: str) -> str:
    """Normalized determinant key, e.g. ('rpoB','S450L') -> 'rpoB:S450L' (case-normalized)."""
    return f"{gene}:{mutation.strip().upper()}"


def catalogue_call(gene: str, mutation: str, resistance_keys: set[str]) -> tuple[str, str]:
    """Regime-A forward call: R iff (gene, mutation) is a catalogued resistance determinant, else S.
    S is 'susceptible-by-absence' (callability unassessed) — the same honest semantics as the shipped
    determinant decoder."""
    if variant_key(gene, mutation) in resistance_keys:
        return "R", "matches a catalogued resistance determinant"
    return "S", "not a catalogued resistance determinant (susceptible-by-absence; callability unassessed)"


def predict_edit(gene: str, mutation: str, *, regime: str | None = None, drug: str | None = None,
                 protein_seq: str | None = None, method: str = "blosum62", esm_table: dict | None = None,
                 am_table: dict | None = None, prosst_table: dict | None = None,
                 resistance_keys: set[str] | None = None,
                 phenotype: str = "molecular fitness (DMS-measured)",
                 determinant_locus: bool = False, molecular_predictor: bool = False,
                 organismal: bool = False) -> dict:
    """Unified forward edit -> phenotype. If `regime` is not given, auto-classify via classify_edit.

    Returns a regime-tagged dict; Regime A -> R/S + catalogue; Regime B -> the ForwardPrediction; Regime C /
    unknown -> ABSTAIN. NEVER routes a resistance edit to the likelihood predictor or a trait to a predictor.
    """
    reason = None
    if regime is None:
        regime, reason = classify_edit(gene, drug=drug, determinant_locus=determinant_locus,
                                       molecular_predictor=molecular_predictor, organismal=organismal)

    base = {"gene": gene, "mutation": mutation, "drug": drug, "regime": regime,
            "routing_reason": reason}

    if regime == REGIME_A:
        keys = resistance_keys or set()
        call, why = catalogue_call(gene, mutation, keys)
        return {**base, "prediction": call, "predictor": "determinant_catalogue",
                "confidence": "high" if call == "R" else "medium", "abstain": False, "notes": [why]}

    if regime == REGIME_B:
        # prefer the VALIDATED ESM2+ProSST hybrid when both per-protein tables are supplied
        # (wiki/prosst_lift_2026-07-18.md: +0.067 vs ESM2 alone); else fall through to the single method.
        if esm_table is not None and prosst_table is not None:
            from .variant_effect import predict_variant_hybrid
            pred = predict_variant_hybrid(protein_seq or "", mutation, esm_table=esm_table,
                                          prosst_table=prosst_table, protein=gene, phenotype_axis=phenotype)
            predictor = "hybrid_esm2_prosst_DMS_validated"
        else:
            pred = predict_effect(protein_seq or "", mutation, protein=gene, phenotype_axis=phenotype,
                                  method=method, esm_table=esm_table, am_table=am_table)
            predictor = f"{method}_DMS_validated"
        return {**base, "prediction": pred.predicted_effect, "raw_score": pred.raw_score,
                "predictor": predictor, "confidence": pred.confidence,
                "abstain": False, "notes": pred.notes}

    # Regime C or unknown -> abstain (honest)
    why = ("organism-level polygenic trait — closed-negative regime; abstaining"
           if regime == REGIME_C else "no forward predictor available for this edit; abstaining")
    return {**base, "prediction": "ABSTAIN", "predictor": "none", "confidence": "low",
            "abstain": True, "notes": [why]}
