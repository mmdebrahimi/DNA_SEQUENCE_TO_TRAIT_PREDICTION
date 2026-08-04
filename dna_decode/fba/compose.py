"""Compose `forward` + `fba`: a POINT MUTATION -> a downstream cell-level trait.

The FBA cell alone takes a full gene KNOCKOUT. Real genotype edits are usually POINT MUTATIONS. This
module chains the two shipped decoders to close that gap:

    missense edit  --forward-->  does it break the enzyme? (LOF)  --fba-->  cell-level trait

- `forward` (DMS-validated, Regime B) scores the missense and returns its OWN method-aware call
  `predicted_effect in {preserved, damaging, uncertain, abstain}` (BLOSUM62 -2 / ESM2 -5 / hybrid 0.33
  thresholds — we INHERIT them, we do not invent a new one).
- damaging  -> model the gene as knocked out -> FBA growth/essentiality (the cell-level consequence).
- preserved -> the variant is tolerated -> no metabolic change -> wild-type growth.
- uncertain / abstain -> DO NOT force a call. Report the CONDITIONAL both-ways (if-LOF vs if-tolerated),
  because forward is a validated RANKER, not a calibrated LOF probability — a forced binary here would be
  the theater the honest rails guard against.

HONEST inheritance: the chain's upstream link (missense->LOF) carries forward's DMS validation; the
downstream link (LOF->cell trait) carries fba's Keio validation (accuracy 0.954). What is NEW and
UNVALIDATED is the binarization threshold on a ranker -> the `lof_call` is labelled a heuristic, and the
`uncertain` pass-through is what keeps it honest.
"""
from __future__ import annotations

# forward.predicted_effect -> what we do in FBA
_LOF_EFFECTS = {"damaging"}
_TOLERATED_EFFECTS = {"preserved"}
_UNCERTAIN_EFFECTS = {"uncertain", "abstain"}


def decide_fba_action(predicted_effect: str) -> str:
    """PURE: map forward's predicted_effect -> an FBA action.

    -> "knockout"   (damaging: model the enzyme as lost)
    -> "wildtype"   (preserved: variant tolerated, no metabolic change)
    -> "conditional" (uncertain/abstain: report BOTH ways, force nothing)
    """
    if predicted_effect in _LOF_EFFECTS:
        return "knockout"
    if predicted_effect in _TOLERATED_EFFECTS:
        return "wildtype"
    if predicted_effect in _UNCERTAIN_EFFECTS:
        return "conditional"
    return "conditional"  # unknown label -> stay honest, don't force


def variant_to_cell_trait(
    model,
    gene_id: str,
    protein_seq: str,
    mutation: str,
    *,
    method: str = "blosum62",
    frac: float = 0.01,
    **forward_kwargs,
) -> dict:
    """missense (gene_id + protein_seq + mutation) -> cell-level trait via forward -> fba.

    `model` is a loaded iML1515 cobra model; `gene_id` is its model gene id (resolved by the caller).
    Returns a `fba-variant-trait-v1` record. Lazy-imports forward + fba engine.
    """
    from dna_decode.forward.variant_effect import predict_effect

    from .model import call_essential, knockout_growth, wildtype_growth

    fwd = predict_effect(
        protein_seq=protein_seq, mutation=mutation, method=method, **forward_kwargs
    ).as_dict()
    action = decide_fba_action(fwd["predicted_effect"])

    wt = wildtype_growth(model)
    ko = knockout_growth(model, gene_id)
    ko_essential = call_essential(ko, wt, frac)

    rec = {
        "record": "fba-variant-trait-v1",
        "gene": gene_id,
        "mutation": mutation,
        "wildtype_growth_per_h": round(wt, 4),
        # forward (upstream link, DMS-validated ranker)
        "forward": {
            "method": fwd["method"],
            "raw_score": fwd["raw_score"],
            "predicted_effect": fwd["predicted_effect"],
            "confidence": fwd["confidence"],
        },
        "lof_call": {
            "damaging": "LOF",
            "preserved": "TOLERATED",
            "uncertain": "UNCERTAIN",
            "abstain": "UNCERTAIN",
        }.get(fwd["predicted_effect"], "UNCERTAIN"),
        "fba_action": action,
    }

    if action == "knockout":
        rec["ko_growth_per_h"] = round(ko, 4)
        rec["growth_fraction_of_wt"] = round(ko / wt, 4) if wt else 0.0
        rec["cell_trait"] = "NON-VIABLE (essential gene lost)" if ko_essential else "viable, altered flux"
        rec["essential_if_lost"] = ko_essential
    elif action == "wildtype":
        rec["ko_growth_per_h"] = None
        rec["cell_trait"] = "viable (variant tolerated -> no metabolic change)"
        rec["essential_if_lost"] = ko_essential  # informational: what a KO WOULD do
    else:  # conditional
        rec["cell_trait_if_LOF"] = "NON-VIABLE (essential gene lost)" if ko_essential else "viable, altered flux"
        rec["cell_trait_if_tolerated"] = "viable (no metabolic change)"
        rec["ko_growth_per_h_if_LOF"] = round(ko, 4)
        rec["essential_if_lost"] = ko_essential

    rec["honesty"] = (
        "chain inherits forward's DMS validation (missense->LOF, a validated RANK) + fba's Keio "
        "validation (LOF->cell trait, accuracy 0.954). The LOF binarization is forward's own method-aware "
        "threshold (heuristic, NOT a calibrated LOF probability); uncertain -> reported both-ways, not forced."
    )
    rec["scope"] = "METABOLIC traits only; missense in an iML1515 metabolic gene. NOT clinical."
    return rec
