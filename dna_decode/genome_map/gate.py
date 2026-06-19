"""Differentiation audit + G1/G2 gate — the make-or-break (Step 6).

Computes the prevent-wrong-inference gate from the genome map ALONE (it uses the
retained `raw_product` / `classification_reason` / `primary_tier` / `phenotype` —
brainstorm catch C4, no re-reading of the raw tools):

  G1 (prevent-wrong-inference): concrete features where the tiered map PREVENTS a
    wrong inference a reader would draw from the raw Bakta TSV —
      (a) demote-homology : a raw product taken as fact (`putative`/`by similarity`/
          DUF…) that the map DEMOTED to `homology-only-hypothesis`.
      (b) surface-determinant : a HIGH-confidence determinant the map SURFACED as
          `determinant-phenotype` that the raw TSV lists as a plain gene.
    "The DB-labelled unknown rate exists" alone does NOT count (the relabelling-
    only / busywork guard) — unknown-tier features are NEVER counted as G1.

  G2 (no-tier-confusion): zero non-`determinant-phenotype` features carry a
    phenotype claim (the phenotype wall, audited from the map).

Per-genome `evaluate_gate` returns a per-genome verdict; `aggregate_spike_verdict`
applies the cross-genome spike rule (G1 needs >=3 on >=1 genome; ANY genome whose
determinant joins are all symbol-fallback -> NO_GO).
"""
from __future__ import annotations

from dna_decode.genome_map import (
    TIER_DETERMINANT_PHENOTYPE,
    TIER_HOMOLOGY_HYPOTHESIS,
)
from dna_decode.genome_map.tier_vocab import LOW_CONFIDENCE_PATTERNS, matches_any

G1_MIN_FEATURES = 3
VERDICT_GO = "GO"
VERDICT_NO_GO = "NO_GO"


def _g1_features(genome_map: dict) -> list[dict]:
    """Extract the prevent-wrong-inference G1 features (type a + type b only).

    Unknown-tier features are intentionally NOT counted (the busywork guard).
    """
    out: list[dict] = []
    for f in genome_map["features"]:
        tier = f["primary_tier"]
        if tier == TIER_HOMOLOGY_HYPOTHESIS and matches_any(f["raw_product"], LOW_CONFIDENCE_PATTERNS):
            # (a) the map demoted a weak-wording product a naive reader would take as fact
            out.append({
                "type": "demote-homology",
                "feature_index": f["feature_index"],
                "raw_product": f["raw_product"],
                "primary_tier": tier,
                "classification_reason": f["classification_reason"],
            })
        elif tier == TIER_DETERMINANT_PHENOTYPE and f["phenotype"]:
            # (b) the map surfaced a determinant the raw TSV lists as a plain gene
            out.append({
                "type": "surface-determinant",
                "feature_index": f["feature_index"],
                "raw_product": f["raw_product"],
                "raw_gene_symbol": f["raw_gene_symbol"],
                "phenotype": f["phenotype"],
            })
    return out


def _g2_violations(genome_map: dict) -> list[dict]:
    """Find non-determinant-phenotype features that carry a phenotype claim (wall breach)."""
    out: list[dict] = []
    for f in genome_map["features"]:
        if f["primary_tier"] != TIER_DETERMINANT_PHENOTYPE and f["phenotype"]:
            out.append({"feature_index": f["feature_index"], "primary_tier": f["primary_tier"]})
    return out


def evaluate_gate(genome_map: dict) -> dict:
    """Per-genome G1/G2 gate result + a per-genome verdict.

    g1_pass = >=3 prevent-wrong-inference features AND >=1 of type (a)/(b)
    (all G1 features are a/b by construction, so the second clause is the
    unknown-rate guard made explicit). The per-genome verdict is NO_GO if the
    genome's determinant joins are ALL symbol-fallback, OR g1 fails, OR g2 fails.
    """
    g1 = _g1_features(genome_map)
    n_demote = sum(1 for x in g1 if x["type"] == "demote-homology")
    n_surface = sum(1 for x in g1 if x["type"] == "surface-determinant")
    g1_pass = len(g1) >= G1_MIN_FEATURES and (n_demote + n_surface) >= 1

    g2_violations = _g2_violations(genome_map)
    g2_pass = len(g2_violations) == 0

    all_symbol_fallback = bool(genome_map["metrics"].get("all_joins_symbol_fallback"))

    verdict = VERDICT_GO
    if all_symbol_fallback or not g1_pass or not g2_pass:
        verdict = VERDICT_NO_GO

    return {
        "genome_accession": genome_map.get("genome_accession"),
        "g1_features": g1,
        "g1_demote_count": n_demote,
        "g1_surface_count": n_surface,
        "g1_pass": g1_pass,
        "g2_spotcheck": {"pass": g2_pass, "violations": g2_violations},
        "all_joins_symbol_fallback": all_symbol_fallback,
        "unknown_under_bakta_db_light": genome_map["metrics"].get("unknown_under_bakta_db_light"),
        "verdict": verdict,
    }


def aggregate_spike_verdict(gate_results: list[dict]) -> dict:
    """Apply the cross-genome spike rule over per-genome gate results.

    Spike GO iff:
      - G1: >=3 prevent-wrong-inference features on >=1 genome (any genome g1_pass), AND
      - G2: no phenotype-wall breach on ANY genome, AND
      - no genome whose determinant joins are ALL symbol-fallback (the trap guard).
    Any failing condition -> NO_GO. (An honest NO-GO is a valid spike outcome.)
    """
    spike_g1_pass = any(g["g1_pass"] for g in gate_results)
    spike_g2_pass = all(g["g2_spotcheck"]["pass"] for g in gate_results)
    any_all_symbol_fallback = any(g["all_joins_symbol_fallback"] for g in gate_results)

    verdict = VERDICT_GO
    reasons: list[str] = []
    if any_all_symbol_fallback:
        verdict = VERDICT_NO_GO
        reasons.append("a genome's determinant joins are ALL symbol-fallback (gene-symbol-trap guard)")
    if not spike_g1_pass:
        verdict = VERDICT_NO_GO
        reasons.append("no genome has >=3 prevent-wrong-inference features (G1 fail / relabelling-only)")
    if not spike_g2_pass:
        verdict = VERDICT_NO_GO
        reasons.append("phenotype-wall breach: a non-determinant feature carried a phenotype (G2 fail)")

    return {
        "verdict": verdict,
        "spike_g1_pass": spike_g1_pass,
        "spike_g2_pass": spike_g2_pass,
        "any_all_symbol_fallback": any_all_symbol_fallback,
        "reasons": reasons or ["G1 satisfied on >=1 genome; G2 clean; no all-symbol-fallback genome"],
    }
