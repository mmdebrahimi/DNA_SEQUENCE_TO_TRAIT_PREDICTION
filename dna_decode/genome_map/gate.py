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

    Two DISTINCT signals (do NOT conflate — a single GO scalar would let the
    cheap half stand for the expensive one):

    - `g1_pass` (TIERING honesty): >=3 prevent-wrong-inference features. Every G1
      feature is type (a) demote-homology or (b) surface-determinant by
      construction, so a genome with ONLY demotions (no determinants, e.g. a
      determinant-free homology-heavy genome) legitimately passes — that is the
      tier vocabulary working, NOT the determinant overlay. (The old "AND >=1 of
      type a/b" clause was a tautology — the qualifying set IS only a/b — and is
      removed; it never constrained anything.)
    - `overlay_go` (OVERLAY integrity, scoped to AMR-bearing genomes): this genome
      has determinant rows AND >=1 high-confidence join AND >=1 surfaced
      determinant-phenotype feature AND G2 clean AND not all-symbol-fallback. This
      is the signal that the determinant->feature join actually works; it is
      None (not False) for a determinant-free genome (n_main_rows == 0), which has
      nothing to demonstrate.

    The per-genome `verdict` is the TIERING verdict (NO_GO if all-symbol-fallback,
    g1 fails, or g2 fails) — back-compat; `overlay_go` is reported alongside it.
    """
    g1 = _g1_features(genome_map)
    n_demote = sum(1 for x in g1 if x["type"] == "demote-homology")
    n_surface = sum(1 for x in g1 if x["type"] == "surface-determinant")
    g1_pass = len(g1) >= G1_MIN_FEATURES

    g2_violations = _g2_violations(genome_map)
    g2_pass = len(g2_violations) == 0

    metrics = genome_map["metrics"]
    all_symbol_fallback = bool(metrics.get("all_joins_symbol_fallback"))

    # Overlay-integrity evidence (the expensive half, separate from tiering).
    jq = metrics.get("join_quality", {})
    n_main_rows = jq.get("n_main_rows", 0)
    n_high = jq.get("n_high_confidence_join", 0)
    n_det_pheno = metrics.get("determinant_phenotype_feature_count", 0)
    overlay_evidence = {
        "n_main_rows": n_main_rows,
        "n_high_confidence_join": n_high,
        "determinant_phenotype_feature_count": n_det_pheno,
        "surface_determinant_count": n_surface,
    }
    # None for a determinant-free genome (nothing to demonstrate); else a real bool.
    if n_main_rows == 0:
        overlay_go = None
    else:
        overlay_go = bool(
            n_high > 0 and n_det_pheno > 0 and n_surface > 0
            and g2_pass and not all_symbol_fallback
        )

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
        "unknown_under_bakta_db_light": metrics.get("unknown_under_bakta_db_light"),
        "overlay_evidence": overlay_evidence,
        "overlay_go": overlay_go,
        "verdict": verdict,
    }


def aggregate_spike_verdict(gate_results: list[dict]) -> dict:
    """Apply the cross-genome spike rule — TWO distinct verdicts, not one scalar.

    `tiering_go` (the headline `verdict`, back-compat) — the Bakta honesty
    re-tiering works:
      - G1: >=3 prevent-wrong-inference features on >=1 genome, AND
      - G2: no phenotype-wall breach on ANY genome, AND
      - no genome whose determinant joins are ALL symbol-fallback.

    `overlay_go` (the SEPARATE determinant-overlay-integrity claim) — at least one
    AMR-bearing genome demonstrated a working determinant->feature join (per-genome
    `overlay_go is True`). A `tiering_go` GO carried entirely by demotions (no
    genome with a real overlay) is honest about re-tiering but does NOT license
    "the overlay is validated" — keep the two claims separate so neither word
    overstates the other.
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

    # overlay-integrity: at least one AMR-bearing genome demonstrated a real join.
    overlay_go = any(g.get("overlay_go") is True for g in gate_results)
    n_overlay_genomes = sum(1 for g in gate_results if g.get("overlay_go") is True)
    overlay_reason = (
        f"determinant-overlay integrity demonstrated on {n_overlay_genomes} AMR-bearing genome(s)"
        if overlay_go else
        "determinant-overlay integrity NOT demonstrated (no AMR-bearing genome with a high-confidence "
        "join + surfaced determinant) — tiering GO does not license an overlay claim"
    )

    return {
        "verdict": verdict,
        "tiering_go": verdict == VERDICT_GO,
        "overlay_go": overlay_go,
        "overlay_reason": overlay_reason,
        "n_overlay_genomes": n_overlay_genomes,
        "spike_g1_pass": spike_g1_pass,
        "spike_g2_pass": spike_g2_pass,
        "any_all_symbol_fallback": any_all_symbol_fallback,
        "reasons": reasons or ["G1 satisfied on >=1 genome; G2 clean; no all-symbol-fallback genome"],
    }
