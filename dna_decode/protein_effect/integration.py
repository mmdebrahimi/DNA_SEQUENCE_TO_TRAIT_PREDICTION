"""Rung-2 + rung-3 fused query (2026-07-13).

A single query for a protein point mutation returns BOTH:
  * rung-2 MOLECULAR effect  — the ESM2-650M zero-shot damage rank (from `predictor`), for ANY protein.
  * rung-3 KNOWN PHENOTYPE    — whether a curated AMR catalog documents this mutation as a resistance
                               determinant, routed by organism/gene to the right in-repo catalog.

The routing is honest about coverage:
  * HIV-1 (genes RT / PR / IN / CA)  -> the in-repo `hiv_amr` mutation-level catalog (Stanford HIVDB DRMs):
    reports which drugs it confers resistance to (the deployed decoder's call).
  * M. tuberculosis (rpoB/katG/...)  -> the WHO TB catalogue (`tb_who_catalogue`) grade-1/2 determinants
    (degrades to "catalogue unavailable" when the gitignored master file is absent).
  * anything else (e.g. E. coli)      -> NO in-repo mutation-string catalog: bacterial resistance point
    mutations resolve via AMRFinder on a FULL GENOME (the deterministic decoder), not a single-mutation
    lookup. This is stated, not faked.

The molecular rung is UNIVERSAL (works on any sequence); the phenotype rung is CATALOG-BOUNDED. Keeping the
two ranks separate + honestly labelled is the point — a molecular-effect score is NOT a resistance call, and
a catalog hit is NOT a molecular-effect claim. Frozen decoder surface + hiv_amr catalog untouched (READ-only).
"""
from __future__ import annotations

import json
from pathlib import Path

from . import predictor as P

HIV_GENES = ("RT", "PR", "IN", "CA")


def _hiv_drugs_for_gene(gene: str) -> list[str]:
    import dna_decode.data.hiv_amr as H
    drugs = list(H.HIV_NNRTI_DRUGS) + list(H.NRTI_DRUGS) + list(H.PI_CLASS.drugs) \
        + list(H.INSTI_CLASS.drugs) + list(H.CAI_CLASS.drugs)
    return sorted({d for d in drugs if H.gene_for_hiv_drug(d) == gene})


def _hiv_phenotype(mutation: str, gene: str) -> dict:
    import dna_decode.data.hiv_amr as H
    tok = mutation.upper()
    resistant, determinants = [], set()
    for drug in _hiv_drugs_for_gene(gene):
        call = H.call_hiv_observed(drug, {gene: {tok}})
        if call.prediction == "R":
            resistant.append(drug)
            determinants.update(call.determinants)
    return {
        "catalog": "Stanford HIVDB major-DRM catalog (in-repo dna_decode.data.hiv_amr)",
        "is_known_resistance_mutation": bool(resistant),
        "determinant": sorted(determinants) or None,
        "resistant_drugs": sorted(resistant) or None,
        "caveat": ("v0 catalogs are CLASS-level for NNRTI/position-based for NRTI/PI/INSTI — a known-DRM call "
                   "may over-call per-drug differential resistance (e.g. K103N spares ETR/RPV)."),
    }


def _tb_phenotype(mutation: str, gene: str) -> dict:
    try:
        from dna_decode.data import tb_who_catalogue as W
        W.verify_pins()
        tok = mutation.upper()
        hit_drugs = []
        for drug in ("rifampicin", "isoniazid", "ethambutol", "levofloxacin", "moxifloxacin"):
            try:
                dets = W.load_determinants(drug)
            except Exception:
                continue
            if any(gene.lower() in (d.gene or "").lower() and tok in (d.variant or "").upper() for d in dets):
                hit_drugs.append(drug)
        return {"catalog": "WHO TB mutation catalogue v2 (grade-1/2)", "is_known_resistance_mutation": bool(hit_drugs),
                "determinant": (f"{gene}:{mutation}" if hit_drugs else None),
                "resistant_drugs": sorted(hit_drugs) or None,
                "caveat": "WHO catalogue is knowledge-baseline (built partly from CRyPTIC); grade-1/2 only."}
    except Exception as e:
        return {"catalog": "WHO TB mutation catalogue v2", "is_known_resistance_mutation": None,
                "determinant": None, "resistant_drugs": None,
                "caveat": f"catalogue unavailable ({type(e).__name__}); the 37 MB master file is gitignored."}


def known_phenotype(mutation: str, gene: str, organism: str) -> dict:
    """Route the rung-3 phenotype lookup by organism/gene. Never fabricates a catalog for an unsupported one."""
    org = (organism or "").strip().lower()
    g = (gene or "").strip().upper()
    if org in ("hiv-1", "hiv", "hiv1") and g in HIV_GENES:
        return _hiv_phenotype(mutation, g)
    if org in ("m. tuberculosis", "mtb", "tuberculosis", "mycobacterium tuberculosis"):
        return _tb_phenotype(mutation, gene)
    return {
        "catalog": None, "is_known_resistance_mutation": None, "determinant": None, "resistant_drugs": None,
        "caveat": (f"no in-repo mutation-level AMR catalog for organism={organism!r} gene={gene!r}. Bacterial "
                   "(e.g. E. coli) resistance point mutations resolve via AMRFinder on a FULL GENOME "
                   "(the deterministic decoder) — not a single-mutation-string lookup."),
    }


def load_logp(cache_path: str | Path) -> dict:
    """Load an ESM masked-marginal cache into {pos(int): {aa: log-prob}}, tolerating either the predictor's
    {sequence, logp} wrapper OR a raw {pos: {aa}} map (the hiv_esm_vs_catalog cache shape)."""
    d = json.loads(Path(cache_path).read_text(encoding="utf-8"))
    raw = d.get("logp", d)
    return {int(k): v for k, v in raw.items()}


def integrated_query(mutation: str, gene: str, organism: str, sequence: str, logp: dict) -> dict:
    """Fuse rung-2 (molecular effect from `logp`) + rung-3 (known phenotype). `logp` = {pos: {aa: logp}}."""
    molecular = P.predict(sequence, mutation, logp)
    phenotype = known_phenotype(mutation, gene, organism)
    return {
        "artifact": "integrated_mutation_query", "schema": "integrated-mutation-query-v1",
        "query": {"organism": organism, "gene": gene, "mutation": mutation.upper()},
        "molecular_effect_rung2": {k: molecular[k] for k in
                                   ("damage_llr", "position_percentile", "direction_hint", "sequence_change")},
        "known_phenotype_rung3": phenotype,
        "honest_framing": ("The molecular rung is a UNIVERSAL zero-shot ESM rank score (~0.49-0.52 Spearman "
                           "face-validity on stability DMS); the phenotype rung is a CATALOG lookup (present "
                           "only where a curated mutation-level catalog exists). A molecular-effect score is "
                           "NOT a resistance call and a catalog hit is NOT a molecular-effect claim — they are "
                           "reported as two separate, honestly-labelled ranks."),
        "molecular_caveat": molecular["honest_caveat"],
    }
