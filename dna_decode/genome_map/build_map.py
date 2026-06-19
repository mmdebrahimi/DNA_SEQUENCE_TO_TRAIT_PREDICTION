"""Map assembler — raw-field JSON/table + DB-labelled unknown rate (Step 5).

Assembles the per-feature genome map from the parsed Bakta features (Step 1) +
the joined determinant hits (Step 4) + per-drug genome-level verdicts. Each
feature gets ONE primary tier:

    determinant-phenotype  — IFF a HIGH-confidence determinant join hit it
                             (symbol-fallback does NOT qualify — the phenotype wall)
    else                   — classify_feature_tier(product, gene_symbol)

The map RETAINS the raw fields (raw_product / raw_gene_symbol / raw_locus_tag /
raw_feature_type / source_tool / classification_reason / secondary_evidence) so
Step 6's G1 audit computes from the map alone (brainstorm catch C4).

PURE: no Docker / no file IO. The caller (Step 7 spike) parses the GFF, runs the
overlay, computes `drug_verdicts` via `amr_rules.call_resistance`, and passes
them in — so this module is fully unit-testable on synthetic inputs.

Phenotype WALL (R1/AC3): the `phenotype` field is non-empty ONLY on
`determinant-phenotype` features. Unknown rate (R2/AC4) is DB-labelled
`unknown_under_bakta_db_light`.
"""
from __future__ import annotations

from collections import defaultdict

from dna_decode.genome_map import (
    PHENOTYPE_TIER,
    TIER_DETERMINANT_PHENOTYPE,
    TIER_PRECEDENCE,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.phenotype_overlay import (
    SYMBOL_FALLBACK,
    determinant_phenotype_field,
)
from dna_decode.genome_map.tiers import classify_feature_tier

SCHEMA_VERSION = "genome-map-v1"


def _determinant_counts_for_drug(hit, drug: str) -> bool:
    """True iff this determinant would be COUNTED toward `drug`'s DEPLOYED R/S rule.

    Mirrors the per-drug inclusion in ``amr_rules.call_resistance`` (via
    ``rule_for``), NOT the broad AMRFinder class match. The per-feature drug label
    MUST use the same refined inclusion as the genome-level rule — otherwise the
    map labels a feature with a drug the deployed caller intentionally excludes
    (e.g. a narrow ``blaTEM-1`` -> ceftriaxone, or a ``qnrB``/``oqxAB`` efflux gene
    -> ciprofloxacin), re-introducing exactly the over-calling the frozen rule was
    refined away from. READ-ONLY against the frozen surface (imports only).

    Inclusion per the deployed config:
      - counter='qrdr_point' (cipro): a QRDR target-alteration POINT mutation
        (gyrA/gyrB/parC/parE with AMRFinder Method POINT*).
      - subclass_any (cef/gent/mero/oxacillin): broad drug-class match AND the
        Subclass contains a refinement token.
      - gene_prefixes (tet): broad class match AND the symbol starts with a prefix.
      - else: the broad drug-class match.
    """
    from dna_decode.data.mic_tiers import amrfinder_classes_for
    from dna_decode.eval.amr_rules import QRDR_GENES, rule_for

    cfg = rule_for(drug)
    sym = hit.symbol or ""
    cls = (hit.cls or "").upper()
    sub = (hit.subclass or "").upper()

    if cfg.get("counter") == "qrdr_point":
        meth = (hit.method or "").upper()
        return "POINT" in meth and any(sym == g or sym.startswith(g + "_") for g in QRDR_GENES)

    classes = {c.upper() for c in amrfinder_classes_for(drug)}
    if not any(c in cls or c in sub for c in classes):
        return False
    refine = cfg.get("subclass_any")
    if refine is not None and not any(t.upper() in sub for t in refine):
        return False
    prefixes = cfg.get("gene_prefixes")
    if prefixes and not sym.lower().startswith(tuple(p.lower() for p in prefixes)):
        return False
    return True


def _phenotype_for_feature(high_hits, drug_verdicts: dict, drugs: list[str]) -> list[dict]:
    """Build the phenotype list for a determinant-phenotype feature (the wall).

    For each HIGH-confidence determinant hit, attach a drug-specific phenotype
    entry ONLY for drugs whose DEPLOYED rule actually counts that determinant
    (`_determinant_counts_for_drug`, mirroring ``call_resistance``); the entry
    carries `drug_rule_counted=True` + the SEPARATE genome-level prediction
    (ABSTAIN-aware) — it reports "this determinant counts toward the drug's rule
    and the genome-level call is X", never "this feature IS resistant". A
    determinant that NO requested drug's refined rule counts still surfaces with
    its mechanism/class (drug=None, DETERMINANT_PRESENT) — it IS a curated
    determinant, but the map asserts no drug for it (the over-call guard).
    """
    out: list[dict] = []
    for jh in high_hits:
        counted = [d for d in drugs if _determinant_counts_for_drug(jh.hit, d)]
        if counted:
            for d in counted:
                field = determinant_phenotype_field(jh, d, drug_verdicts.get(d))
                if field is not None:
                    out.append(field)
        else:
            # Curated determinant present but no requested drug's REFINED rule counts it:
            # show the mechanism/class, assert no drug (do NOT over-call via broad class).
            out.append({
                "drug": None,
                "determinant_symbol": jh.hit.symbol,
                "amrfinder_class": jh.hit.cls,
                "amrfinder_subclass": jh.hit.subclass,
                "method": jh.hit.method,
                "join_confidence": jh.join_confidence,
                "drug_rule_counted": False,
                "phenotype": "DETERMINANT_PRESENT",
                "provenance": "amrfinder_curated_determinant",
                "abstain": False,
            })
    return out


def build_genome_map(
    accession: str,
    organism: str | None,
    features,
    joined_hits: list,
    join_counts: dict,
    *,
    drug_verdicts: dict | None = None,
    drugs: list[str] | None = None,
    degraded: bool = False,
) -> dict:
    """Assemble the per-feature genome map dict.

    `joined_hits` is the output of phenotype_overlay.join_hits (its feature_index
    points into `features`' row order). `drug_verdicts` maps drug ->
    call_resistance verdict (for the per-feature ABSTAIN + the genome-level
    calls block). `degraded` flags the offline path (AC12).
    """
    drug_verdicts = drug_verdicts or {}
    drugs = drugs or sorted(drug_verdicts.keys())

    # Group joined hits by feature index.
    by_feature: dict[int, list] = defaultdict(list)
    for jh in joined_hits:
        if jh.feature_index is not None:
            by_feature[jh.feature_index].append(jh)

    rows = list(features.to_dict("records"))
    feature_maps: list[dict] = []
    per_tier = {t: 0 for t in TIER_PRECEDENCE}

    for i, r in enumerate(rows):
        product = str(r.get("product") or "")
        gene_symbol = str(r.get("gene_symbol") or "")
        hits_here = by_feature.get(i, [])
        high_hits = [jh for jh in hits_here if jh.is_high_confidence]
        fallback_hits = [jh for jh in hits_here if jh.join_confidence == SYMBOL_FALLBACK]

        func_tier, func_reason = classify_feature_tier(product, gene_symbol)
        secondary: list[dict] = []
        phenotype: list[dict] = []

        if high_hits:
            primary_tier = TIER_DETERMINANT_PHENOTYPE
            classification_reason = (
                f"high-confidence determinant join "
                f"({', '.join(sorted({jh.join_confidence for jh in high_hits}))})"
            )
            phenotype = _phenotype_for_feature(high_hits, drug_verdicts, drugs)
            # retain the molecular-function read as secondary evidence
            secondary.append({"type": "molecular_function", "tier": func_tier, "reason": func_reason})
        else:
            primary_tier = func_tier
            classification_reason = func_reason

        # Symbol-fallback determinants are VISIBLE secondary evidence — never phenotype.
        for jh in fallback_hits:
            secondary.append({
                "type": "determinant_symbol_fallback",
                "symbol": jh.hit.symbol,
                "amrfinder_class": jh.hit.cls,
                "join_confidence": SYMBOL_FALLBACK,
                "note": "gene-symbol-only match (the 0%-overlap trap) — NOT a phenotype claim",
            })

        per_tier[primary_tier] += 1
        source_tool = "amrfinder+bakta" if (high_hits or fallback_hits) else "bakta"

        feature_maps.append({
            "feature_index": i,
            "seqid": str(r.get("seqid") or ""),
            "start": int(r.get("start") or 0),
            "end": int(r.get("end") or 0),
            "strand": str(r.get("strand") or ""),
            "primary_tier": primary_tier,
            "classification_reason": classification_reason,
            "raw_product": product,
            "raw_gene_symbol": gene_symbol,
            "raw_locus_tag": str(r.get("locus_tag") or ""),
            "raw_feature_type": str(r.get("type") or ""),
            "source_tool": source_tool,
            "secondary_evidence": secondary,
            # phenotype is non-empty ONLY on determinant-phenotype features (the wall).
            "phenotype": phenotype,
        })

    total = len(feature_maps)
    determinant_features = [f for f in feature_maps if f["primary_tier"] == PHENOTYPE_TIER]

    # genome-level R/S calls, tier-tagged + reported SEPARATELY from the features (AC7).
    genome_level_calls = {
        d: {
            "prediction": v.get("prediction"),
            "confidence": v.get("confidence"),
            "n_determinants": v.get("n_determinants"),
            "rule": v.get("rule"),
        }
        for d, v in drug_verdicts.items()
    }

    metrics = {
        "total_features": total,
        "per_tier_counts": per_tier,
        # R2/AC4: the unknown rate carries the DB/version coverage caveat IN its name.
        "unknown_under_bakta_db_light": (per_tier[TIER_UNKNOWN] / total) if total else 0.0,
        "determinant_phenotype_feature_count": len(determinant_features),
        "determinant_phenotype_features": [
            {
                "feature_index": f["feature_index"],
                "seqid": f["seqid"], "start": f["start"], "end": f["end"],
                "raw_product": f["raw_product"], "raw_gene_symbol": f["raw_gene_symbol"],
                "phenotype": f["phenotype"],
            }
            for f in determinant_features
        ],
        "join_quality": dict(join_counts),
        "all_joins_symbol_fallback": (
            join_counts.get("n_main_rows", 0) > 0
            and join_counts.get("n_high_confidence_join", 0) == 0
        ),
        "genome_level_calls": genome_level_calls,
    }

    return {
        "artifact": "genome_map",
        "schema_version": SCHEMA_VERSION,
        "genome_accession": accession,
        "amrfinder_organism": organism,
        "degraded_coverage": degraded,
        "features": feature_maps,
        "metrics": metrics,
    }


def build_feature_table(genome_map: dict) -> list[dict]:
    """Flatten the genome map to a flat per-feature table (one row per feature).

    Phenotype is collapsed to a compact string (the matched drugs / property) for
    the table view; the full structured phenotype stays in the JSON map.
    """
    table: list[dict] = []
    for f in genome_map["features"]:
        pheno = f["phenotype"]
        if pheno:
            drugs = sorted({(p.get("drug") or p.get("phenotype") or "") for p in pheno if p})
            pheno_str = "; ".join(s for s in drugs if s)
        else:
            pheno_str = ""
        table.append({
            "feature_index": f["feature_index"],
            "seqid": f["seqid"],
            "start": f["start"],
            "end": f["end"],
            "strand": f["strand"],
            "primary_tier": f["primary_tier"],
            "raw_gene_symbol": f["raw_gene_symbol"],
            "raw_locus_tag": f["raw_locus_tag"],
            "raw_product": f["raw_product"],
            "classification_reason": f["classification_reason"],
            "phenotype": pheno_str,
            "source_tool": f["source_tool"],
        })
    return table
