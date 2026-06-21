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
    TIER_VIRULENCE_DETERMINANT,
)
from dna_decode.genome_map.phenotype_overlay import (
    SYMBOL_FALLBACK,
    determinant_phenotype_field,
)
from dna_decode.genome_map.tiers import classify_feature_tier
from dna_decode.genome_map.virulence_overlay import cluster_pathotype_context

SCHEMA_VERSION = "genome-map-v1"


def _counted_symbols_by_drug(drug_verdicts: dict, drugs: list[str]) -> dict[str, set[str]]:
    """Per-drug set of determinant SYMBOLS the DEPLOYED ``call_resistance`` counted.

    The single source of truth for "does this determinant count toward drug d":
    read it from the verdict the deployed caller already produced, rather than
    RE-deriving inclusion in this module. The verdict's ``determinants`` list is
    exactly the determinants that passed the deployed rule — DEFAULT or the
    organism-CALIBRATED registry branch (which `run_genome_map_for` activates via
    ``call_resistance(organism=...)``), and ``[]`` on ABSTAIN/EXPRESSION_FLOOR/
    INDETERMINATE. Matching against this set mirrors the deployed inclusion by
    construction — no divergence possible (calibrated, default, or abstain).

    Residual (accepted for v1): symbol-only matching marks ALL same-symbol
    high-confidence rows as counted; fine for current rule shapes (duplicate
    tandem copies all correctly count; QRDR symbols carry the mutation suffix).
    """
    out: dict[str, set[str]] = {}
    for d in drugs:
        dets = (drug_verdicts.get(d) or {}).get("determinants") or []
        out[d] = {(det.get("symbol") or "").strip() for det in dets if det.get("symbol")}
    return out


def _phenotype_for_feature(high_hits, drug_verdicts: dict, drugs: list[str]) -> list[dict]:
    """Build the phenotype list for a determinant-phenotype feature (the wall).

    Three honest per-hit branches (first match wins, per drug):
      1. COUNTED — the determinant is in the drug's deployed verdict `determinants`
         (mirrors ``call_resistance`` exactly, calibrated or default) -> a drug
         entry with `drug_rule_counted=True` + the SEPARATE genome-level prediction
         + `threshold_met` (never "this feature IS resistant").
      2. ABSTAIN-relevant (AC8) — the drug's deployed verdict is ABSTAIN/SUSPEND
         (e.g. an EXPRESSION_FLOOR organism×drug returns `determinants=[]`) AND the
         determinant is broadly relevant to the drug's class -> an explicit ABSTAIN
         entry (`drug_rule_counted=False`, NOT a forced R call). The genome-level
         overlay propagates the abstain in `genome_level_calls`.
      3. DETERMINANT_PRESENT — a curated determinant no requested drug counted and
         none abstain on -> mechanism/class only, drug=None (the over-call guard).
    """
    counted_symbols = _counted_symbols_by_drug(drug_verdicts, drugs)
    out: list[dict] = []
    for jh in high_hits:
        sym = (jh.hit.symbol or "").strip()
        counted = [d for d in drugs if sym in counted_symbols[d]]
        abstain_drugs = [
            d for d in drugs
            if d not in counted
            and (drug_verdicts.get(d) or {}).get("prediction") in {"ABSTAIN", "INDETERMINATE", "SUSPEND"}
            and _hit_broad_class_matches(jh.hit, d)
        ]
        if counted:
            for d in counted:
                field = determinant_phenotype_field(jh, d, drug_verdicts.get(d))
                if field is not None:
                    out.append(field)
        elif abstain_drugs:
            for d in abstain_drugs:
                out.append({
                    "drug": d,
                    "determinant_symbol": jh.hit.symbol,
                    "amrfinder_class": jh.hit.cls,
                    "amrfinder_subclass": jh.hit.subclass,
                    "method": jh.hit.method,
                    "join_confidence": jh.join_confidence,
                    "drug_rule_counted": False,
                    "threshold_met": False,
                    "phenotype": "ABSTAIN",
                    "genome_prediction": "ABSTAIN",
                    "provenance": (drug_verdicts.get(d) or {}).get("rule", "amrfinder_curated_determinant"),
                    "abstain": True,
                })
        else:
            # Curated determinant present but no requested drug counted it and none abstain:
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


def _hit_broad_class_matches(hit, drug: str) -> bool:
    """True iff the determinant's Class/Subclass is broadly relevant to `drug`.

    Used ONLY to decide whether an ABSTAINing drug's overlay should surface an
    ABSTAIN annotation on this determinant (AC8 relevance) — NEVER to make a
    resistance call (that is verdict-derived). An over-broad match here only adds
    an honest ABSTAIN note, never an R claim.
    """
    from dna_decode.data.mic_tiers import amrfinder_classes_for

    try:
        classes = {c.upper() for c in amrfinder_classes_for(drug)}
    except Exception:  # noqa: BLE001 — unknown drug -> not relevant (CLI validates --drugs)
        return False
    cls = (hit.cls or "").upper()
    sub = (hit.subclass or "").upper()
    return any(c in cls or c in sub for c in classes)


def _virulence_field(jh, db_sha: str | None) -> dict:
    """The per-feature virulence wall entry for ONE high-confidence VF join.

    Presence of a curated VF determinant — vf_gene + cluster + the pathotype(s) the
    cluster contributes to (clustered only) + the non-independence caveat + DB SHA.
    NEVER a learned pathogenicity claim (the virulence analog of the AMR phenotype wall).
    """
    from dna_decode.pathotype.vf_runner import NON_INDEPENDENCE_CAVEAT

    cluster = jh.hit.subclass or None
    return {
        "vf_gene": jh.hit.symbol,
        "allele_id": jh.hit.name,
        "cluster": cluster,
        "pathotype_context": cluster_pathotype_context(cluster),
        "join_confidence": jh.join_confidence,
        "caveat": NON_INDEPENDENCE_CAVEAT,
        "db_sha": db_sha,
        "claim": "curated_virulence_determinant_present",  # presence, never "pathogenic"
    }


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
    virulence_joined_hits: list | None = None,
    virulence_join_counts: dict | None = None,
    pathotype_call: dict | None = None,
    virulence_db_sha: str | None = None,
) -> dict:
    """Assemble the per-feature genome map dict.

    `joined_hits` is the output of phenotype_overlay.join_hits (its feature_index
    points into `features`' row order). `drug_verdicts` maps drug ->
    call_resistance verdict (for the per-feature ABSTAIN + the genome-level
    calls block). `degraded` flags the offline path (AC12).

    Virulence overlay (v2, optional — all default None so the AMR-only path is
    byte-unchanged): `virulence_joined_hits` (from virulence_overlay.join_virulence),
    `virulence_join_counts` (VF-namespaced quality counts), `pathotype_call` (the
    genome-level overlay), `virulence_db_sha`. Per-feature precedence is AMR
    determinant-phenotype > virulence-determinant > function tier (Open Question B).
    The virulence metrics are SEPARATE keys (M2) — they never feed the AMR
    `all_joins_symbol_fallback` nor the AMR spike gate.
    """
    drug_verdicts = drug_verdicts or {}
    drugs = drugs or sorted(drug_verdicts.keys())

    # Group AMR joined hits by feature index.
    by_feature: dict[int, list] = defaultdict(list)
    for jh in joined_hits:
        if jh.feature_index is not None:
            by_feature[jh.feature_index].append(jh)

    # Group virulence joined hits by feature index (independent of the AMR join).
    by_virulence: dict[int, list] = defaultdict(list)
    for jh in (virulence_joined_hits or []):
        if jh.feature_index is not None:
            by_virulence[jh.feature_index].append(jh)

    rows = list(features.to_dict("records"))
    feature_maps: list[dict] = []
    per_tier = {t: 0 for t in TIER_PRECEDENCE}
    n_virulence_features = 0

    for i, r in enumerate(rows):
        product = str(r.get("product") or "")
        gene_symbol = str(r.get("gene_symbol") or "")
        hits_here = by_feature.get(i, [])
        high_hits = [jh for jh in hits_here if jh.is_high_confidence]
        fallback_hits = [jh for jh in hits_here if jh.join_confidence == SYMBOL_FALLBACK]

        vir_here = by_virulence.get(i, [])
        high_vir = [jh for jh in vir_here if jh.is_high_confidence]
        vir_fallback = [jh for jh in vir_here if jh.join_confidence == SYMBOL_FALLBACK]

        func_tier, func_reason = classify_feature_tier(product, gene_symbol)
        secondary: list[dict] = []
        phenotype: list[dict] = []
        virulence: list[dict] = []

        if high_hits:
            # AMR wins (Open Question B): determinant-phenotype, even if a VF hit also lands here.
            primary_tier = TIER_DETERMINANT_PHENOTYPE
            classification_reason = (
                f"high-confidence determinant join "
                f"({', '.join(sorted({jh.join_confidence for jh in high_hits}))})"
            )
            phenotype = _phenotype_for_feature(high_hits, drug_verdicts, drugs)
            secondary.append({"type": "molecular_function", "tier": func_tier, "reason": func_reason})
            # a co-located VF determinant is honest secondary evidence, never the AMR phenotype.
            for jh in high_vir:
                secondary.append({"type": "virulence_secondary", "vf_gene": jh.hit.symbol,
                                  "cluster": jh.hit.subclass or None, "join_confidence": jh.join_confidence})
        elif high_vir:
            # virulence-determinant tier (its own presence wall, the `virulence` field).
            primary_tier = TIER_VIRULENCE_DETERMINANT
            classification_reason = (
                f"high-confidence virulence determinant join "
                f"({', '.join(sorted({jh.join_confidence for jh in high_vir}))})"
            )
            virulence = [_virulence_field(jh, virulence_db_sha) for jh in high_vir]
            secondary.append({"type": "molecular_function", "tier": func_tier, "reason": func_reason})
            n_virulence_features += 1
        else:
            primary_tier = func_tier
            classification_reason = func_reason

        # Symbol-fallback AMR determinants are VISIBLE secondary evidence — never phenotype.
        for jh in fallback_hits:
            secondary.append({
                "type": "determinant_symbol_fallback",
                "symbol": jh.hit.symbol,
                "amrfinder_class": jh.hit.cls,
                "join_confidence": SYMBOL_FALLBACK,
                "note": "gene-symbol-only match (the 0%-overlap trap) — NOT a phenotype claim",
            })
        # Symbol-fallback VF hits are VISIBLE secondary evidence — never the virulence tier.
        for jh in vir_fallback:
            secondary.append({
                "type": "virulence_symbol_fallback",
                "vf_gene": jh.hit.symbol,
                "cluster": jh.hit.subclass or None,
                "join_confidence": SYMBOL_FALLBACK,
                "note": "gene-symbol-only VF match — NOT a located virulence determinant",
            })

        per_tier[primary_tier] += 1
        source_tool = "amrfinder+bakta" if (high_hits or fallback_hits) else "bakta"
        if high_vir or vir_fallback:
            source_tool = source_tool + "+virulencefinder" if source_tool != "bakta" else "virulencefinder+bakta"

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
            # phenotype is non-empty ONLY on determinant-phenotype features (the AMR wall).
            "phenotype": phenotype,
            # virulence is non-empty ONLY on virulence-determinant features (the VF wall).
            "virulence": virulence,
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
        # --- virulence overlay (M2: SEPARATE keys; never feed the AMR keys/gate above) ---
        "virulence_determinant_feature_count": n_virulence_features,
        "virulence_join_quality": dict(virulence_join_counts) if virulence_join_counts else {},
        "all_virulence_joins_symbol_fallback": bool(
            (virulence_join_counts or {}).get("all_virulence_joins_symbol_fallback", False)
        ),
        # the genome-level pathotype overlay (the virulence analog of the AMR R/S call),
        # reported SEPARATELY from both the per-feature tiers and the AMR calls.
        "genome_pathotype_call": pathotype_call,
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
