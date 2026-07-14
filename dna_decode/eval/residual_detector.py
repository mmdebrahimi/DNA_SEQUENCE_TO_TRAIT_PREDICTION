"""Residual-signal detector — the de-confounding machinery promoted to a first-class product (2026-07-13).

`/innovate` round-1 survivor **g8-residual-detector**: the leave-one-clade-out (Mash-clade GroupKFold)
de-confound is currently a hidden guard buried in `scripts/crossaxis_lineage_deconfound.py`. This module
turns it into a REUSABLE decoder output: given a cross-axis lineage-deconfound artifact, report which
genome features carry mechanism signal *that survives phylogeny removal* vs which are clade-mediated
(possible clonal-structure artifacts).

THE HONESTY WALL (load-bearing): a residual-signal tier is a SIGNAL-PROVENANCE statement, NOT a phenotype
prediction. `GENERALIZES` = the feature's genotype->axis association survives leave-one-clade-out CV (real
mechanism signal beyond lineage). `LINEAGE_MEDIATED` = clade-concentrated: the association MAY be clonal
population structure, not a real linkage. `UNTESTED` = no cross-axis entry. The report never emits an R/S
call — it reports where the de-confounded signal lives.

Pure + dependency-free. This is the canonical home for the tier vocabulary; `dna_decode.viz.network_adapter`
re-uses `classify_residual` (the same classification the network browser draws as solid/dashed borders).
"""
from __future__ import annotations

import json
from pathlib import Path

# residual-signal tiers (the de-confound verdict per feature)
GENERALIZES = "generalizes"            # survives leave-one-clade-out CV -> real signal beyond lineage
LINEAGE_MEDIATED = "lineage_mediated"  # clade-concentrated -> may be clonal population structure
UNTESTED = "untested"                  # no cross-axis entry for this feature


def classify_residual(entry: dict | None) -> str:
    """Map one cross-axis per_gene record -> a residual-signal tier. Pure.

    A record that is clade-concentrated is LINEAGE_MEDIATED regardless of its naive AUC. A record that
    generalizes-beyond-lineage is GENERALIZES. A record that was TESTED but neither flag is set (it
    dropped below the generalization bar under clade-grouping without being flagged clade-concentrated)
    is treated as LINEAGE_MEDIATED — the conservative call, since its signal did not survive de-confounding.
    """
    if not entry:
        return UNTESTED
    if entry.get("clade_concentrated"):
        return LINEAGE_MEDIATED
    if entry.get("generalizes_beyond_lineage"):
        return GENERALIZES
    return LINEAGE_MEDIATED


def _family_of(feature_id: str) -> str:
    """Coarse gene-family bucket for the region rollup (grouping aid, NOT a claim)."""
    fid = feature_id
    for pref, fam in (("gyrA", "QRDR"), ("gyrB", "QRDR"), ("parC", "QRDR"), ("parE", "QRDR"),
                      ("bla", "beta-lactam"), ("CTX", "beta-lactam"), ("TEM", "beta-lactam"),
                      ("SHV", "beta-lactam"), ("KPC", "beta-lactam"), ("NDM", "beta-lactam"),
                      ("OXA", "beta-lactam"), ("aac", "aminoglycoside"), ("aad", "aminoglycoside"),
                      ("aph", "aminoglycoside"), ("ant", "aminoglycoside"), ("arm", "aminoglycoside"),
                      ("rmt", "aminoglycoside"), ("sul", "sulfa/trimethoprim"), ("dfr", "sulfa/trimethoprim"),
                      ("tet", "tet/macrolide"), ("erm", "tet/macrolide"), ("mph", "tet/macrolide"),
                      ("mef", "tet/macrolide"), ("qnr", "quinolone-plasmid"), ("oqx", "quinolone-plasmid"),
                      ("catch", "phenicol"), ("cat", "phenicol"), ("fos", "fosfomycin"), ("mcr", "colistin")):
        if fid.startswith(pref):
            return fam
    return "other"


def build_residual_report(artifact: dict) -> dict:
    """Turn a cross-axis lineage-deconfound artifact into a residual-signal report. Pure.

    Report = {meta, per_feature[], tier_counts, family_rollup, honest_caveats}. per_feature is sorted
    GENERALIZES-first then by de-confounded AUC (clade-grouped) descending, so the strongest residual
    signals lead.
    """
    per_gene = artifact.get("per_gene", {}) or {}
    rows = []
    for fid, entry in per_gene.items():
        tier = classify_residual(entry)
        rows.append({
            "feature_id": fid,
            "tier": tier,
            "family": _family_of(fid),
            "n_present": entry.get("n_present"),
            "auc_naive": entry.get("auc_naive"),
            "auc_clade_grouped": entry.get("auc_clade_grouped"),
            "drop": entry.get("drop"),
        })
    # GENERALIZES first, then by de-confounded strength (clade-grouped AUC) descending
    _rank = {GENERALIZES: 0, LINEAGE_MEDIATED: 1, UNTESTED: 2}
    rows.sort(key=lambda r: (_rank.get(r["tier"], 3), -(r["auc_clade_grouped"] or 0.0), r["feature_id"]))

    tier_counts = {GENERALIZES: 0, LINEAGE_MEDIATED: 0, UNTESTED: 0}
    for r in rows:
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1

    # family rollup: how much residual (GENERALIZES) signal each gene family carries
    family_rollup: dict[str, dict] = {}
    for r in rows:
        fam = family_rollup.setdefault(r["family"], {GENERALIZES: 0, LINEAGE_MEDIATED: 0, UNTESTED: 0, "total": 0})
        fam[r["tier"]] += 1
        fam["total"] += 1

    return {
        "meta": {
            "axis_label": artifact.get("axis_label"),
            "organism": artifact.get("organism"),
            "source_verdict": artifact.get("verdict"),
            "median_auc_naive": artifact.get("median_auc_naive"),
            "median_auc_clade_grouped": artifact.get("median_auc_clade_grouped"),
            "mash_threshold": (artifact.get("mash", {}) or {}).get("threshold"),
            "n_clades": (artifact.get("mash", {}) or {}).get("n_clades"),
            "largest_clade_frac": (artifact.get("mash", {}) or {}).get("largest_clade_frac"),
            "n_features": len(rows),
        },
        "per_feature": rows,
        "tier_counts": tier_counts,
        "family_rollup": family_rollup,
        "honest_caveats": [
            "A residual-signal tier is a SIGNAL-PROVENANCE statement, NOT a phenotype prediction — the "
            "report says WHERE de-confounded mechanism signal lives, never R/S.",
            "GENERALIZES = association survives leave-one-clade-out CV (real signal beyond lineage). "
            "LINEAGE_MEDIATED = clade-concentrated: MAY be clonal population structure, not a real linkage.",
            "De-confounding quality is bounded by the Mash-clade partition (threshold + largest-clade "
            "fraction in meta); a coarse partition under-removes lineage.",
        ] + list(artifact.get("honest_caveats", []) or []),
    }


def load_report(artifact_path: str | Path) -> dict:
    """Read a cross-axis lineage-deconfound JSON from disk and build its residual-signal report."""
    art = json.loads(Path(artifact_path).read_text(encoding="utf-8"))
    return build_residual_report(art)
