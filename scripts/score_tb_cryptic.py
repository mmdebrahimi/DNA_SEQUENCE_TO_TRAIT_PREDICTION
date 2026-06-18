"""Step 5 — v1b lineage-collapsed TB scoring orchestrator + cohort/callability gate (deliverable a).

Pure scoring (`score_cohort`) takes {preds, labels, clusters} and reuses the FROZEN
`clonality.cluster_weighted_confusion` (brainstorm C2a — NOT representative-dedup). It emits the
lineage-collapsed sens/spec + raw sens/spec + raw->lineage shrinkage + n_discordant +
n_clusters_mixed_prediction (M1) + Wilson CI + effective_lineage_n + n_uncallable, tagged
`WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

Status gate (C3 + honesty):
  - lineage assignment unavailable (no clusters / all UNASSIGNED-degenerate) ->
    `LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL` (never a raw-only headline).
  - scored set is a convenience subset (not the full prevalence-preserving per-drug cohort) ->
    `TB_SUBSET_PLUMBING` (metrics computed but NEVER the baseline label).
  - else -> `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

The actual cohort fetch (full per-drug masked + regeno VCFs, ~1.6 TB -> D:) is the CLI runtime path
(`main`); it is NOT exercised by unit tests, which drive the pure scoring with synthetic inputs.
"""
from __future__ import annotations

from dna_decode.eval.clonality import (
    cluster_members,
    cluster_weighted_confusion,
    effective_lineage_n,
    wilson_ci,
)
from dna_decode.organism_rules import tb_amr, tb_lineage, tb_vcf

BASELINE_LABEL = "WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE"
PLUMBING_LABEL = "TB_SUBSET_PLUMBING"
BLOCKED_LABEL = "LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL"

_R, _S, _ABSTAIN = "R", "S", "ABSTAIN"


def n_clusters_mixed_prediction(preds: dict[str, str], clusters: dict[str, int]) -> int:
    """M1: same-cluster member predictions disagree (within-lineage determinant heterogeneity)."""
    n = 0
    for sids in cluster_members(clusters).values():
        calls = {preds.get(s, _ABSTAIN) for s in sids} & {_R, _S}
        if calls == {_R, _S}:
            n += 1
    return n


def raw_confusion(preds: dict[str, str], labels: dict[str, str]) -> dict:
    """Per-isolate confusion (ABSTAIN excluded) — the un-collapsed, clonality-inflated baseline."""
    tp = fp = tn = fn = 0
    for sid, lab in labels.items():
        p = preds.get(sid, _ABSTAIN)
        L = str(lab).upper()
        if p == _R and L == _R:
            tp += 1
        elif p == _R and L == _S:
            fp += 1
        elif p == _S and L == _S:
            tn += 1
        elif p == _S and L == _R:
            fn += 1
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def _lineage_available(clusters: dict[str, int]) -> bool:
    # degenerate iff every isolate is its own singleton (no real collapse happened)
    if not clusters:
        return False
    return len(set(clusters.values())) < len(clusters)


def score_cohort(
    preds: dict[str, str],
    labels: dict[str, str],
    clusters: dict[str, int],
    *,
    drug: str,
    cohort_complete: bool,
) -> dict:
    """Lineage-collapsed scoring + status gate. Returns a self-describing result dict."""
    n_abstain = sum(1 for s in labels if preds.get(s, _ABSTAIN) == _ABSTAIN)
    base = {
        "drug": drug,
        "rule_status": "deterministic",
        "n_isolates": len(labels),
        "n_uncallable_abstain": n_abstain,
    }

    if not _lineage_available(clusters):
        return {**base, "status": BLOCKED_LABEL,
                "reason": "no non-singleton lineage clusters — lineage assignment unavailable"}

    cw = cluster_weighted_confusion(preds, labels, clusters)
    raw = raw_confusion(preds, labels)
    eff_r = effective_lineage_n(clusters, labels, "R")
    eff_s = effective_lineage_n(clusters, labels, "S")
    raw_r = sum(1 for s, l in labels.items() if str(l).upper() == _R)
    raw_s = sum(1 for s, l in labels.items() if str(l).upper() == _S)

    result = {
        **base,
        "status": BASELINE_LABEL if cohort_complete else PLUMBING_LABEL,
        "honesty": ("In-distribution knowledge-baseline: the WHO catalogue was built partly from "
                    "CRyPTIC. NOT independent validation — see the separate post-2023 gold-set arm."),
        "lineage_collapsed": {
            "sens": cw["sens"], "spec": cw["spec"],
            "tp": cw["tp"], "fp": cw["fp"], "tn": cw["tn"], "fn": cw["fn"],
            "sens_wilson_ci": wilson_ci(cw["tp"], cw["tp"] + cw["fn"]),
            "spec_wilson_ci": wilson_ci(cw["tn"], cw["tn"] + cw["fp"]),
            "n_clusters_R": cw["n_clusters_R"], "n_clusters_S": cw["n_clusters_S"],
            "n_discordant": cw["n_discordant"], "n_cluster_abstain": cw["n_cluster_abstain"],
            "n_clusters_mixed_prediction": n_clusters_mixed_prediction(preds, clusters),
        },
        "raw": raw,
        "effective_lineage_n": {"R": eff_r, "S": eff_s},
        "raw_to_lineage_shrinkage": {"R": [raw_r, eff_r], "S": [raw_s, eff_s]},
    }
    if not cohort_complete:
        result["plumbing_note"] = ("convenience/partial cohort — metrics are NOT the baseline; "
                                   "fetch the full per-drug prevalence-preserving cohort to earn the label")
    return result


def run_v1b(
    strain_masked: dict[str, str],
    strain_label: dict[str, str],
    determinants,
    barcode,
    *,
    drug: str,
    cohort_complete: bool,
    strain_regeno: dict[str, str] | None = None,
) -> dict:
    """End-to-end v1b: masked VCFs -> determinant calls + barcode lineage -> lineage-collapsed score.

    `strain_regeno` (optional) supplies per-isolate regeno text for the callability gate; absent ->
    score_drug runs in callability-unassessed mode (S, not ABSTAIN) and the result flags it.
    """
    strain_regeno = strain_regeno or {}
    calls_by_strain = {sid: tb_vcf.parse_masked_calls(txt) for sid, txt in strain_masked.items()}
    preds = {
        sid: tb_amr.score_drug(drug, calls, determinants,
                               regeno_text=strain_regeno.get(sid)).prediction
        for sid, calls in calls_by_strain.items()
    }
    clusters = tb_lineage.lineage_clusters(calls_by_strain, barcode)
    out = score_cohort(preds, strain_label, clusters, drug=drug, cohort_complete=cohort_complete)
    out["callability_assessed"] = bool(strain_regeno)
    out["n_unassigned_lineage"] = sum(
        1 for v in tb_lineage.lineage_assignments(calls_by_strain, barcode).values()
        if v == tb_lineage.UNASSIGNED
    )
    return out


def _main() -> int:
    """CLI: score a cached per-drug cohort. PoC default = the cached masked subset (-> PLUMBING)."""
    import argparse
    import json
    from datetime import date as _date
    from pathlib import Path

    from dna_decode.data import tb_lineage_barcode, tb_who_catalogue
    from dna_decode.organism_rules import tb_vcf as _vcf

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="rifampicin")
    ap.add_argument("--masked-cache", default="data/raw/cryptic/vcf_cache",
                    help="dir of cached DECOMPRESSED masked VCFs (PoC) — or a D: regeno-aware cohort dir")
    ap.add_argument("--regeno-cache", default=None, help="dir of cached regeno VCFs (callability)")
    ap.add_argument("--cohort-complete", action="store_true",
                    help="assert this IS the full prevalence-preserving per-drug cohort (earns BASELINE)")
    ap.add_argument("--max", type=int, default=0, help="cap isolates (0 = all cached)")
    a = ap.parse_args()

    code = {"rifampicin": "RIF", "isoniazid": "INH"}[a.drug]
    rows = _vcf.reuse_rows()
    masked_dir = Path(a.masked_cache)
    regeno_dir = Path(a.regeno_cache) if a.regeno_cache else None

    def _cached(rel: str, d: Path) -> str | None:
        f = d / (rel.replace("/", "_").replace("..", "").strip("_"))
        return f.read_text(encoding="utf-8", errors="replace") if f.exists() else None

    strain_masked, strain_label, strain_regeno = {}, {}, {}
    for r in rows:
        ph = (r.get(f"{code}_BINARY_PHENOTYPE") or "").strip().upper()
        q = (r.get(f"{code}_PHENOTYPE_QUALITY") or "").strip().upper()
        masked_rel, regeno_rel = _vcf.vcf_paths_for(r)
        if ph not in ("R", "S") or q != "HIGH" or not masked_rel:
            continue
        txt = _cached(masked_rel, masked_dir)
        if txt is None:
            continue
        sid = (r.get("UNIQUEID") or r.get("ENA_RUN") or masked_rel).strip()
        strain_masked[sid] = txt
        strain_label[sid] = ph
        if regeno_dir and regeno_rel:
            rg = _cached(regeno_rel, regeno_dir)
            if rg is not None:
                strain_regeno[sid] = rg
        if a.max and len(strain_masked) >= a.max:
            break

    if not strain_masked:
        print(f"no cached masked VCFs found for {a.drug} under {masked_dir}")
        return 2

    tb_who_catalogue.verify_pins()
    dets = tb_who_catalogue.load_determinants(a.drug)
    barcode = tb_lineage_barcode.load_barcode()
    res = run_v1b(strain_masked, strain_label, dets, barcode, drug=a.drug,
                  cohort_complete=a.cohort_complete, strain_regeno=strain_regeno or None)

    out = Path(f"wiki/tb_{code.lower()}_cryptic_results_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"status={res['status']}  n_isolates={res['n_isolates']}  "
          f"callability_assessed={res.get('callability_assessed')}")
    lc = res.get("lineage_collapsed")
    if lc:
        print(f"  lineage-collapsed: sens={lc['sens']} spec={lc['spec']} "
              f"(R-lineages={lc['n_clusters_R']} S-lineages={lc['n_clusters_S']} "
              f"discordant={lc['n_discordant']} mixed-pred={lc['n_clusters_mixed_prediction']})")
        print(f"  raw: sens={res['raw']['sens']} spec={res['raw']['spec']} "
              f"shrinkage R{res['raw_to_lineage_shrinkage']['R']} S{res['raw_to_lineage_shrinkage']['S']}")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
