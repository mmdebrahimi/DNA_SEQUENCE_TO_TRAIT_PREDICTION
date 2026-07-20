"""AR Isolate Bank N. gonorrhoeae scorer — the NON-FROZEN gono rule vs measured BMD-MIC labels.

Mirrors `external_cohort_revalidate` (the E. coli/Klebsiella external arm) but routes to the curated
NON-FROZEN gonococcal rule `dna_decode/organism_rules/neisseria_amr.call_ng_amr(drug, symbols)` instead
of the frozen `call_resistance`. Detection = AMRFinderPlus `-O Neisseria_gonorrhoeae` (verified present in
ncbi/amr:4.2.7); the determinant symbols are parsed from main.tsv and fed to `call_ng_amr`.

Reuses the shared arm's predictor-agnostic machinery (leakage-manifest drift guard, one-strain smoke,
per-strain records, powering gate, artifact schema). A Kaggle-produced main.tsv scores WITHOUT re-running
Docker (the predictor reuses a cached AMRFinder dir when present).

Scope: the AR Bank gono panel is 7 drugs; 6 are scorable (azithromycin / cefixime / ceftriaxone /
ciprofloxacin / penicillin / tetracycline). **gentamicin ABSTAINS** (no validated determinant ->
call_ng_amr returns INDETERMINATE -> excluded from sens/spec). NAMESPACE-SEPARATE artifacts
(`wiki/ar_bank_gono_validation_<drug>_<run_id>_<date>.json`); CURATED_NONFROZEN, NOT the frozen deployed
surface (byte-unchanged). Endorsement falsifier per drug (mirrors the AMR-Portal cell): spec >= 0.85 on a
powered (>=10 in the informative class) provenance-disjoint set, else UNDERPOWERED/INDETERMINATE.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.organism_rules.neisseria_amr import call_ng_amr  # noqa: E402
from scripts.external_cohort_revalidate import (  # noqa: E402
    assert_manifest_alignment, build_artifact, conf_from_records, gate_ok,
    powering_gate, predict_records, rule_predictor, smoke_predict, _read_selected,
)

AMRFINDER_ORGANISM = "Neisseria_gonorrhoeae"
REGISTRY_ORGANISM = "Neisseria gonorrhoeae"
SCORABLE_DRUGS = ("azithromycin", "cefixime", "ceftriaxone", "ciprofloxacin", "penicillin", "tetracycline")
SPEC_FLOOR = 0.85


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", default="ar_bank_gono")
    ap.add_argument("--drug", required=True, help="one of the scorable gono drugs (gentamicin ABSTAINS)")
    ap.add_argument("--labels-dir", required=True, help="dir with selected_{strict,relaxed}.tsv + buckets_<drug>.json")
    ap.add_argument("--preflight-json", default=None, help="Step-1 preflight artifact (fail-closed if absent)")
    ap.add_argument("--allow-degraded", action="store_true")
    ap.add_argument("--min-per-class", type=int, default=10,
                    help="min informative-class N to count as powered (default 10; lower for a documented pilot)")
    ap.add_argument("--skip-smoke", action="store_true")
    ap.add_argument("--cohort-manifest", default=None,
                    help="cohort_manifest_external_<run_id>.json -> drift-guards selected.tsv subset of manifest")
    ap.add_argument("--run-id", default=_date.today().isoformat())
    a = ap.parse_args()

    if a.drug.strip().lower() not in SCORABLE_DRUGS:
        print(f"REFUSE: {a.drug!r} is not a scorable gono drug {SCORABLE_DRUGS} "
              f"(gentamicin ABSTAINS -- no validated determinant)")
        return 2

    preflight = json.loads(Path(a.preflight_json).read_text(encoding="utf-8")) if a.preflight_json else None
    ok, reason = gate_ok(preflight, a.allow_degraded)
    print(f"GATE: {reason}")
    if not ok:
        return 2

    labels_dir = Path(a.labels_dir)
    strict_labels = _read_selected(labels_dir / "selected_strict.tsv")
    relaxed_labels = _read_selected(labels_dir / "selected_relaxed.tsv")
    buckets_path = labels_dir / f"buckets_{a.drug}.json"
    buckets = json.loads(buckets_path.read_text(encoding="utf-8")) if buckets_path.exists() else {}

    if a.cohort_manifest:
        man = json.loads(Path(a.cohort_manifest).read_text(encoding="utf-8"))
        man_bs = man.get("biosamples") or sorted({r["biosample"] for r in man.get("rows", [])})
        assert_manifest_alignment(set(strict_labels) | set(relaxed_labels), man_bs)

    from dna_decode.data.external_cohort_genomes import resolve_cohort_genomes
    from dna_decode.eval.biosample_resolver import BioSampleResolver
    resolver = BioSampleResolver()
    genomes = resolve_cohort_genomes(set(strict_labels) | set(relaxed_labels), resolver)
    resolver.save_cache()

    base = Path(f"data/raw/{a.cohort}_extval_{a.drug}")
    own_runs, gcache = base / "amrfinder_runs", base / "refseq"
    reuse_glob = f"data/raw/{a.cohort}_*/amrfinder_runs"
    predict = rule_predictor(a.drug, own_runs, gcache, reuse_glob, call_ng_amr, AMRFINDER_ORGANISM)

    if not a.skip_smoke:
        smoke = smoke_predict(genomes["free"], predict)
        print(f"SMOKE: {smoke['gca']} -> {smoke['prediction']} (ok={smoke['ok']})")
        if not smoke["ok"] and not a.allow_degraded:
            print(f"SMOKE FAILED (--skip-smoke/--allow-degraded to bypass): {smoke['reason']}")
            return 3

    strict_records, strict_excl = predict_records(genomes["free"], strict_labels, predict)
    strict_conf = conf_from_records(strict_records, strict_excl)
    relaxed_records, relaxed_excl = predict_records(genomes["free"], relaxed_labels, predict)
    relaxed_conf = conf_from_records(relaxed_records, relaxed_excl)

    base.mkdir(parents=True, exist_ok=True)
    (base / "predictions_strict.json").write_text(json.dumps(strict_records, indent=2), encoding="utf-8")

    n_attempted = len(strict_records)
    n_indet = sum(1 for r in strict_records if str(r["prediction"]).upper() not in ("R", "S"))
    pg = powering_gate(strict_conf, n_attempted, n_indet, min_per_class=a.min_per_class)

    leakage_control = (f"BioSample+assembly-level preflight: "
                       f"{('PASS' if preflight and preflight.get('verdict') == 'PASS' else 'DEGRADED')}; "
                       f"{genomes['n_assembly_required']} reads-only excluded (ASSEMBLY-REQUIRED)")
    artifact = build_artifact(a.cohort, a.drug, strict=strict_conf, relaxed=relaxed_conf, buckets=buckets,
                              leakage_control=leakage_control,
                              degraded=not (preflight and preflight.get("verdict") == "PASS"),
                              powering=pg, run_degraded=pg["degraded"], run_id=a.run_id,
                              amrfinder_organism=AMRFINDER_ORGANISM, registry_organism=REGISTRY_ORGANISM)
    # Gono cell provenance overrides (non-frozen curated rule; its own honesty tier).
    artifact.update({
        "_schema": "ar-bank-gono-validation-v1",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "rule_text": call_ng_amr(a.drug, [])["rule"],
        "spec_floor": SPEC_FLOOR,
        "evidence_tier": "external_clinical_curated_nonfrozen",
        "independence_tier": ("CDC AR Isolate Bank measured BMD-MIC, provenance-disjoint (BioSample + "
                              "resolution-free assembly-base) vs the tuning cohorts; genotype = AMRFinderPlus "
                              "-O Neisseria_gonorrhoeae determinant calls; NOT methodology-independent (the "
                              "rule is curated from Pathogenwatch 485.toml + literature)."),
        "frozen_surface_changed": False,
    })

    out = Path(f"wiki/ar_bank_gono_validation_{a.drug}_{a.run_id}_{_date.today().isoformat()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"RESULT {a.drug}: strict acc={strict_conf['acc']} sens={strict_conf['sens']} "
          f"spec={strict_conf['spec']} (n={strict_conf['n_scored']}, scored {pg['scored_R']}R/{pg['scored_S']}S, "
          f"indeterminate {pg['indeterminate_fraction']}); artifact: {out}")
    if pg["hard_fail"]:
        print(f"POWERING HARD FAIL: {pg['reasons']}")
        return 3
    if pg["degraded"] and not a.allow_degraded:
        print(f"RUN DEGRADED (--allow-degraded to accept): {pg['reasons']}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
