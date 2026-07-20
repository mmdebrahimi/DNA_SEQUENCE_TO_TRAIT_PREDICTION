"""Generalized AR Isolate Bank scorer — any organism in `ar_bank_registry`.

Generalizes `ar_bank_gono_validate`: the registry supplies the AMRFinder `-O` organism + a uniform
`rule_fn(drug, symbols)` dispatching to that organism's NON-FROZEN cell. Reuses the shared arm's
predictor-agnostic machinery (rule_predictor / powering gate / manifest drift-guard / smoke / artifact).
A Kaggle-produced main.tsv scores WITHOUT re-running Docker. Frozen surface untouched.

  uv run python -m scripts.ar_bank_organism_validate --organism enterococcus_faecium --drug levofloxacin \
     --labels-dir data/raw/ar_bank_enterococcus_faecium_extval_levofloxacin --cohort-manifest <manifest> \
     --allow-degraded --min-per-class 5
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.organism_rules.ar_bank_registry import config_for, rule_fn_for
from scripts.external_cohort_revalidate import (
    assert_manifest_alignment, build_artifact, conf_from_records, gate_ok,
    powering_gate, predict_records, rule_predictor, smoke_predict, _read_selected,
)

SPEC_FLOOR = 0.85


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, help="registry key (e.g. enterococcus_faecium)")
    ap.add_argument("--drug", required=True)
    ap.add_argument("--labels-dir", required=True)
    ap.add_argument("--preflight-json", default=None)
    ap.add_argument("--allow-degraded", action="store_true")
    ap.add_argument("--min-per-class", type=int, default=5)
    ap.add_argument("--skip-smoke", action="store_true")
    ap.add_argument("--cohort-manifest", default=None)
    ap.add_argument("--run-id", default=_date.today().isoformat())
    a = ap.parse_args()
    cfg = config_for(a.organism)
    if a.drug.strip().lower() not in {c.lower() for c in cfg["drug_map"]}:
        print(f"REFUSE: {a.drug!r} not scorable for {a.organism} (scorable: {sorted(cfg['drug_map'])})")
        return 2
    amrfinder_organism = cfg["amrfinder_organism"]
    registry_organism = cfg["registry_organism"]
    rule_fn = rule_fn_for(a.organism)

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

    base = Path(f"data/raw/ar_bank_{a.organism}_extval_{a.drug}")
    own_runs, gcache = base / "amrfinder_runs", base / "refseq"
    reuse_glob = f"data/raw/ar_bank_{a.organism}_*/amrfinder_runs"
    predict = rule_predictor(a.drug, own_runs, gcache, reuse_glob, rule_fn, amrfinder_organism)

    if not a.skip_smoke:
        smoke = smoke_predict(genomes["free"], predict)
        print(f"SMOKE: {smoke['gca']} -> {smoke['prediction']} (ok={smoke['ok']})")
        if not smoke["ok"] and not a.allow_degraded:
            print(f"SMOKE FAILED: {smoke['reason']}")
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
                       f"{genomes['n_assembly_required']} reads-only excluded")
    artifact = build_artifact(f"ar_bank_{a.organism}", a.drug, strict=strict_conf, relaxed=relaxed_conf,
                              buckets=buckets, leakage_control=leakage_control,
                              degraded=not (preflight and preflight.get("verdict") == "PASS"),
                              powering=pg, run_degraded=pg["degraded"], run_id=a.run_id,
                              amrfinder_organism=amrfinder_organism, registry_organism=registry_organism)
    artifact.update({
        "_schema": "ar-bank-organism-validation-v1", "organism_key": a.organism,
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "rule_text": rule_fn(a.drug, [])["rule"], "spec_floor": SPEC_FLOOR,
        "registry_note": cfg["note"], "frozen_surface_changed": False,
    })
    out = Path(f"wiki/ar_bank_{a.organism}_validation_{a.drug}_{a.run_id}_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"RESULT {a.organism}/{a.drug}: strict acc={strict_conf['acc']} sens={strict_conf['sens']} "
          f"spec={strict_conf['spec']} (n={strict_conf['n_scored']}, {pg['scored_R']}R/{pg['scored_S']}S, "
          f"indet {pg['indeterminate_fraction']}); artifact: {out}")
    if pg["hard_fail"]:
        print(f"POWERING HARD FAIL: {pg['reasons']}")
        return 3
    if pg["degraded"] and not a.allow_degraded:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
