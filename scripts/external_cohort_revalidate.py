"""External-cohort scorer — re-validate the FROZEN decoder on an independent cohort.

Mirrors the canonical `scripts/provenance_disjoint_validate.py` per-strain loop
(ensure_run -> call_resistance -> independent_cohort_validate._conf), but sources
the cohort from an external measured-MIC project + writes to a SEPARATE
`wiki/external_validation_*` namespace so the FROZEN report card / lineage layer
(which key cells by canonical_cell_key(organism,drug)) are NOT collided into.

The decoder is reused UNCHANGED. The organism triple is taken VERBATIM from the
frozen E. coli provdisjoint cells so external numbers are comparable:

  AMRFINDER_ORGANISM = "Escherichia"             (AMRFinder -O; what the frozen cells used)
  REGISTRY_ORGANISM  = "Escherichia_coli_Shigella" (call_resistance organism=; falls through
                       to the validated DRUG_RULE default — Escherichia has NO calibrated entry)

Source of these strings: wiki/provenance_disjoint_validation_escherichia_coli_shigella_*.json
(amrfinder_organism / registry_organism fields). Do NOT invent — see plan Risk Flags.

Both STRICT (HIGH_R/HIGH_S; primary) and RELAXED (+DECISIVE; secondary) label sets
are scored. Fail-closed unless the preflight PASSED (or --allow-degraded).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path
from typing import Callable

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.independent_cohort_validate import _conf

# Verbatim from the frozen E. coli provdisjoint cells (read, NOT invented).
AMRFINDER_ORGANISM = "Escherichia"
REGISTRY_ORGANISM = "Escherichia_coli_Shigella"

# Powering floor — an external CLINICAL claim must be powered, not just non-empty.
# Mirrors independent_cohort_validate.py's promotion gate (MIN_PER_CLASS=10).
MIN_PER_CLASS = 10
MAX_INDETERMINATE_FRACTION = 0.20   # of attempted-FREE strains (reads-only exclusions don't count)

EVIDENCE_TIER = "external_clinical"
INDEPENDENCE_TIER = (
    "external clinical cohort — different country / lab / AST method than the decoder's "
    "US-NCBI-PD tuning provenance; clonality/lineage disclosed in the roll-up. NOT the same "
    "as methodology-independent: comparability rests on the same AMRFinder -O + DRUG_RULE."
)

Predict = Callable[[str], str]


def real_predictor(drug: str, own_runs: Path, gcache: Path, reuse_glob: str,
                   amrfinder_organism: str = AMRFINDER_ORGANISM,
                   registry_organism: str = REGISTRY_ORGANISM) -> Predict:
    """Build the real predict(gca) -> prediction using the SHIPPED toolchain.

    Mirrors provenance_disjoint_validate: reuse an AMRFinder cache dir if present,
    else download+run via ensure_run, then apply the frozen call_resistance with the
    pinned organism. Network + Docker; exercised only in the manual full run.

    The organism triple defaults to the frozen E. coli cell's values; pass the target
    organism's VERBATIM triple (read from its frozen provdisjoint cell) for Klebsiella etc.
    """
    from dna_decode.eval.amr_rules import call_resistance
    from scripts.organism_drug_validate import _run_dir, ensure_run

    def predict(gca: str) -> str:
        mt = _run_dir(gca, own_runs, reuse_glob)
        if mt is None:
            ensure_run(gca, own_runs, gcache, amrfinder_organism, reuse_glob)
            mt = _run_dir(gca, own_runs, reuse_glob)
        if mt is None:
            return "INDETERMINATE"
        return call_resistance(mt / "main.tsv", drug, organism=registry_organism)["prediction"]

    return predict


def parse_determinant_symbols(main_tsv: Path) -> list[str]:
    """Read an AMRFinderPlus main.tsv -> the list of determinant symbols (the 'Element symbol' column;
    'Gene symbol' on older AMRFinder). These feed the organism_rules `call_<org>_amr(drug, symbols)`
    rules (POINT mutations like gyrA_S91F + acquired genes like tet(M)/blaTEM). Empty if the file is
    absent/empty."""
    import csv
    if not main_tsv.exists():
        return []
    out: list[str] = []
    with main_tsv.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym = (row.get("Element symbol") or row.get("Gene symbol")
                   or row.get("amr_element_symbol") or "").strip()
            if sym:
                out.append(sym)
    return out


def rule_predictor(drug: str, own_runs: Path, gcache: Path, reuse_glob: str, rule_fn,
                   amrfinder_organism: str) -> Predict:
    """Like `real_predictor`, but instead of the frozen `call_resistance` it parses the AMRFinder
    determinant symbols and applies a NON-FROZEN organism rule `rule_fn(drug, symbols) -> {"prediction": ...}`
    (e.g. `organism_rules.neisseria_amr.call_ng_amr`). Reuses a cached AMRFinder dir if present (so a
    Kaggle-produced main.tsv scores without re-running Docker), else downloads+runs with -O
    `amrfinder_organism`. A missing run or an ABSTAIN rule -> INDETERMINATE (excluded from sens/spec)."""
    from scripts.organism_drug_validate import _run_dir, ensure_run

    def predict(gca: str) -> str:
        mt = _run_dir(gca, own_runs, reuse_glob)
        if mt is None:
            ensure_run(gca, own_runs, gcache, amrfinder_organism, reuse_glob)
            mt = _run_dir(gca, own_runs, reuse_glob)
        if mt is None:
            return "INDETERMINATE"
        return rule_fn(drug, parse_determinant_symbols(mt / "main.tsv"))["prediction"]

    return predict


def predict_records(free_genomes: dict[str, str], labels: dict[str, str],
                    predict: Predict) -> tuple[list[dict], int]:
    """ONE predict() call per FREE strain -> per-strain records (the Step-5 substrate).

    Each record: {biosample, gca, prediction, label, y}. BioSamples not in
    `free_genomes` (ASSEMBLY-REQUIRED) are excluded; returns (records, n_excluded).
    """
    records: list[dict] = []
    n_excluded = 0
    for bs, rs in sorted(labels.items()):
        gca = free_genomes.get(bs)
        if gca is None:
            n_excluded += 1
            continue
        label = str(rs).strip().upper()
        records.append({"biosample": bs, "gca": gca, "prediction": predict(gca),
                        "label": label, "y": 1 if label == "R" else 0})
    return records, n_excluded


def conf_from_records(records: list[dict], n_excluded: int = 0) -> dict:
    """`independent_cohort_validate._conf` over per-strain records + exclusion count."""
    conf = _conf([(r["prediction"], r["y"]) for r in records])
    conf["n_excluded_no_assembly"] = n_excluded
    return conf


def smoke_predict(free_genomes: dict[str, str], predict: Predict) -> dict:
    """Fail-fast: run ONE FREE strain through the live AMRFinder/decoder path BEFORE
    the full loop, so a broken Docker mount / missing DB / NCBI failure trips on
    strain 1 (seconds) instead of after N AMRFinder runs (hours). ok iff the one
    prediction is a real R/S call (INDETERMINATE means the live path is broken).
    """
    if not free_genomes:
        return {"ok": False, "reason": "no FREE genomes to smoke", "gca": None, "prediction": None}
    gca = sorted(free_genomes.values())[0]
    pred = predict(gca)
    ok = str(pred).upper() in ("R", "S")
    return {"ok": ok, "gca": gca, "prediction": pred,
            "reason": "" if ok else f"smoke strain {gca} -> {pred} (live AMRFinder/decoder path broken)"}


def score_label_set(free_genomes: dict[str, str], labels: dict[str, str],
                    predict: Predict) -> dict:
    """Score one label set (thin wrapper): predict each FREE strain -> _conf.

    BioSamples not in `free_genomes` (ASSEMBLY-REQUIRED) are excluded + counted.
    """
    records, n_excluded = predict_records(free_genomes, labels, predict)
    return conf_from_records(records, n_excluded)


def powering_gate(strict_conf: dict, n_attempted_free: int, n_indeterminate: int, *,
                  min_per_class: int = MIN_PER_CLASS,
                  max_indeterminate_frac: float = MAX_INDETERMINATE_FRACTION) -> dict:
    """Decide whether the strict-tier result is a POWERED, non-fail-open external claim.

    - HARD FAIL (not overridable): n_scored == 0, OR strict scored R/S below min_per_class
      (an underpowered claim — lower --min-per-class for a documented pilot, don't override).
    - DEGRADED (overridable by --allow-degraded): indeterminate fraction (of attempted-FREE
      strains, so reads-only ASSEMBLY-REQUIRED exclusions don't count) exceeds the threshold —
      a sign of broken Docker/NCBI rather than a clean run.
    """
    n_scored = strict_conf.get("n_scored", 0)
    scored_R = strict_conf.get("tp", 0) + strict_conf.get("fn", 0)
    scored_S = strict_conf.get("tn", 0) + strict_conf.get("fp", 0)
    frac = round(n_indeterminate / n_attempted_free, 4) if n_attempted_free else 0.0
    reasons: list[str] = []
    hard_fail = False
    if n_scored == 0:
        hard_fail = True
        reasons.append("n_scored == 0 (no FREE strain produced an R/S call — likely broken Docker/NCBI)")
    if scored_R < min_per_class:
        hard_fail = True
        reasons.append(f"strict scored R {scored_R} < min_per_class {min_per_class} (underpowered)")
    if scored_S < min_per_class:
        hard_fail = True
        reasons.append(f"strict scored S {scored_S} < min_per_class {min_per_class} (underpowered)")
    degraded = frac > max_indeterminate_frac
    if degraded:
        reasons.append(f"indeterminate fraction {frac} > {max_indeterminate_frac} "
                       f"(of {n_attempted_free} attempted-FREE)")
    return {"hard_fail": hard_fail, "degraded": degraded, "reasons": reasons,
            "n_attempted_free": n_attempted_free, "n_indeterminate": n_indeterminate,
            "indeterminate_fraction": frac, "scored_R": scored_R, "scored_S": scored_S,
            "min_per_class": min_per_class, "max_indeterminate_fraction": max_indeterminate_frac}


def build_artifact(cohort: str, drug: str, *, strict: dict, relaxed: dict, buckets: dict,
                   leakage_control: str, degraded: bool = False,
                   powering: dict | None = None, run_degraded: bool = False,
                   run_id: str | None = None,
                   amrfinder_organism: str = AMRFINDER_ORGANISM,
                   registry_organism: str = REGISTRY_ORGANISM) -> dict:
    """Assemble the external-validation-v1 artifact (separate namespace)."""
    return {
        "_schema": "external-validation-v1",
        "date": _date.today().isoformat(),
        "run_id": run_id,
        "cohort": cohort,
        "organism": registry_organism,
        "amrfinder_organism": amrfinder_organism,
        "registry_organism": registry_organism,
        "drug": drug,
        "evidence_tier": EVIDENCE_TIER,
        "independence_tier": INDEPENDENCE_TIER,
        "strict": strict,
        "relaxed": relaxed,
        "buckets": buckets,
        "leakage_control": leakage_control,
        "independence_degraded": bool(degraded),
        "powering": powering or {},
        "run_degraded": bool(run_degraded),
        "primary_metric": "strict",
    }


def gate_ok(preflight: dict | None, allow_degraded: bool) -> tuple[bool, str]:
    """Decide whether scoring may proceed. Fail-closed unless preflight PASSED."""
    if preflight is None:
        if allow_degraded:
            return True, "no preflight artifact; proceeding DEGRADED (--allow-degraded)"
        return False, "no preflight artifact — refusing to score (fail-closed; pass --allow-degraded)"
    if preflight.get("verdict") == "PASS":
        return True, "preflight PASS"
    if allow_degraded:
        return True, f"preflight {preflight.get('verdict')}; proceeding DEGRADED (--allow-degraded)"
    return False, f"preflight {preflight.get('verdict')} ({preflight.get('reasons')}) — fail-closed"


class ManifestDriftError(RuntimeError):
    """selected.tsv BioSamples are not a subset of the cohort manifest (drift)."""


def assert_manifest_alignment(selected_biosamples, manifest_biosamples) -> None:
    """Drift guard: every scored BioSample MUST be in the cohort manifest.

    Guarantees the leakage/availability verdict (computed by the exact-set preflight
    over the manifest) covers EXACTLY what gets scored — the split that the brainstorm
    flagged cannot reappear.
    """
    extra = sorted(set(selected_biosamples) - set(manifest_biosamples))
    if extra:
        raise ManifestDriftError(
            f"{len(extra)} scored BioSample(s) absent from the cohort manifest: {extra[:5]}"
            f"{'...' if len(extra) > 5 else ''}")


def _read_selected(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if path.exists():
        for ln in path.read_text(encoding="utf-8").splitlines():
            if "\t" in ln:
                a, rs = ln.split("\t")
                out[a.strip()] = rs.strip()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", required=True, help="cohort slug (e.g. spain_probac)")
    ap.add_argument("--drug", required=True, help="canonical pilot drug")
    ap.add_argument("--labels-dir", required=True,
                    help="dir with selected_strict.tsv + selected_relaxed.tsv + buckets_<drug>.json (Step 2)")
    ap.add_argument("--preflight-json", default=None, help="Step 1 preflight artifact (fail-closed if absent)")
    ap.add_argument("--allow-degraded", action="store_true",
                    help="proceed even if preflight FAILED/absent OR the run is indeterminate-degraded")
    ap.add_argument("--min-per-class", type=int, default=MIN_PER_CLASS,
                    help=f"min strict-scored R and S to count as powered (default {MIN_PER_CLASS}; "
                         f"lower for a documented small pilot)")
    ap.add_argument("--skip-smoke", action="store_true",
                    help="skip the one-strain fail-fast smoke before the full scoring loop")
    ap.add_argument("--cohort-manifest", default=None,
                    help="cohort_manifest_external_<run_id>.json — drift-guards selected.tsv "
                         "BioSamples subset of the manifest (required for a live run)")
    ap.add_argument("--run-id", default=_date.today().isoformat(),
                    help="run id stamped into the artifact + filename for run-scoped roll-up")
    ap.add_argument("--amrfinder-organism", default=AMRFINDER_ORGANISM,
                    help="AMRFinder -O organism (VERBATIM from the target's frozen provdisjoint cell; "
                         f"default {AMRFINDER_ORGANISM})")
    ap.add_argument("--registry-organism", default=REGISTRY_ORGANISM,
                    help="call_resistance organism= (VERBATIM from the target's frozen cell; "
                         f"default {REGISTRY_ORGANISM})")
    a = ap.parse_args()

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

    # Drift guard: selected.tsv BioSamples MUST be a subset of the cohort manifest.
    if a.cohort_manifest:
        man = json.loads(Path(a.cohort_manifest).read_text(encoding="utf-8"))
        man_bs = man.get("biosamples") or sorted({r["biosample"] for r in man.get("rows", [])})
        assert_manifest_alignment(set(strict_labels) | set(relaxed_labels), man_bs)

    # Resolve genomes for the union of labeled BioSamples, then score each set.
    from dna_decode.data.external_cohort_genomes import resolve_cohort_genomes
    from dna_decode.eval.biosample_resolver import BioSampleResolver

    resolver = BioSampleResolver()
    all_bs = set(strict_labels) | set(relaxed_labels)
    genomes = resolve_cohort_genomes(all_bs, resolver)
    resolver.save_cache()

    base = Path(f"data/raw/{a.cohort}_extval_{a.drug}")
    own_runs = base / "amrfinder_runs"
    gcache = base / "refseq"
    reuse_glob = f"data/raw/{a.cohort}_*/amrfinder_runs"
    predict = real_predictor(a.drug, own_runs, gcache, reuse_glob,
                             amrfinder_organism=a.amrfinder_organism,
                             registry_organism=a.registry_organism)

    # Fail-fast smoke: one strain through the live path before the full loop.
    if not a.skip_smoke:
        smoke = smoke_predict(genomes["free"], predict)
        print(f"SMOKE: {smoke['gca']} -> {smoke['prediction']} (ok={smoke['ok']})")
        if not smoke["ok"] and not a.allow_degraded:
            print(f"SMOKE FAILED (pass --skip-smoke or --allow-degraded to bypass): {smoke['reason']}")
            return 3

    strict_records, strict_excl = predict_records(genomes["free"], strict_labels, predict)
    strict_conf = conf_from_records(strict_records, strict_excl)
    relaxed_records, relaxed_excl = predict_records(genomes["free"], relaxed_labels, predict)
    relaxed_conf = conf_from_records(relaxed_records, relaxed_excl)

    # Persist per-strain records — the substrate for Step 5's inline clonality recompute.
    base.mkdir(parents=True, exist_ok=True)
    (base / "predictions_strict.json").write_text(json.dumps(strict_records, indent=2), encoding="utf-8")
    (base / "predictions_relaxed.json").write_text(json.dumps(relaxed_records, indent=2), encoding="utf-8")

    # Powering gate — fail-open guard (n_scored==0 / underpowered / high-indeterminate).
    n_attempted_free = len(strict_records)
    n_indeterminate = sum(1 for r in strict_records if str(r["prediction"]).upper() not in ("R", "S"))
    pg = powering_gate(strict_conf, n_attempted_free, n_indeterminate, min_per_class=a.min_per_class)

    leakage_control = (
        f"BioSample-level preflight: {('PASS' if preflight and preflight.get('verdict') == 'PASS' else 'DEGRADED')}; "
        f"{genomes['n_assembly_required']} reads-only BioSamples excluded (ASSEMBLY-REQUIRED)"
    )
    artifact = build_artifact(a.cohort, a.drug, strict=strict_conf, relaxed=relaxed_conf,
                              buckets=buckets, leakage_control=leakage_control,
                              degraded=not (preflight and preflight.get("verdict") == "PASS"),
                              powering=pg, run_degraded=pg["degraded"], run_id=a.run_id,
                              amrfinder_organism=a.amrfinder_organism,
                              registry_organism=a.registry_organism)

    base.mkdir(parents=True, exist_ok=True)
    (base / "selected_strict.tsv").write_text(
        "".join(f"{a}\t{rs}\n" for a, rs in sorted(strict_labels.items())), encoding="utf-8")
    (base / "selected_relaxed.tsv").write_text(
        "".join(f"{a}\t{rs}\n" for a, rs in sorted(relaxed_labels.items())), encoding="utf-8")

    out = Path(f"wiki/external_validation_{a.cohort}_{a.drug}_{a.run_id}_{_date.today().isoformat()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"RESULT strict acc={strict_conf['acc']} sens={strict_conf['sens']} spec={strict_conf['spec']} "
          f"(n={strict_conf['n_scored']}, scored {pg['scored_R']}R/{pg['scored_S']}S, "
          f"indeterminate {pg['indeterminate_fraction']}); artifact: {out}")
    # Fail-open guard: a hard-fail / un-overridden degraded run must NOT exit 0.
    if pg["hard_fail"]:
        print(f"POWERING HARD FAIL: {pg['reasons']}")
        return 3
    if pg["degraded"] and not a.allow_degraded:
        print(f"RUN DEGRADED (pass --allow-degraded to accept): {pg['reasons']}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
