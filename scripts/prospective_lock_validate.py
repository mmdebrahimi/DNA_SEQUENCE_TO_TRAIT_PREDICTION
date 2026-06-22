"""Prospective-lock scorer for the frozen AMR decoder (LA-3, 2026-06-22).

Scores the FROZEN decoder on isolates that became public STRICTLY AFTER the prospective-lock date — a
leakage-free test BY CONSTRUCTION (the decoder could not have been tuned to data that did not yet exist).
Mirrors `scripts/external_cohort_revalidate.py` (verify-then-predict-then-`_conf`); the ONLY new gate is the
TEMPORAL eligibility filter (`prospective_lock.is_prospective_eligible`) instead of accession/BioSample
overlap.

Two-layer honesty:
  - `verify_lock` MUST pass before any scoring — the decoder scoring the post-lock data must be byte-identical
    to the locked decoder (else the "prospective" claim is void). HARD FAIL (exit 2) on lock drift.
  - The cohort is ACCRUING: post-lock measured-AST genomes appear over time. `n_scored == 0` is NOT an error
    — it is the honest "nothing has accrued yet" state (exit 0, status ACCRUING). A real number requires
    enough post-lock isolates per class (powering gate mirrors MIN_PER_CLASS=10).

The pure parts (partition / conf / artifact) are importable + offline-tested (`tests/test_prospective_lock.py`);
the LIVE network+Docker fetch (download genome → AMRFinder → call_resistance) is a manual accruing step
behind a lazy import, exactly like the external-cohort arm.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval.prospective_lock import (  # noqa: E402
    LOCK_SCHEMA, is_prospective_eligible, verify_lock,
)

MIN_PER_CLASS = 10   # mirrors independent_cohort_validate / external_cohort_revalidate powering gate


def partition_by_eligibility(cohort: list[dict], lock_date: str) -> tuple[list[dict], list[dict]]:
    """Split cohort rows into prospective-eligible vs excluded (with a per-row reason).

    Each row needs `first_public_date` (INSDC). Eligible = provably post-lock; everything else is excluded
    fail-closed (undatable / on-or-before lock). Returns (eligible, excluded)."""
    eligible, excluded = [], []
    for row in cohort:
        v = is_prospective_eligible(row.get("first_public_date"), lock_date)
        if v.eligible:
            eligible.append(row)
        else:
            excluded.append({**row, "exclusion_reason": v.reason})
    return eligible, excluded


def conf_from_records(records: list[dict]) -> dict:
    """`independent_cohort_validate._conf` over per-isolate {prediction, y} records."""
    from scripts.independent_cohort_validate import _conf
    return _conf([(r["prediction"], r["y"]) for r in records])


def powering_gate(conf: dict, *, min_per_class: int = MIN_PER_CLASS) -> dict:
    """ACCRUING (n_scored==0) / UNDERPOWERED (a class < min) / POWERED. Never a hard fail — a prospective
    cohort legitimately starts empty and fills over time."""
    n = conf.get("n_scored", 0)
    scored_R = conf.get("tp", 0) + conf.get("fn", 0)
    scored_S = conf.get("tn", 0) + conf.get("fp", 0)
    if n == 0:
        return {"status": "ACCRUING", "n_scored": 0, "detail": "no post-lock isolates have accrued yet"}
    if scored_R < min_per_class or scored_S < min_per_class:
        return {"status": "UNDERPOWERED", "n_scored": n, "scored_R": scored_R, "scored_S": scored_S,
                "detail": f"need >= {min_per_class}/class; have R={scored_R} S={scored_S}"}
    return {"status": "POWERED", "n_scored": n, "scored_R": scored_R, "scored_S": scored_S}


def build_artifact(manifest: dict, organism: str, drug: str, *, lock_ok: bool,
                   eligible: list[dict], excluded: list[dict], conf: dict, powering: dict,
                   generated: str) -> dict:
    """The prospective result, STAMPED with the lock manifest hashes + date so it is provably prospective."""
    return {
        "artifact": "prospective_lock_validation",
        "schema": "prospective-lock-result-v1",
        "generated": generated,
        "organism": organism, "drug": drug,
        "prospective_lock_verified": lock_ok,
        "lock_manifest": {
            "schema": manifest.get("schema"),
            "lock_date": manifest.get("lock_date"),
            "frozen_commit": manifest.get("frozen_commit"),
            "surface_sha256": manifest.get("surface_sha256"),
        },
        "n_eligible_post_lock": len(eligible),
        "n_excluded": len(excluded),
        "exclusion_reasons": _count_reasons(excluded),
        "confusion": conf,
        "powering": powering,
        "honest_note": ("leakage-free BY CONSTRUCTION (every scored isolate became public after lock_date); "
                        "this is a temporal prospective stress test of the FROZEN decoder, accruing over "
                        "time. Not lineage-independent clinical validation."),
    }


def _count_reasons(excluded: list[dict]) -> dict:
    out: dict[str, int] = {}
    for r in excluded:
        out[r.get("exclusion_reason", "?")] = out.get(r.get("exclusion_reason", "?"), 0) + 1
    return out


def _load_manifest(path: Path) -> dict:
    m = json.loads(path.read_text(encoding="utf-8"))
    if m.get("schema") != LOCK_SCHEMA:
        raise ValueError(f"not a {LOCK_SCHEMA} manifest: {path}")
    return m


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path,
                    default=REPO / "wiki" / "prospective_lock_manifest_2026-06-22.json")
    ap.add_argument("--cohort-tsv", type=Path, default=None,
                    help="post-lock cohort TSV: columns biosample, first_public_date, gca, drug, label "
                         "(measured R/S). Omit to run a lock-only verification (no scoring).")
    ap.add_argument("--organism", default="Escherichia_coli_Shigella")
    ap.add_argument("--amrfinder-organism", default="Escherichia")
    ap.add_argument("--drug", default="ciprofloxacin")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    if not args.manifest.exists():
        print(f"ERROR: lock manifest not found at {args.manifest}", file=sys.stderr)
        return 2
    manifest = _load_manifest(args.manifest)

    # HARD GATE: the decoder must be byte-identical to the locked one, or the prospective claim is void.
    v = verify_lock(manifest)
    if not v.ok:
        print(f"ERROR: LOCK BROKEN — frozen decoder surface drifted/missing. "
              f"drifted={v.drifted} missing={v.missing}. A prospective score is only valid against the "
              f"locked decoder; refusing.", file=sys.stderr)
        return 2
    print(f"[prospective-lock] verify_lock OK — decoder byte-identical to lock {manifest['lock_date']} "
          f"({len(manifest['surface_sha256'])} files pinned)")

    today = _date.today().isoformat()
    if args.cohort_tsv is None:
        print("[prospective-lock] lock-only verification (no --cohort-tsv) — lock is valid and ready; "
              "the prospective cohort accrues as post-lock measured-AST genomes appear.")
        return 0

    if not args.cohort_tsv.exists():
        print(f"ERROR: cohort TSV not found at {args.cohort_tsv}", file=sys.stderr)
        return 2
    cohort = _read_cohort(args.cohort_tsv)
    eligible, excluded = partition_by_eligibility(cohort, manifest["lock_date"])
    print(f"[prospective-lock] {len(eligible)} eligible (post-lock) / {len(excluded)} excluded "
          f"({_count_reasons(excluded)})")

    # LIVE scoring (lazy — network + Docker), mirroring external_cohort_revalidate.real_predictor.
    records = _predict_eligible(eligible, args)
    conf = conf_from_records(records) if records else {"n_scored": 0}
    powering = powering_gate(conf)
    art = build_artifact(manifest, args.organism, args.drug, lock_ok=True, eligible=eligible,
                         excluded=excluded, conf=conf, powering=powering, generated=today)
    out = args.out or (REPO / "wiki" / f"prospective_lock_validation_{args.organism}_{args.drug}_{today}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"[prospective-lock] status={powering['status']} -> wrote {out}")
    return 0


def _read_cohort(path: Path) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            vals = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, vals)))
    return rows


def _predict_eligible(eligible: list[dict], args) -> list[dict]:
    """Frozen call_resistance over each eligible isolate (lazy network+Docker import; accruing step)."""
    if not eligible:
        return []
    from dna_decode.eval.amr_rules import call_resistance
    from scripts.organism_drug_validate import _run_dir, ensure_run
    own_runs = REPO / "data" / "amrfinder_runs"
    gcache = REPO / "data" / "refseq_cache"
    records = []
    for row in eligible:
        gca = row.get("gca", "")
        if not gca:
            continue
        ensure_run(gca, own_runs, gcache, args.amrfinder_organism, "*")
        mt = _run_dir(own_runs, gca)
        pred = call_resistance(mt / "main.tsv", args.drug, organism=args.organism)["prediction"]
        records.append({"gca": gca, "prediction": pred, "y": row.get("label", "")})
    return records


if __name__ == "__main__":
    raise SystemExit(main())
